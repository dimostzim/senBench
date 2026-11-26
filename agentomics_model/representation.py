"""Frozen abundance-local ordinal representation for bulk RNA-seq."""
from __future__ import annotations

import csv
import json
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import ndtri
from scipy.stats import rankdata
from threadpoolctl import threadpool_limits

LABEL_MAPPING = {"not_senescent": 0, "senescent": 1}
DEFAULT_BATCH_ROWS = 96
SCHEMA = "abundance_local_ordinal_v1"
MIN_REFERENCES = 32
MAX_REFERENCES = 128
SD_FLOOR = 0.5
DEFAULT_BLOCK_TARGETS = 128
DEFAULT_WORKERS = 8
MAX_CONTEXT_GENES = 1024


def _header(path: Path) -> tuple[str, ...]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        header = tuple(next(csv.reader(f), ()))
    if "id" not in header or any(not col for col in header) or len(set(header)) != len(header):
        raise ValueError(f"Missing id, empty column name or duplicate header in {path}")
    return header


def _check_ids(ids, context: str) -> None:
    if not len(ids) or any(not isinstance(i, str) or not i.strip() for i in ids):
        raise ValueError(f"Missing rows or blank id in {context}")
    if len(set(ids)) != len(ids):
        raise ValueError(f"Duplicate id in {context}")


def _batches(path: Path, genes: tuple[str, ...], batch_rows: int):
    if batch_rows < 1:
        raise ValueError("batch_rows must be positive")
    # No NA or numeric coercion for opaque IDs like 'NA' and '001'.
    with pd.read_csv(path, usecols=["id", *genes], converters={"id": str},
                     dtype={g: np.float64 for g in genes}, keep_default_na=False,
                     na_values=["NA", ""], float_precision="round_trip",
                     chunksize=batch_rows, encoding="utf-8-sig", skip_blank_lines=False) as chunks:
        for chunk in chunks:
            ids = chunk.pop("id").tolist()
            if tuple(chunk.columns) != genes:
                chunk = chunk.reindex(columns=genes)
            counts = chunk.to_numpy(dtype=np.float64, copy=False)
            if np.isinf(counts).any() or (counts < 0).any():
                raise ValueError(f"Infinite or negative count in {path}")
            yield ids, counts


def load_labels(split_dir: str | Path, ids: tuple[str, ...]) -> np.ndarray:
    """Supervised steps only: align numeric labels to the expression row order."""
    path = Path(split_dir) / "labels.csv"
    _check_ids(ids, "expression rows")
    _header(path)
    table = pd.read_csv(path, usecols=["id", "numeric_label"], dtype=str,
                        keep_default_na=False, encoding="utf-8-sig", skip_blank_lines=False)
    keys = table["id"].tolist()
    _check_ids(keys, str(path))
    if set(keys) != set(ids):
        raise ValueError(f"Expression/label ID sets differ in {split_dir}")
    if not table["numeric_label"].isin(["0", "1"]).all():
        raise ValueError(f"Expected numeric_label of 0 or 1 in {path}")
    return table.set_index("id").loc[list(ids), "numeric_label"].to_numpy(dtype=np.int64)


@dataclass(frozen=True)
class Features:
    ids: tuple[str, ...]
    x: np.ndarray
    observed_genes: np.ndarray  # Number of retained target genes with observed raw counts.
    supported_targets: np.ndarray  # Observed targets with at least MIN_REFERENCES observed references.
    observed_reference_comparisons: np.ndarray  # Observed target/reference pairs, even below threshold.
    missing_counts: np.ndarray | None = None
    zero_counts: np.ndarray | None = None
    n_input_genes: int | None = None
    ordinal_support: np.ndarray | None = None  # Explicit support independent of the feature value.

    @property
    def raw_target_coverage(self):
        return self.observed_genes / self.x.shape[1]

    @property
    def contrast_support(self):
        return self.supported_targets / self.x.shape[1]

    @property
    def table_gene_counts(self) -> np.ndarray | None:
        """Full header gene count per row (before any fitted filtering)."""
        if self.n_input_genes is None:
            return None
        return np.full(len(self.ids), self.n_input_genes, dtype=np.int64)

    @property
    def table_missingness(self):
        if self.missing_counts is None:
            return None
        if not self.n_input_genes:
            return np.full(len(self.ids), np.nan)
        return self.missing_counts / self.n_input_genes

    @property
    def observed_zero_fraction(self):
        if self.zero_counts is None:
            return None
        observed = self.n_input_genes - self.missing_counts
        return np.divide(self.zero_counts, observed, out=np.full(len(observed), np.nan), where=observed > 0)


def quality_strata(features: Features) -> np.ndarray:
    """Diagnostic/regularization cell IDs, not model inputs.

    IDs 0–3 are <1% missing, 4–7 are >=1% missing; within each group
    the zero-fraction bins start at 0%, 30%, 40%, and 50%, respectively.
    Rows with no observed genes have undefined zero fraction and get -1.
    """
    total = features.n_input_genes
    missing, zeros = features.missing_counts, features.zero_counts
    if total is None or missing is None or zeros is None:
        raise ValueError("Quality strata require full-table diagnostics")
    observed = total - missing
    cell = (100 * missing >= total).astype(np.int64) * 4
    cell += 10 * zeros >= 3 * observed
    cell += 10 * zeros >= 4 * observed
    cell += 2 * zeros >= observed
    cell[observed == 0] = -1
    return cell


def _standardize(values, means, scales):
    x = (values - means) / scales
    x[np.isnan(x)] = 0.
    return x.astype(np.float32)


def _table_counts(counts):
    return np.isnan(counts).sum(axis=1), (counts == 0).sum(axis=1)


def _quantiles(values):
    return np.quantile(values, [0, .25, .5, .75, 1]).tolist() if np.size(values) else []

def _align(counts, source, names):
    indices = np.fromiter((source.get(g, counts.shape[1]) for g in names), dtype=np.intp, count=len(names))
    if (indices == counts.shape[1]).any():
        counts = np.column_stack((counts, np.full(len(counts), np.nan)))
    return counts[:, indices]



@dataclass(frozen=True)
class OrdinalState:
    output_genes: tuple[str, ...]
    reference_genes: tuple[str, ...]
    reference_indices: np.ndarray  # Into reference_genes; -1 denotes padding.
    means: np.ndarray
    scales: np.ndarray


def validate_state(state: OrdinalState) -> None:
    p, q = len(state.output_genes), len(state.reference_genes)
    for names in (state.output_genes, state.reference_genes):
        if len(set(names)) != len(names) or any(not isinstance(g, str) or not g or g == 'id' for g in names):
            raise ValueError('Invalid frozen gene names')
    idx = state.reference_indices
    if (not p or idx.ndim != 2 or idx.shape[0] != p or not MIN_REFERENCES <= idx.shape[1] <= MAX_REFERENCES
            or idx.dtype.kind != 'i' or (idx < -1).any() or (idx >= q).any()
            or ((idx >= 0).sum(axis=1) < MIN_REFERENCES).any()):
        raise ValueError('Invalid frozen reference indices or inadequate reference lists')
    ordered = np.sort(idx, axis=1)
    if ((ordered[:, 1:] == ordered[:, :-1]) & (ordered[:, 1:] >= 0)).any():
        raise ValueError('Duplicate reference in frozen list')
    lookup = {g: j for j, g in enumerate(state.reference_genes)}
    own = np.array([lookup.get(g, -2) for g in state.output_genes])
    if (idx == own[:, None]).any():
        raise ValueError('Self reference in frozen list')
    if (state.means.shape != (p,) or state.scales.shape != (p,) or not np.isfinite(state.means).all()
            or not np.isfinite(state.scales).all() or (state.scales < SD_FLOOR).any()):
        raise ValueError('Invalid frozen normalization statistics')


def select_context_columns(x: np.ndarray, y: np.ndarray, *, max_context: int = MAX_CONTEXT_GENES,
                           block_targets: int = DEFAULT_BLOCK_TARGETS) -> tuple[np.ndarray, np.ndarray]:
    """Return (X indices, B indices), both in frozen ordinal column order.

    Fit on training rows only: choose weak marginal label correlates from the
    actual float32 ordinal matrix, including its neutral-filled entries.
    """
    x, y = np.asarray(x), np.asarray(y)
    if (x.ndim != 2 or not x.shape[0] or not x.shape[1] or y.shape != (len(x),)
            or not np.isin(y, (0, 1)).all() or max_context < 0 or block_targets < 1):
        raise ValueError('Expected a nonempty ordinal matrix, aligned binary labels, and valid limits')
    yc = y.astype(np.float64) - y.mean(dtype=np.float64)
    ynorm = np.linalg.norm(yc)
    scores = np.full(x.shape[1], np.inf, dtype=np.float64)
    for start in range(0, x.shape[1], block_targets):
        sl = slice(start, start + block_targets)
        a = x[:, sl].astype(np.float64, order='F')
        if not np.isfinite(a).all():
            raise ValueError('Expected finite ordinal features')
        varies = np.any(a != a[0], axis=0)
        a -= a.mean(axis=0)
        norm = np.sqrt(np.einsum('ij,ij->j', a, a))
        corr = np.zeros(a.shape[1], dtype=np.float64)
        if ynorm:
            np.divide(np.abs(np.einsum('i,ij->j', yc, a)), ynorm * norm,
                      where=varies, out=corr)
        scores[sl][varies] = np.minimum(corr[varies], 1.0)
    eligible = np.flatnonzero(np.isfinite(scores))
    n_context = min(max_context, max(0, len(eligible) - 1))
    b = np.sort(eligible[np.argsort(scores[eligible], kind='stable')[:n_context]])
    keep = np.ones(x.shape[1], dtype=bool)
    keep[b] = False
    return np.flatnonzero(keep), b


def project_state(state: OrdinalState, x_indices: np.ndarray) -> OrdinalState:
    """Produce an ephemeral X-only view; do not delete raw reference genes.

    Use with iter_features/transform for scoring and iter_coverage for an audit.
    The quadratic artifact stores this projected view, not the training-only B columns.
    """
    x_indices = np.asarray(x_indices)
    if (x_indices.ndim != 1 or not x_indices.size or x_indices.dtype.kind not in 'iu'
            or x_indices[0] < 0 or x_indices[-1] >= len(state.output_genes)
            or np.any(x_indices[1:] <= x_indices[:-1])):
        raise ValueError('Expected strictly increasing X indices within the frozen ordinal columns')
    return OrdinalState(tuple(state.output_genes[i] for i in x_indices), state.reference_genes,
                        state.reference_indices[x_indices], state.means[x_indices], state.scales[x_indices])


def _midrank_percentiles(raw):
    n = (~np.isnan(raw)).sum(axis=1, keepdims=True)
    return (rankdata(raw, axis=1, method='average', nan_policy='omit') - 0.5) / np.maximum(n, 1)


def _absolute_correlation(percentiles, y, candidates, block_targets):
    """Pairwise-complete Pearson; zero for constant rank or constant label."""
    score = np.empty(len(candidates), dtype=np.float64)
    labels = y[:, None]
    for start in range(0, len(candidates), block_targets):
        stop = min(start + block_targets, len(candidates))
        rank = percentiles[:, candidates[start:stop]]
        obs = ~np.isnan(rank)
        n = obs.sum(axis=0)
        rank_mean = np.nansum(rank, axis=0) / n
        label_mean = np.sum(obs * labels, axis=0) / n
        a = np.where(obs, rank - rank_mean, 0.)
        b = np.where(obs, labels - label_mean, 0.)
        denom = np.sqrt(np.sum(a * a, axis=0) * np.sum(b * b, axis=0))
        varies = np.max(rank, axis=0, where=obs, initial=-np.inf) > np.min(rank, axis=0, where=obs, initial=np.inf)
        score[start:stop] = np.divide(np.abs(np.sum(a * b, axis=0)), denom,
                                      out=np.zeros_like(denom), where=varies & (denom > 0))
    return score


def _retain_references(coordinate, candidates, correlation):
    ranked = np.lexsort((candidates, coordinate[candidates]))
    retain = []
    for band in np.array_split(ranked, 10):
        ordered = band[np.lexsort((candidates[band], correlation[band]))]
        retain.extend(candidates[ordered[:(len(band) + 1) // 2]])
    return np.sort(np.asarray(retain, dtype=np.intp))  # Header order.


def _nearest_references(coordinate, retained):
    width = min(len(retained), MAX_REFERENCES)
    idx = np.full((len(coordinate), width), -1, dtype=np.int32)
    ref_coord = coordinate[retained]
    for j, center in enumerate(coordinate):
        distance = np.abs(ref_coord - center)
        self_pos = np.searchsorted(retained, j)
        is_ref = self_pos < len(retained) and retained[self_pos] == j
        if is_ref:
            distance[self_pos] = np.inf
        k = min(width, len(retained) - int(is_ref))
        if k == 0:
            continue
        # Partition to find the cutoff; resolve boundary ties by header order.
        cutoff = np.partition(distance, k - 1)[k - 1]
        closer = np.flatnonzero(distance < cutoff)
        ties = np.flatnonzero(distance == cutoff)[:k - len(closer)]
        subset = np.concatenate((closer, ties))
        idx[j, :k] = subset[np.lexsort((subset, distance[subset]))]
    return idx


def _ordinal_values(targets, refs, indices, block_targets):
    """For a bounded row block, bound the comparison tensor along targets too."""
    n, p = targets.shape
    values = np.full((n, p), np.nan, dtype=np.float64)
    pair_totals = np.zeros(n, dtype=np.int64)
    refs = np.column_stack((refs, np.full(n, np.nan)))  # -1 gathers missing.
    for start in range(0, p, block_targets):
        cols = slice(start, start + block_targets)
        selected = refs[:, indices[cols]]
        target = targets[:, cols]
        observed = ~np.isnan(target)
        m = np.count_nonzero(~np.isnan(selected), axis=2)
        pair_totals += np.sum(m * observed, axis=1)
        lo = np.count_nonzero(selected < target[:, :, None], axis=2)
        eq = np.count_nonzero(selected == target[:, :, None], axis=2)
        supported = observed & (m >= MIN_REFERENCES)
        u = (lo[supported] + 0.5 * eq[supported] + 0.5) / (m[supported] + 1)
        values[:, cols][supported] = ndtri(np.clip(u, 0.02, 0.98))
    return values, pair_totals


def iter_coverage(input_dir: str | Path, state: OrdinalState, batch_rows: int = DEFAULT_BATCH_ROWS,
                  *, block_targets: int = DEFAULT_BLOCK_TARGETS, include_support: bool = False):
    """Single-scan audit without calculating features or scoring. Yields Features-like
    diagnostics via CoverageBatch; counts refer to the final retained output genes.
    """
    validate_state(state)
    if block_targets < 1:
        raise ValueError('block_targets must be positive')
    path = Path(input_dir) / 'data.csv'
    genes = tuple(g for g in _header(path) if g != 'id')
    source = {g: j for j, g in enumerate(genes)}
    seen = set()
    for keys, counts in _batches(path, genes, batch_rows):
        _check_ids(keys, str(path))
        if seen.intersection(keys):
            raise ValueError(f'Duplicate id in {path}')
        seen.update(keys)
        targets = ~np.isnan(_align(counts, source, state.output_genes))
        refs = ~np.isnan(_align(counts, source, state.reference_genes))
        refs = np.column_stack((refs, np.zeros(len(keys), dtype=bool)))
        support = np.zeros(len(keys), dtype=np.int64)
        mask = np.empty_like(targets) if include_support else None
        for start in range(0, targets.shape[1], block_targets):
            cols = slice(start, start + block_targets)
            m = np.count_nonzero(refs[:, state.reference_indices[cols]], axis=2)
            supported = targets[:, cols] & (m >= MIN_REFERENCES)
            support += np.count_nonzero(supported, axis=1)
            if mask is not None:
                mask[:, cols] = supported
        missing, zeros = _table_counts(counts)
        yield CoverageBatch(tuple(keys), targets.sum(axis=1), support, missing, zeros, len(genes), mask)
    if not seen:
        raise ValueError(f'Missing rows in {path}')


@dataclass(frozen=True)
class CoverageBatch:
    ids: tuple[str, ...]
    observed_genes: np.ndarray
    supported_targets: np.ndarray
    missing_counts: np.ndarray
    zero_counts: np.ndarray
    n_input_genes: int
    ordinal_support: np.ndarray | None = None


def fit_transform(input_dir: str | Path, *, split_dir: str | Path, batch_rows: int = DEFAULT_BATCH_ROWS,
                  block_targets: int = DEFAULT_BLOCK_TARGETS, workers: int = DEFAULT_WORKERS,
                  report: bool = True) -> tuple[OrdinalState, Features]:
    """Fit from raw training counts and aligned training labels, then emit train features."""
    if min(batch_rows, block_targets, workers) < 1:
        raise ValueError('batch_rows, block_targets and workers must be positive')
    path = Path(input_dir) / 'data.csv'
    genes = tuple(g for g in _header(path) if g != 'id')
    ids, blocks, missing, zeros = [], [], [], []
    observed = np.zeros(len(genes), dtype=np.int64)
    positive = np.zeros_like(observed)
    for keys, counts in _batches(path, genes, batch_rows):
        ids.extend(keys)
        blocks.append(counts)
        m, z = _table_counts(counts)
        missing.append(m)
        zeros.append(z)
        observed += (~np.isnan(counts)).sum(axis=0)
        positive += (counts > 0).sum(axis=0)
    ids = tuple(ids)
    y = load_labels(split_dir, ids)
    keep = (1000 * observed >= 997 * len(ids)) & (5 * positive >= 2 * observed)
    target_names = np.asarray(genes)[keep]
    candidates = np.flatnonzero((5 * positive >= 3 * observed)[keep])
    if len(candidates) < MIN_REFERENCES:
        raise ValueError('Too few eligible reference candidates')
    raw = np.concatenate([a[:, keep] for a in blocks], axis=0)
    del blocks, counts

    row_blocks = [slice(i, i + batch_rows) for i in range(0, len(raw), batch_rows)]
    with ThreadPoolExecutor(max_workers=min(workers, len(row_blocks))) as pool:
        ranks = np.empty_like(raw)

        def rank_block(rows):
            ranks[rows] = _midrank_percentiles(raw[rows])

        list(pool.map(rank_block, row_blocks))
        coordinate = np.nanmedian(ranks, axis=0)
        correlation = _absolute_correlation(ranks, y, candidates, block_targets)
        del ranks
        retained = _retain_references(coordinate, candidates, correlation)
        if len(retained) < MIN_REFERENCES:
            raise ValueError('Too few references after abundance-band filtering')
        indices = _nearest_references(coordinate, retained)
        values = np.empty_like(raw)

        def ordinal_block(rows):
            values[rows], _ = _ordinal_values(raw[rows], raw[rows][:, retained], indices, block_targets)

        list(pool.map(ordinal_block, row_blocks))

    available = ~np.isnan(values)
    lo = np.min(values, axis=0, where=available, initial=np.inf)
    hi = np.max(values, axis=0, where=available, initial=-np.inf)
    varies = hi > lo
    if not varies.any():
        raise ValueError('No target has nonconstant available training ordinal values')
    no_support = int((~available.any(axis=0)).sum())
    del available
    final_raw, values, indices = raw[:, varies], values[:, varies], indices[varies]
    used = np.unique(indices[indices >= 0])
    indices = np.where(indices >= 0, np.searchsorted(used, indices), -1).astype(np.int32)
    means, sd = np.nanmean(values, axis=0), np.nanstd(values, axis=0, ddof=0)
    state = OrdinalState(tuple(target_names[varies].tolist()), tuple(target_names[retained[used]].tolist()),
                         indices, means, np.maximum(sd, SD_FLOOR))
    validate_state(state)
    observed_targets = (~np.isnan(final_raw)).sum(axis=1)
    support = ~np.isnan(values)
    supported_targets = support.sum(axis=1)
    # Training diagnostics need not repeat comparisons. Count reference observations using masks.
    pair_totals = np.zeros(len(ids), dtype=np.int64)
    if np.isnan(raw[:, retained[used]]).any():
        for rows in row_blocks:
            refs_observed = np.column_stack((~np.isnan(raw[rows][:, retained[used]]),
                                             np.zeros(len(raw[rows]), dtype=bool)))
            for start in range(0, len(state.output_genes), block_targets):
                cols = slice(start, start + block_targets)
                m = np.count_nonzero(refs_observed[:, indices[cols]], axis=2)
                pair_totals[rows] += np.sum(m * (~np.isnan(final_raw[rows, cols])), axis=1)
    else:
        pair_totals = (~np.isnan(final_raw)) @ (indices >= 0).sum(axis=1)
    train = Features(ids, _standardize(values, state.means, state.scales), observed_targets,
                     supported_targets, pair_totals, np.concatenate(missing), np.concatenate(zeros), len(genes), support)
    if report:
        print(json.dumps({'ordinal_fit': {
            'rows': len(ids), 'input_genes': len(genes), 'target_candidates': len(target_names),
            'reference_candidates': len(candidates), 'selected_references': len(retained),
            'used_references': len(used), 'output_features': len(state.output_genes),
            'dropped_no_support': no_support, 'dropped_constant': int((~varies).sum()) - no_support,
            'sd_floor_count': int((sd < SD_FLOOR).sum()),
            'observed_targets_quantiles': _quantiles(observed_targets),
            'supported_targets_quantiles': _quantiles(supported_targets),
            'zero_support_rows': int((supported_targets == 0).sum()),
        }}), flush=True)
    return state, train


def iter_features(input_dir: str | Path, state: OrdinalState, batch_rows: int = DEFAULT_BATCH_ROWS,
                  *, diagnostics: bool = False, block_targets: int = DEFAULT_BLOCK_TARGETS):
    """Replay frozen comparisons in bounded batches, with no labels or batch statistics."""
    validate_state(state)
    if block_targets < 1:
        raise ValueError('block_targets must be positive')
    path = Path(input_dir) / 'data.csv'
    needed = set(state.output_genes) | set(state.reference_genes)
    genes = tuple(g for g in _header(path) if g != 'id' and (diagnostics or g in needed))
    source = {g: i for i, g in enumerate(genes)}
    seen = set()
    for keys, counts in _batches(path, genes, batch_rows):
        if not keys:
            continue
        _check_ids(keys, str(path))
        if not seen.isdisjoint(keys):
            raise ValueError(f'Duplicate id in {path}')
        seen.update(keys)
        targets = _align(counts, source, state.output_genes)
        refs = _align(counts, source, state.reference_genes)
        values, pair_counts = _ordinal_values(targets, refs, state.reference_indices, block_targets)
        m, z = _table_counts(counts) if diagnostics else (None, None)
        support = ~np.isnan(values)
        yield Features(tuple(keys), _standardize(values, state.means, state.scales),
                       (~np.isnan(targets)).sum(axis=1), support.sum(axis=1),
                       pair_counts, m, z, len(genes) if diagnostics else None, support)


def transform(input_dir: str | Path, state: OrdinalState, batch_rows: int = DEFAULT_BATCH_ROWS,
              *, diagnostics: bool = False, block_targets: int = DEFAULT_BLOCK_TARGETS) -> Features:
    """Materialized interface for compatibility with the training evaluator."""
    batches = list(iter_features(input_dir, state, batch_rows, diagnostics=diagnostics,
                                 block_targets=block_targets))
    if not batches:
        raise ValueError(f'Missing rows in {input_dir}')
    def cat(name):
        return np.concatenate([getattr(batch, name) for batch in batches])
    return Features(tuple(i for batch in batches for i in batch.ids), cat('x'),
                    cat('observed_genes'), cat('supported_targets'), cat('observed_reference_comparisons'),
                    cat('missing_counts') if diagnostics else None,
                    cat('zero_counts') if diagnostics else None, batches[0].n_input_genes,
                    cat('ordinal_support'))
PANEL_LIMIT = 128


def select_stable_panel(x, y, cells, *, blas_threads=8, block_columns=256):
    """Rank on quality-cell stability, then greedily prune pooled Pearson redundancy."""
    x, y, cells = np.asarray(x), np.asarray(y), np.asarray(cells)
    if (x.ndim != 2 or y.shape != (len(x),) or cells.shape != y.shape
            or not np.isin(y, [0, 1]).all() or not (np.any(y == 0) and np.any(y == 1))
            or not np.isin(cells, np.arange(-1, 8)).all()):
        raise ValueError("Invalid aligned binary labels, features or quality cells")
    if block_columns < 1 or blas_threads < 1:
        raise ValueError("Block size and thread count must be positive")
    counts = np.zeros((9, 2), dtype=np.int64)
    np.add.at(counts, (cells.astype(np.intp) + 1, y.astype(np.intp)), 1)
    eligible = np.flatnonzero((counts[1:] >= 20).all(axis=1))
    fallback = len(eligible) < 2
    groups = [(y == 0, y == 1)] if fallback else [((cells == g) & (y == 0), (cells == g) & (y == 1)) for g in eligible]
    p = x.shape[1]
    scores, nonconstant = np.empty(p), np.empty(p, dtype=bool)
    for start in range(0, p, block_columns):
        sl = slice(start, start + block_columns)
        block = np.asarray(x[:, sl], dtype=np.float64)
        if not np.isfinite(block).all():
            raise ValueError("Panel selection requires finite, neutral-encoded features")
        nonconstant[sl] = np.max(block, axis=0) > np.min(block, axis=0)
        mu0 = np.stack([block[control].mean(axis=0) for control, _ in groups])
        mu1 = np.stack([block[case].mean(axis=0) for _, case in groups])
        d = mu1 - mu0
        scores[sl] = (np.abs(d[0]) if fallback else
                      np.abs(d.mean(axis=0)) - 0.5 * d.std(axis=0)
                      - 0.25 * (mu0.std(axis=0) + mu1.std(axis=0)))
    candidates = np.flatnonzero(nonconstant)
    ranking = candidates[np.argsort(-scores[candidates], kind="stable")]
    selected = []
    rejected = 0
    # Cache only accepted unit vectors and the current rank block. No p-by-p matrix.
    units = np.empty((min(PANEL_LIMIT, len(ranking)), len(x)), dtype=np.float64)
    with threadpool_limits(limits=blas_threads, user_api="blas"):
        for start in range(0, len(ranking), block_columns):
            cols = ranking[start:start + block_columns]
            z = np.ascontiguousarray(x[:, cols].T, dtype=np.float64)
            z -= z.mean(axis=1, keepdims=True)
            z /= np.linalg.norm(z, axis=1, keepdims=True)
            within = np.abs(z @ z.T)
            blocked = np.zeros(len(cols), dtype=bool)
            if selected:
                blocked[:] = np.any(np.abs(z @ units[:len(selected)].T) >= .85, axis=1)
            for i, col in enumerate(cols):
                if blocked[i]:
                    rejected += 1
                    continue
                units[len(selected)] = z[i]
                selected.append(int(col))
                if len(selected) == PANEL_LIMIT:
                    break
                blocked[i + 1:] |= within[i, i + 1:] >= .85
            if len(selected) == PANEL_LIMIT:
                break
    indices = np.asarray(selected, dtype=np.intp)
    audit = dict(quality_cell_counts={str(g): counts[g + 1].tolist() for g in range(-1, 8)},
                 quality_cell_count_order=[0, 1], eligible_quality_cells=eligible.tolist(),
                 pooled_fallback=fallback, constant_columns=int((~nonconstant).sum()),
                 score_quantiles=_quantiles(scores[candidates]),
                 selected_score_quantiles=_quantiles(scores[indices]),
                 redundancy_rejections=rejected, ranked_candidates_examined=len(indices) + rejected,
                 panel_size=len(indices))
    return indices, audit



QUADRATIC_SCHEMA = 'abundance_local_ordinal_quadratic_v1'
PAIR_WEIGHT = 0.35
PAIR_BLOCK = 256


@dataclass(frozen=True)
class QuadraticState:
    ordinal: OrdinalState  # Already projected to X, retaining required raw references.
    selected_indices: np.ndarray  # Indices in X, in frozen column order.
    pair_indices: np.ndarray  # (k, 2) direct indices in X, lexicographically ordered.
    product_means: np.ndarray
    product_scales: np.ndarray

    @property
    def n_features(self):
        return len(self.ordinal.output_genes) + len(self.pair_indices)


def validate_quadratic_state(state: QuadraticState) -> None:
    validate_state(state.ordinal)
    p = len(state.ordinal.output_genes)
    selected, pairs = state.selected_indices, state.pair_indices
    if (selected.ndim != 1 or selected.dtype.kind not in 'iu' or len(selected) > PANEL_LIMIT
            or (selected < 0).any() or (selected >= p).any()
            or np.any(selected[1:] <= selected[:-1])):
        raise ValueError('Invalid frozen panel indices')
    if (pairs.ndim != 2 or pairs.shape[1] != 2 or pairs.dtype.kind not in 'iu'
            or not np.isin(pairs, selected).all() or np.any(pairs[:, 0] >= pairs[:, 1])):
        raise ValueError('Invalid frozen interaction indices')
    key = pairs[:, 0].astype(np.int64) * p + pairs[:, 1]
    if np.any(key[1:] <= key[:-1]):
        raise ValueError('Interaction indices must be unique and in frozen order')
    mean, scale = state.product_means, state.product_scales
    if (mean.shape != (len(pairs),) or scale.shape != mean.shape
            or not np.isfinite(mean).all() or not np.isfinite(scale).all() or (scale < 1).any()):
        raise ValueError('Invalid frozen product moments')


def _product_inputs(x, support):
    if x.ndim != 2 or support.shape != x.shape or support.dtype != np.bool_:
        raise ValueError('Expected ordinal X and its aligned Boolean support matrix')


def _products(x, support, pairs):
    i, j = pairs.T
    supported = support[:, i] & support[:, j]
    products = np.clip(x[:, i].astype(np.float64), -3., 3.)
    products *= np.clip(x[:, j].astype(np.float64), -3., 3.)
    return products, supported


def fit_products(x, support, selected):
    """Return retained X-index pairs and supported-only float64 moments."""
    _product_inputs(x, support)
    selected = np.sort(np.asarray(selected, dtype=np.intp))
    a, b = np.triu_indices(len(selected), k=1)
    candidates = np.column_stack((selected[a], selected[b]))
    retained, means, scales = [], [], []
    for start in range(0, len(candidates), PAIR_BLOCK):
        pairs = candidates[start:start + PAIR_BLOCK]
        prod, mask = _products(x, support, pairs)
        varies = (np.max(prod, axis=0, where=mask, initial=-np.inf)
                  > np.min(prod, axis=0, where=mask, initial=np.inf))
        prod, mask = prod[:, varies], mask[:, varies]
        n = mask.sum(axis=0)
        mean = np.sum(prod, axis=0, where=mask) / n
        delta = np.where(mask, prod - mean, 0.)
        sd = np.sqrt(np.sum(delta * delta, axis=0) / n)
        retained.append(pairs[varies])
        means.append(mean)
        scales.append(np.maximum(sd, 1.))
    if not retained:
        return candidates, np.empty(0), np.empty(0)
    return np.concatenate(retained), np.concatenate(means), np.concatenate(scales)


def build_design(x, support, state: QuadraticState):
    """Return float64 Z and per-row usable pair counts; diagnostics stay outside Z."""
    _product_inputs(x, support)
    p = x.shape[1]
    if p != len(state.ordinal.output_genes):
        raise ValueError('Ordinal X width does not match projected state')
    z = np.empty((len(x), state.n_features), dtype=np.float64)
    z[:, :p] = x
    usable = np.zeros(len(x), dtype=np.int64)
    for start in range(0, len(state.pair_indices), PAIR_BLOCK):
        sl = slice(start, start + PAIR_BLOCK)
        product, mask = _products(x, support, state.pair_indices[sl])
        z[:, p + start:p + start + product.shape[1]] = np.where(
            mask, PAIR_WEIGHT * ((product - state.product_means[sl]) / state.product_scales[sl]), 0.)
        usable += mask.sum(axis=1)
    return z, usable


@dataclass(frozen=True)
class QuadraticTraining:
    state: QuadraticState
    ids: tuple[str, ...]
    y: np.ndarray
    cells: np.ndarray
    z: np.ndarray
    b: np.ndarray
    audit: dict


def fit_quadratic_design(input_dir: str | Path, *, split_dir: str | Path,
                         batch_rows: int = DEFAULT_BATCH_ROWS,
                         workers: int = DEFAULT_WORKERS, report: bool = True) -> QuadraticTraining:
    """Fit the partition and interactions once on aligned training rows only."""
    ordinal, train = fit_transform(input_dir, split_dir=split_dir, batch_rows=batch_rows,
                                   workers=workers, report=report)
    y = load_labels(split_dir, train.ids)
    cells = quality_strata(train)
    x_idx, b_idx = select_context_columns(train.x, y)
    x, support = train.x[:, x_idx], train.ordinal_support[:, x_idx]
    selected, audit = select_stable_panel(x, y, cells)
    selected.sort()
    pairs, means, scales = fit_products(x, support, selected)
    state = QuadraticState(project_state(ordinal, x_idx), selected, pairs, means, scales)
    validate_quadratic_state(state)
    z, usable = build_design(x, support, state)
    audit.update(linear_features=x.shape[1], context_features=len(b_idx),
                 candidate_pairs=len(selected) * (len(selected) - 1) // 2,
                 retained_pairs=len(pairs), supported_ordinal_count_quantiles=_quantiles(support.sum(axis=1)),
                 usable_pair_count_quantiles=_quantiles(usable), no_usable_pair_rows=int((usable == 0).sum()),
                 usable_pair_count_quantiles_by_quality_cell={
                     str(g): _quantiles(usable[cells == g]) for g in range(-1, 8)})
    if report:
        print(json.dumps({'quadratic_fit': audit}), flush=True)
    return QuadraticTraining(state, train.ids, y, cells, z, train.x[:, b_idx], audit)


def iter_quadratic_features(input_dir: str | Path, state: QuadraticState,
                            batch_rows: int = DEFAULT_BATCH_ROWS, *, diagnostics: bool = False):
    """Yield (ordinal batch, Z, usable pair counts), without labels or fitting."""
    validate_quadratic_state(state)
    any_rows = False
    for batch in iter_features(input_dir, state.ordinal, batch_rows, diagnostics=diagnostics):
        any_rows = True
        z, usable = build_design(batch.x, batch.ordinal_support, state)
        yield batch, z, usable
    if not any_rows:
        raise ValueError(f'Missing rows in {input_dir}')


def iter_quadratic_coverage(input_dir: str | Path, state: QuadraticState,
                            batch_rows: int = DEFAULT_BATCH_ROWS):
    """Yield (table and ordinal coverage diagnostics, usable pair counts) without forming Z."""
    validate_quadratic_state(state)
    for batch in iter_coverage(input_dir, state.ordinal, batch_rows, include_support=True):
        usable = np.zeros(len(batch.ids), dtype=np.int64)
        for start in range(0, len(state.pair_indices), PAIR_BLOCK):
            i, j = state.pair_indices[start:start + PAIR_BLOCK].T
            usable += np.count_nonzero(batch.ordinal_support[:, i] & batch.ordinal_support[:, j], axis=1)
        yield batch, usable
