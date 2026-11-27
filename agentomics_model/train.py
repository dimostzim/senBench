"""Fit one penalized joint ordinal–quadratic logistic model and evaluate it frozen."""
import argparse
import csv
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize
from scipy.special import expit
from sklearn.metrics import average_precision_score, confusion_matrix, matthews_corrcoef, roc_auc_score
from sklearn.utils.extmath import randomized_svd
from threadpoolctl import threadpool_limits

from training_reporter import TrainingReporter
from pipeline import LinearClassifier, MODEL_NAME, load_pipeline, save_pipeline
from representation import fit_quadratic_design, load_labels, iter_quadratic_features
from inference import _scorer, COLUMNS
from matching import audit_graph, build_graph, test_graph

MIN_CELL = 20
STRATUM_WEIGHT = 0.20
CONTEXT_RANK_CAP = 384
CONTEXT_FULL_WEIGHT_RANK = 96
CONTEXT_AUDIT_BANDS = ((0, 96), (96, 192), (192, 384))
NONLINEAR_WEIGHT = 0.35
NONLINEAR_RANK_CAP = 96
NONLINEAR_FULL_WEIGHT_RANK = 24
NONLINEAR_AUDIT_BANDS = ((0, 24), (24, 48), (48, 96))


def penalty_design(x, y, strata):
    """Return centered means and weights for equally averaged within-class penalties."""
    masks = [[(y == c) & (strata == g) for g in range(-1, 8)] for c in (0, 1)]
    counts = [[int(m.sum()) for m in row] for row in masks]
    eligible = [[g for g in range(8) if counts[c][g + 1] >= MIN_CELL] for c in (0, 1)]
    active = [g if len(g) >= 2 else [] for g in eligible]
    blocks, weights = [], []
    k = sum(bool(cells) for cells in active)
    for c, cells in enumerate(active):
        if cells:
            means = np.stack([x[masks[c][g + 1]].mean(axis=0, dtype=np.float64) for g in cells])
            blocks.append(means - means.mean(axis=0))
            weights.extend([1.0 / (k * len(cells))] * len(cells))
    d = np.concatenate(blocks, axis=0) if blocks else np.empty((0, x.shape[1]), dtype=np.float64)
    audit = {'cells': list(range(-1, 8)), 'class_counts': counts,
             'eligible_cells_by_class': eligible, 'active_cells_by_class': active}
    return d, np.asarray(weights, dtype=np.float64), audit


def context_scale(rank):
    """Square roots of the one-based context weights."""
    return np.minimum(1.0, CONTEXT_FULL_WEIGHT_RANK / np.arange(1, rank + 1, dtype=np.float64))


def requested_context_rank(n_controls, n_context):
    return max(0, min(CONTEXT_RANK_CAP, n_controls - 1, n_context))


def nonlinear_scale(rank):
    return np.minimum(1.0, NONLINEAR_FULL_WEIGHT_RANK / np.arange(1, rank + 1, dtype=np.float64))


def control_basis(c):
    """Iteration 17's retained left singular vectors, in descending order."""
    rank = requested_context_rank(*c.shape)
    if rank < 1:
        return np.empty((len(c), 0)), np.empty(0)
    u, s, _ = randomized_svd(c, n_components=rank, n_iter=7, random_state=2026)
    keep = s > (1e-8 * s[0])
    return u[:, keep], s[keep]


def nonlinear_basis(c, u):
    """Leading directions of the residualized, centered quadratic kernel."""
    m = len(c)
    empty = (np.empty((m, 0)), np.empty(0))
    if m < 2:
        return empty
    s = np.einsum('ij,ij->', c, c) / m
    if s == 0:
        return empty
    k = (c @ c.T) / s
    np.square(k, out=k)
    mean = k.mean(axis=0)
    k -= mean[None, :]
    k -= mean[:, None]
    k += mean.mean()
    k = (k + k.T) * 0.5
    top = eigh(k, subset_by_index=(m - 1, m - 1), eigvals_only=True)[0]
    if top <= 0:
        return empty
    if u.shape[1]:
        k -= u @ (u.T @ k)
        k -= (k @ u) @ u.T
    k = (k + k.T) * 0.5
    spectrum, v = eigh(k, subset_by_index=(max(0, m - NONLINEAR_RANK_CAP), m - 1))
    spectrum, v = spectrum[::-1], v[:, ::-1]
    keep = spectrum > 1e-8 * top
    return v[:, keep], spectrum[keep]


def context_design(x, y, b, include_nonlinear=False):
    controls = y == 0
    n = int(controls.sum())
    c = b[controls].astype(np.float64)
    xc = x[controls].astype(np.float64, copy=False)
    if n:
        c -= c.mean(axis=0)
        xc -= xc.mean(axis=0)
    u, spectrum = control_basis(c)
    scale = np.sqrt(max(n, 1))
    d = (u.T @ xc) / scale
    d *= context_scale(len(spectrum))[:, None]
    if not include_nonlinear:
        return d, spectrum
    v, eig = nonlinear_basis(c, u)
    e = (v.T @ xc) / scale
    e *= nonlinear_scale(len(eig))[:, None]
    overlap = float(np.abs(u.T @ v).max()) if u.shape[1] and v.shape[1] else None
    return d, spectrum, e, eig, overlap


class Objective:
    def __init__(self, x, y, d, a, context=None, nonlinear=None, gpu=False, graph=None):
        self.x, self.y, self.d, self.a = x, y, d, a
        self.graph = graph
        self.binary_weight = 0.75 if graph else 1.0
        self.context = np.empty((0, x.shape[1])) if context is None else context
        self.nonlinear = np.empty((0, x.shape[1])) if nonlinear is None else nonlinear
        self.sign = 1.0 - 2.0 * y
        self.device_x = None
        if gpu:
            try:
                import torch
            except ImportError:
                pass
            else:
                if torch.cuda.is_available():
                    self.torch = torch
                    try:
                        self.device_x = torch.from_numpy(x).to('cuda')
                    except torch.cuda.OutOfMemoryError:
                        pass

    def margins(self, w):
        if self.device_x is None:
            return self.x @ w
        return self.device_x.mv(self.torch.from_numpy(w).to('cuda')).cpu().numpy()

    def backproject(self, residual):
        if self.device_x is None:
            return self.x.T @ residual
        return self.device_x.T.mv(self.torch.from_numpy(residual).to('cuda')).cpu().numpy()

    def terms(self, theta):
        w, b = theta[:-1], theta[-1]
        margin = self.margins(w) + b
        loss = np.mean(np.logaddexp(0.0, self.sign * margin))
        proj = self.d @ w
        return margin, proj, loss, w @ w, self.a @ (proj * proj)

    def pair_terms(self, scores):
        g = self.graph
        if not g:
            return 0.0, 0.0
        delta = scores[g.control] - scores[g.positive]
        force = 0.25 * g.weight * expit(delta)
        n = len(self.y)
        residual = (np.bincount(g.control, weights=force, minlength=n)
                    - np.bincount(g.positive, weights=force, minlength=n))
        return g.weight @ np.logaddexp(0.0, delta), residual

    def __call__(self, theta):
        margin, proj, loss, ridge, strata = self.terms(theta)
        w = theta[:-1]
        residual = self.binary_weight * self.sign * expit(self.sign * margin) / len(self.y)
        pair_loss, pair_residual = self.pair_terms(margin)
        gw = self.backproject(residual + pair_residual) + 2.0 * w + (2.0 * STRATUM_WEIGHT) * (self.d.T @ (self.a * proj))
        context_proj = self.context @ w
        gw += self.context.T @ context_proj
        nonlinear_proj = self.nonlinear @ w
        gw += (2.0 * NONLINEAR_WEIGHT) * (self.nonlinear.T @ nonlinear_proj)
        total = (self.binary_weight * loss + 0.25 * pair_loss + ridge + STRATUM_WEIGHT * strata
                 + 0.50 * (context_proj @ context_proj)
                 + NONLINEAR_WEIGHT * (nonlinear_proj @ nonlinear_proj))
        return total, np.r_[gw, residual.sum()]


def test_objective(gpu):
    rng = np.random.default_rng(238)
    ordinal = rng.normal(size=(240, 3))
    raw_pairs = np.prod(np.clip(ordinal[:, [[0, 1], [1, 2]]], -3, 3), axis=2)
    pairs = (raw_pairs - raw_pairs.mean(axis=0)) / np.maximum(raw_pairs.std(axis=0), 1.)
    x = np.column_stack((ordinal, 0.35 * pairs))
    y = np.repeat([0., 1.], [180, 60])
    g = np.repeat([0, 4, 5, 0, 6], [60, 60, 60, 20, 40])
    d, a, audit = penalty_design(x, y, g)
    assert audit['active_cells_by_class'] == [[0, 4, 5], [0, 6]]
    np.testing.assert_allclose(a, [1/6] * 3 + [1/4] * 2)
    theta = rng.normal(size=6)
    assert np.all(theta[3:5] != 0)
    b = rng.normal(size=(len(y), 108))
    context, spectrum = context_design(x, y, b)
    assert context.shape == (108, x.shape[1]) and context.dtype == np.float64
    assert np.linalg.norm(context[96:] @ theta[:-1]) > 0
    full_weights = context_scale(384)**2
    np.testing.assert_array_equal(full_weights[:96], np.ones(96))
    np.testing.assert_array_equal(full_weights[[0, 95, 191, 383]], [1.0, 1.0, 0.25, 0.0625])
    weights = full_weights[:len(spectrum)]
    for n_controls, expected_rank in ((0, 0), (1, 0), (5, 4), (410, 384)):
        by = np.r_[np.zeros(n_controls), 1.0]
        bx = rng.normal(size=(len(by), x.shape[1]))
        wide = rng.normal(size=(len(by), 400))
        assert context_design(bx, by, wide)[0].shape == (expected_rank, x.shape[1])
    bc = b[y == 0] - b[y == 0].mean(axis=0)
    u, s, _ = np.linalg.svd(bc, full_matrices=False)
    logits = x[y == 0] @ theta[:-1] + theta[-1]
    expected = (u.T @ (logits - logits.mean()))**2 / len(logits)
    np.testing.assert_allclose((context @ theta[:-1])**2, weights * expected, atol=1e-13, rtol=1e-11)
    np.testing.assert_allclose(spectrum, s, rtol=1e-12)
    for empty in (np.empty((len(y), 0)), np.ones((len(y), 2))):
        assert context_design(x, y, empty)[0].shape == (0, x.shape[1])
    rank_one = np.column_stack([b[:, 0], 1e-10 * b[:, 1]])
    assert context_design(x, y, rank_one)[0].shape == (1, x.shape[1])
    # An even, centered quadratic in a symmetric one-dimensional control fixture
    # must survive removal of the odd linear direction. Its projection has no
    # eigenvalue multiplier.
    q = np.array([-2., -1., 0., 1., 2.])[:, None]
    q_squared = q**2 - (q**2).mean()
    _, _, qe, qeig, q_overlap = context_design(q_squared, np.zeros(len(q)), q, True)
    assert len(qeig) == 1 and q_overlap < 1e-12
    np.testing.assert_allclose(qeig, [np.sum(q_squared**2) / np.mean(q**2)**2], rtol=1e-13)
    np.testing.assert_allclose(np.abs(qe[0, 0]), np.linalg.norm(q_squared) / np.sqrt(len(q)), rtol=1e-13)
    assert nonlinear_basis(q[:1], np.empty((1, 0)))[0].shape == (1, 0)
    assert nonlinear_basis(np.zeros_like(q), np.empty((len(q), 0)))[0].shape == (len(q), 0)
    assert nonlinear_basis(q[[0, 4]], np.empty((2, 0)))[0].shape == (2, 0)

    d_linear, _, e, eig, overlap = context_design(x, y, b, True)
    np.testing.assert_array_equal(d_linear, context)
    assert len(eig) > 24 and overlap < 1e-11
    assert np.linalg.norm(e[24:] @ theta[:-1]) > 1e-4
    v_scale = nonlinear_scale(len(eig))
    np.testing.assert_array_equal(v_scale[[0, 23, 47]], [1., 1., .5])
    v, _ = nonlinear_basis(bc, control_basis(bc)[0])
    np.testing.assert_allclose((e @ theta[:-1])**2,
                               (v.T @ (logits - logits.mean()))**2 / len(logits) * v_scale**2,
                               atol=1e-13, rtol=1e-11)
    # Check the full objective with both taper tails active at once.
    graph, rank = build_graph(b, y)
    assert rank == 24 and len(graph) > 0
    cpu = Objective(x, y, d, a, context, e, graph=graph)
    check = Objective(x, y, d, a, context, e, gpu=gpu, graph=graph)
    unpaired = Objective(x, y, d, a, context, e)
    scores = x @ theta[:-1] + theta[-1]
    pair_loss = sum(w * np.logaddexp(0., scores[c] - scores[p])
                    for p, c, w in zip(graph.positive, graph.control, graph.weight))
    np.testing.assert_allclose(cpu(theta)[0] - unpaired(theta)[0],
                               0.25 * (pair_loss - unpaired.terms(theta)[2]), atol=1e-13)
    np.testing.assert_allclose(cpu(theta)[1][-1], 0.75 * unpaired(theta)[1][-1], atol=1e-14)
    empty, _ = build_graph(np.ones((len(y), 2)), y)
    fallback = Objective(x, y, d, a, context, e, graph=empty)
    np.testing.assert_array_equal(fallback(theta)[1], unpaired(theta)[1])
    assert fallback(theta)[0] == unpaired(theta)[0]
    h = 1e-5 * np.eye(len(theta))
    numeric = np.array([(cpu(theta + step)[0] - cpu(theta - step)[0]) / (2e-5) for step in h])
    np.testing.assert_allclose(check(theta)[1], numeric, atol=1e-9, rtol=1e-8)
    # Compare the isolated nonlinear contribution against an independently calculated value.
    increment = 0.35 * np.sum(v_scale**2 * (v.T @ (logits - logits.mean()))**2) / len(logits)
    np.testing.assert_allclose(cpu(theta)[0] - Objective(x, y, d, a, context, graph=graph)(theta)[0], increment,
                               atol=1e-13, rtol=1e-11)
    equal = np.repeat(y[:, None], x.shape[1], axis=1)
    d, a, _ = penalty_design(equal, y, g)
    assert Objective(equal, y, d, a).terms(theta)[-1] == 0.0
    d, a, _ = penalty_design(x, y, np.zeros_like(g))
    assert d.shape == (0, x.shape[1]) and Objective(x, y, d, a).terms(theta)[-1] == 0.0


def evaluate(s, y):
    p = expit(s)
    decision = s >= 0.0
    tn, fp, fn, tp = confusion_matrix(y, decision, labels=[0, 1]).ravel()
    both_classes = np.unique(y).size == 2
    return {'n': len(y), 'auroc': float(roc_auc_score(y, s)) if both_classes else None,
            'ap': float(average_precision_score(y, s)) if both_classes else None,
            'log_loss': float(np.mean(np.logaddexp(0.0, (1.0 - 2.0 * y) * s))),
            'accuracy': float(np.mean(decision == y)),
            'mcc': float(matthews_corrcoef(y, decision)),
            'prevalence': float(y.mean()), 'mean_probability': float(p.mean()),
            'tn': int(tn), 'fp': int(fp), 'fn': int(fn), 'tp': int(tp)}


def main():
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument('--train-data', type=Path, required=True)
    parser.add_argument('--validation-data', type=Path, required=True)
    parser.add_argument('--artifacts-dir', type=Path, required=True)
    args = parser.parse_args()
    tr, va, out = (p.resolve() for p in (args.train_data, args.validation_data, args.artifacts_dir))
    if any(out.is_relative_to(d / 'input') for d in (tr, va)):
        parser.error('--artifacts-dir cannot be inside a split input folder')
    same_split = tr == va
    if not same_split:
        def label_ids(split):
            with (split / 'labels.csv').open(encoding='utf-8-sig', newline='') as f:
                return {r['id'] for r in csv.DictReader(f)}
        if not label_ids(tr).isdisjoint(label_ids(va)):
            raise ValueError('Distinct training/validation splits contain overlapping IDs')
    mode = 'in-sample (same split)' if same_split else 'held-out validation'
    reporter = TrainingReporter()
    reporter.report_unavailable('SciPy L-BFGS-B exposes solver iterations, not epochs or batches; only the converged fit is evaluated.')
    with threadpool_limits(limits=8, user_api='blas'):
        test_graph()
        test_objective(gpu=True)
        train = fit_quadratic_design(tr / 'input', split_dir=tr)
        state, ids = train.state, train.ids
        y = train.y.astype(np.float64)
        if np.unique(y).size != 2:
            raise ValueError('Both classes are required for training')
        x = train.z
        n_linear = len(state.ordinal.output_genes)
        n_controls = int((y == 0).sum())
        requested_rank = requested_context_rank(n_controls, train.b.shape[1])
        context, spectrum, nonlinear, eigenvalues, overlap = context_design(x, y, train.b, True)
        assert context.shape == (len(spectrum), x.shape[1]) and len(spectrum) <= requested_rank
        assert nonlinear.shape == (len(eigenvalues), x.shape[1]) and len(eigenvalues) <= NONLINEAR_RANK_CAP
        assert np.all(spectrum[:-1] >= spectrum[1:])
        assert np.all(eigenvalues[:-1] >= eigenvalues[1:])
        d, a, audit = penalty_design(x, y, train.cells)
        print(json.dumps({'train_quality_strata': audit, 'n_ordinal': n_linear,
                          'n_context': train.b.shape[1], 'n_predictors': state.n_features,
                          'n_controls': n_controls,
                          'requested_control_rank': requested_rank,
                          'retained_control_rank': len(spectrum),
                          'retained_singular_value_range': [float(spectrum[-1]), float(spectrum[0])]
                          if len(spectrum) else []}), flush=True)
        print(json.dumps({'nonlinear_context': {
            'retained_rank': len(eigenvalues),
            'retained_eigenvalue_range': [float(eigenvalues[-1]), float(eigenvalues[0])]
            if len(eigenvalues) else [],
            'max_abs_linear_overlap': overlap}}), flush=True)
        graph, rank = build_graph(train.b, y)
        print(json.dumps({'context_matching': audit_graph(graph, rank, y, train.cells)}), flush=True)
        del train
        objective = Objective(x, y, d, a, context, nonlinear, gpu=True, graph=graph)
        res = minimize(objective, np.zeros(x.shape[1] + 1), method='L-BFGS-B', jac=True,
                       options={'maxiter': 1000, 'gtol': 1e-6, 'ftol': 1e-11})
        s, _, loss, ridge, strata = objective.terms(res.x)
        context_proj = context @ res.x[:-1]
        unweighted_sq = (context_proj / context_scale(len(spectrum)))**2
        rank_sums = {f'{lo + 1}-{hi}': {
                         'unweighted': float(unweighted_sq[lo:hi].sum()),
                         'weighted': float(np.square(context_proj[lo:hi]).sum())}
                     for lo, hi in CONTEXT_AUDIT_BANDS}
        q = float(context_proj @ context_proj)
        nonlinear_sq = np.square(nonlinear @ res.x[:-1])
        nonlinear_unweighted_sq = nonlinear_sq / nonlinear_scale(len(eigenvalues))**2
        nonlinear_sums = {f'{lo + 1}-{hi}': {
                              'unweighted': float(nonlinear_unweighted_sq[lo:hi].sum()),
                              'weighted': float(nonlinear_sq[lo:hi].sum())}
                          for lo, hi in NONLINEAR_AUDIT_BANDS}
        q_nonlinear = float(nonlinear_sq.sum())
        block_audit = {}
        for name, sl in [('ordinal', slice(0, n_linear)), ('pair', slice(n_linear, None))]:
            coef = res.x[:-1][sl]
            block_audit[name] = {'coefficient_l2': float(np.linalg.norm(coef)),
                                 'logit_sd': float((x[:, sl] @ coef).std())}
        print(json.dumps({'fit': {'converged': bool(res.success), 'message': str(res.message),
                                 'iterations': res.nit, 'function_evals': res.nfev,
                                 'gradient_inf_norm': float(np.abs(res.jac).max()),
                                 'backend': 'CUDA' if objective.device_x is not None else 'CPU',
                                 'n_train': len(y), 'n_features': x.shape[1], 'log_loss': loss,
                                 'binary_weight': objective.binary_weight,
                                 'pair_log_loss': objective.pair_terms(s)[0] if graph else None,
                                 'ridge': ridge, 'stratum_penalty': strata,
                                 'weighted_stratum_penalty': STRATUM_WEIGHT * strata,
                                 'score_projection_sq_by_rank': rank_sums,
                                 'tapered_projection_sq': q, 'weighted_context_penalty': 0.50 * q,
                                 'nonlinear_projection_sq_by_rank': nonlinear_sums,
                                 'tapered_nonlinear_projection_sq': q_nonlinear,
                                 'weighted_nonlinear_penalty': NONLINEAR_WEIGHT * q_nonlinear,
                                 'blocks': block_audit, 'training_logit_sd': float(s.std()),
                                 'objective': res.fun}}), flush=True)
        if not res.success or not np.isfinite(res.fun):
            raise RuntimeError(f'L-BFGS-B did not converge: {res.message}')
        model = LinearClassifier(res.x[:-1].copy(), float(res.x[-1]))
        del objective
        out.mkdir(parents=True, exist_ok=True)
        path = out / MODEL_NAME
        save_pipeline(path, model, state)
        saved, frozen = load_pipeline(path)
        def check_state(expected, actual):
            for key, value in expected.items():
                if isinstance(value, dict):
                    check_state(value, actual[key])
                else:
                    np.testing.assert_array_equal(value, actual[key])
        check_state(asdict(state), asdict(frozen))
        np.testing.assert_array_equal(saved.w, model.w)
        assert saved.b == model.b
        np.testing.assert_allclose(saved.decision_function(x), s, atol=1e-12, rtol=1e-12)
        empty = np.zeros((1, x.shape[1]))
        np.testing.assert_array_equal(saved.decision_function(empty), [saved.b])
        np.testing.assert_array_equal(saved.predict_proba(empty)[0, 1], expit(saved.b))
        if not same_split:
            del x
        val_ids, score_parts = [], []
        score = _scorer(saved, use_cuda=True)
        for batch, z, _ in iter_quadratic_features(va / 'input', frozen):
            if same_split:
                sl = slice(len(val_ids), len(val_ids) + len(batch.ids))
                assert batch.ids == ids[sl]
                np.testing.assert_array_equal(z, x[sl])
            val_ids.extend(batch.ids)
            score_parts.append(score(z))
        if same_split:
            assert tuple(val_ids) == ids
        elif not set(ids).isdisjoint(val_ids):
            raise ValueError('Distinct training/validation splits contain overlapping IDs')
        y_val = load_labels(va, tuple(val_ids)).astype(np.float64)
        val_scores = np.concatenate(score_parts)
        metrics = evaluate(val_scores, y_val)
        if not same_split:
            prediction_path = Path(__file__).resolve().parent / 'validation_predictions.csv'
            with prediction_path.open('w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerow(COLUMNS)
                writer.writerows(zip(val_ids, (val_scores >= 0).astype(int),
                                     expit(-val_scores), expit(val_scores)))

    print(json.dumps({'evaluation': mode, **metrics}), flush=True)
    reporter.report_epoch(epoch=1, train_loss=float(loss), val_loss=metrics['log_loss'],
                          val_metric_name='AUROC' if metrics['auroc'] is not None else None,
                          val_metric=metrics['auroc'])


if __name__ == '__main__':
    main()
