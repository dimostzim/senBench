#!/usr/bin/env python3
"""Apply the released SenCat RNA marker weights to the raw counts of one study, with the tool's own scripts."""

import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
from importlib.metadata import version
from pathlib import Path

import anndata
import numpy as np
import pandas as pd

TOOL = "sencat"
TOOLS_TSV_NAME = "SenCat ML score"
SCORE_NAME = "rna_score"
MARKERS_FILE = "analysis/transcriptomics_loose5000f_tuned_common_features.csv"
NORMALIZE_SCRIPT = "workflow/scripts/data/normalize_counts.py"
SCORE_SCRIPT = "workflow/scripts/cls/marker_classifier.py"
# The design the tool's configuration gives for its external dataset (config/config.yml, NORMALIZATION_DESIGN).
DESIGN = "sample"
SCORES_FILE = "scores.tsv"
RUN_INFO_FILE = "run_info.json"
TOOL_RESULTS_FILE = "sencat_results.csv"
SCORE_COLUMNS = ("sample", "tool", "score_name", "score", "call")
REPORTED_PACKAGES = ("pydeseq2", "scanpy", "anndata", "pandas", "numpy", "scipy", "formulaic")


def read_raw_counts(counts_path):
    """Return the table of genes (NCBI gene id) by samples, after checking that it holds raw counts."""
    counts = pd.read_csv(counts_path, sep="\t", index_col="GeneID", dtype={"GeneID": str})
    if counts.index.has_duplicates:
        raise ValueError(f"{counts_path}: a GeneID occurs more than once")
    if counts.columns.has_duplicates:
        raise ValueError(f"{counts_path}: a sample occurs more than once")
    values = counts.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any() or (values % 1 != 0).any():
        raise ValueError(f"{counts_path}: the values are not all whole numbers of zero or more")
    return counts.astype("int64")


def convert_to_ensembl_ids(counts, annotation_path):
    """Return the counts by Ensembl gene id: genes without one are dropped, genes that share one are summed."""
    annotation = pd.read_csv(annotation_path, sep="\t", usecols=["GeneID", "EnsemblGeneID"], dtype=str)
    ensembl_ids = annotation.set_index("GeneID")["EnsemblGeneID"]
    unknown = counts.index.difference(ensembl_ids.index)
    if len(unknown):
        raise ValueError(f"{len(unknown)} GeneIDs of the counts are not in {annotation_path}, first {unknown[0]}")
    gene_ensembl_ids = ensembl_ids.reindex(counts.index)
    has_ensembl_id = gene_ensembl_ids.notna()
    by_ensembl_id = counts[has_ensembl_id].groupby(gene_ensembl_ids[has_ensembl_id]).sum()
    conversion = {
        "genes_in_counts": len(counts),
        "genes_without_ensembl_id_dropped": int((~has_ensembl_id).sum()),
        "ensembl_ids_after_conversion": len(by_ensembl_id),
        "ensembl_ids_summed_from_several_genes": int((gene_ensembl_ids[has_ensembl_id].value_counts() > 1).sum()),
    }
    return by_ensembl_id, conversion


def write_counts_h5ad(by_ensembl_id, h5ad_path):
    """Write samples by genes with the sample name as the column that the design names."""
    adata = anndata.AnnData(by_ensembl_id.T)
    adata.obs[DESIGN] = adata.obs_names
    adata.write_h5ad(h5ad_path)


def run_tool_script(repository, script, arguments, log_path):
    subprocess.run(
        [sys.executable, str(repository / script), *arguments, "--log", str(log_path), "--log-level", "INFO"],
        check=True,
    )


def read_tool_commit(repository):
    return subprocess.run(
        ["git", "-C", str(repository), "rev-parse", "HEAD"], check=True, capture_output=True, text=True
    ).stdout.strip()


def count_markers(markers_path, input_genes, normalized_genes):
    """Count the tool's markers that enter the score and the two reasons for which the others are missing."""
    markers = pd.read_csv(markers_path, index_col=0)
    in_input = markers.index.isin(input_genes)
    in_score = markers.index.isin(normalized_genes)
    absolute_weights = markers["coef"].abs()
    return {
        "markers_total": len(markers),
        "markers_found": int(in_score.sum()),
        "markers_missing": int((~in_score).sum()),
        "markers_missing_without_ensembl_id_in_annotation": int((~in_input).sum()),
        "markers_missing_removed_by_tool_gene_filter": int((in_input & ~in_score).sum()),
        "share_of_absolute_weight_found": float(absolute_weights[in_score].sum() / absolute_weights.sum()),
    }


def read_tool_scores(results_path, samples):
    """Return the tool's score per sample, after checking one finite score for every input sample."""
    results = pd.read_csv(results_path, index_col=0)
    if list(results.index) != list(samples):
        raise ValueError(f"{results_path}: the samples differ from the samples of the counts")
    if not np.isfinite(results["score"]).all():
        raise ValueError(f"{results_path}: a score is not a finite number")
    return results["score"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True, help="must not exist")
    parser.add_argument("--repository", type=Path, required=True, help="unchanged clone of wf-ml-markers-senescence")
    arguments = parser.parse_args()

    counts = read_raw_counts(arguments.counts)
    by_ensembl_id, conversion = convert_to_ensembl_ids(counts, arguments.annotation)
    repository = arguments.repository.resolve()
    markers_path = repository / MARKERS_FILE
    out_dir = arguments.out_dir.resolve()
    out_dir.mkdir(parents=True)
    tool_results_path = out_dir / TOOL_RESULTS_FILE

    with tempfile.TemporaryDirectory(dir=out_dir) as work_dir:
        counts_h5ad = Path(work_dir) / "counts.h5ad"
        normalized_h5ad = Path(work_dir) / "normalized.h5ad"
        write_counts_h5ad(by_ensembl_id, counts_h5ad)
        run_tool_script(
            repository,
            NORMALIZE_SCRIPT,
            ["--input-h5ad", str(counts_h5ad), "--design", DESIGN, "--output-h5ad", str(normalized_h5ad)],
            out_dir / "normalize_counts.log",
        )
        run_tool_script(
            repository,
            SCORE_SCRIPT,
            ["--markers", str(markers_path), "--input-h5ad", str(normalized_h5ad),
             "--output-results-csv", str(tool_results_path)],
            out_dir / "marker_classifier.log",
        )
        normalized_genes = anndata.read_h5ad(normalized_h5ad, backed="r").var_names

    scores = read_tool_scores(tool_results_path, counts.columns)
    pd.DataFrame(
        {"sample": scores.index, "tool": TOOL, "score_name": SCORE_NAME, "score": scores.to_numpy(), "call": ""},
        columns=SCORE_COLUMNS,
    ).to_csv(out_dir / SCORES_FILE, sep="\t", index=False)

    run_info = {
        "tool": TOOL,
        "tools_tsv_name": TOOLS_TSV_NAME,
        "tool_repository": str(repository),
        "tool_commit": read_tool_commit(repository),
        "markers_file": MARKERS_FILE,
        "markers_file_md5": hashlib.md5(markers_path.read_bytes()).hexdigest(),
        "design": DESIGN,
        "samples": len(scores),
        **conversion,
        "genes_after_tool_gene_filter": len(normalized_genes),
        **count_markers(markers_path, by_ensembl_id.index, normalized_genes),
        "package_versions": {
            "python": sys.version.split()[0],
            **{package: version(package) for package in REPORTED_PACKAGES},
        },
    }
    (out_dir / RUN_INFO_FILE).write_text(json.dumps(run_info, indent=2) + "\n")


if __name__ == "__main__":
    main()
