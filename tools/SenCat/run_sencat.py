import subprocess
import sys
import tempfile
from pathlib import Path

import anndata
import numpy as np
import pandas as pd

TOOL = "sencat"
SCORE_NAME = "rna_score"
MARKERS_FILE = "analysis/transcriptomics_loose5000f_tuned_common_features.csv"
NORMALIZE_SCRIPT = "workflow/scripts/data/normalize_counts.py"
SCORE_SCRIPT = "workflow/scripts/cls/marker_classifier.py"
DESIGN = "sample"
SCORE_COLUMNS = ("sample", "tool", "score_name", "score", "call")


def read_raw_counts(counts_path):
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
    annotation = pd.read_csv(annotation_path, sep="\t", usecols=["GeneID", "EnsemblGeneID"], dtype=str)
    ensembl_ids = annotation.set_index("GeneID")["EnsemblGeneID"]
    unknown = counts.index.difference(ensembl_ids.index)
    if len(unknown):
        raise ValueError(f"{len(unknown)} GeneIDs of the counts are not in {annotation_path}, first {unknown[0]}")
    gene_ensembl_ids = ensembl_ids.reindex(counts.index)
    has_ensembl_id = gene_ensembl_ids.notna()
    return counts[has_ensembl_id].groupby(gene_ensembl_ids[has_ensembl_id]).sum()


def write_counts_h5ad(by_ensembl_id, h5ad_path):
    adata = anndata.AnnData(by_ensembl_id.T)
    adata.obs[DESIGN] = adata.obs_names
    adata.write_h5ad(h5ad_path)


def run_tool_script(repository, script, arguments, log_path):
    subprocess.run(
        [sys.executable, str(repository / script), *arguments, "--log", str(log_path), "--log-level", "INFO"],
        check=True,
    )


def read_tool_scores(results_path, samples):
    results = pd.read_csv(results_path, index_col=0)
    if list(results.index) != list(samples):
        raise ValueError(f"{results_path}: the samples differ from the samples of the counts")
    if not np.isfinite(results["score"]).all():
        raise ValueError(f"{results_path}: a score is not a finite number")
    return results["score"]


def main():
    counts_path, annotation_path, repository, out_dir = (Path(argument) for argument in sys.argv[1:])
    counts = read_raw_counts(counts_path)
    by_ensembl_id = convert_to_ensembl_ids(counts, annotation_path)
    repository = repository.resolve()
    out_dir = out_dir.resolve()
    out_dir.mkdir(parents=True)
    tool_results_path = out_dir / "sencat_results.csv"

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
            ["--markers", str(repository / MARKERS_FILE), "--input-h5ad", str(normalized_h5ad),
             "--output-results-csv", str(tool_results_path)],
            out_dir / "marker_classifier.log",
        )

    scores = read_tool_scores(tool_results_path, counts.columns)
    pd.DataFrame(
        {"sample": scores.index, "tool": TOOL, "score_name": SCORE_NAME, "score": scores.to_numpy(), "call": ""},
        columns=SCORE_COLUMNS,
    ).to_csv(out_dir / "scores.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
