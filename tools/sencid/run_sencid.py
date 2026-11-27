"""Run SenCID, unchanged, on the raw counts of one study; write scores.tsv and run_info.json."""

import argparse
import gzip
import json
import os
import platform
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd
import pkg_resources

TOOL = "SenCID"
SID_NUMBERS = [1, 2, 3, 4, 5, 6]
# The score of the model that SenCID's own recommendation model picks for the sample.
MAIN_SCORE_NAME = "RecSID_Score"
CALL_OF_BINARIZATION = {1: "senescent", 0: "not senescent"}
# SenCID names its output files after its input file.
TOOL_INPUT_NAME = "counts_by_symbol"
REPORTED_PACKAGES = [
    "SenCID", "numpy", "pandas", "scipy", "scikit-learn", "joblib", "scanpy", "anndata",
    "dca", "tensorflow", "keras", "h5py",
]


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", required=True,
                        help="tab-separated, gzip; column GeneID, then one column of raw counts per sample")
    parser.add_argument("--annotation", required=True,
                        help="GEO's Human.GRCh38.p13.annot.tsv.gz (columns GeneID and Symbol are used)")
    parser.add_argument("--out-dir", required=True, help="new folder; refused if it exists")
    return parser.parse_args()


def read_counts(counts_path):
    """Return raw counts, genes (GeneID) by samples."""
    with gzip.open(counts_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    if header[0] != "GeneID":
        raise ValueError(f"first column of {counts_path} is {header[0]!r}, expected 'GeneID'")
    if len(set(header)) != len(header):
        raise ValueError(f"{counts_path} has repeated column names")
    counts = pd.read_csv(counts_path, sep="\t", dtype={"GeneID": str}).set_index("GeneID")
    if not counts.index.is_unique:
        raise ValueError(f"{counts_path} has repeated GeneID values")
    if counts.isna().any().any() or (counts < 0).any().any() or (counts % 1 != 0).any().any():
        raise ValueError(f"{counts_path} holds values that are not whole numbers of zero or more")
    return counts.astype("int64")


def convert_gene_ids_to_symbols(counts, annotation_path):
    """Return the counts with the annotation's Symbol as row name; rows without a Symbol are left out."""
    annotation = pd.read_csv(annotation_path, sep="\t", usecols=["GeneID", "Symbol"], dtype=str)
    symbol_of_gene_id = annotation.dropna(subset=["Symbol"]).set_index("GeneID")["Symbol"]
    symbols = counts.index.to_series().map(symbol_of_gene_id)
    counts_by_symbol = counts[symbols.notna().to_numpy()]
    counts_by_symbol.index = pd.Index(symbols.dropna(), name="Symbol")
    return counts_by_symbol


def count_tool_genes(counts_by_symbol):
    """Compare the gene symbols of the input with SenCID's gene list and with the genes of each model."""
    from SenCID.Pred import GetFeatures

    gene_list_path = pkg_resources.resource_filename("SenCID", "resource/seneset.txt")
    tool_genes = pd.read_csv(gene_list_path, header=None)[0]
    is_found = tool_genes.isin(counts_by_symbol.index)
    missing = set(tool_genes[~is_found])
    total_of_found = counts_by_symbol[counts_by_symbol.index.isin(tool_genes)].sum(axis=1)
    models = {}
    for sid_number in SID_NUMBERS:
        features = GetFeatures(sid_number)
        models[f"SID{sid_number}"] = {
            "n_genes": len(features),
            "genes_missing": sorted(missing.intersection(features)),
        }
    return {
        "n_tool_genes": len(tool_genes),
        "n_tool_genes_found": int(is_found.sum()),
        "n_tool_genes_missing": int((~is_found).sum()),
        "tool_genes_missing": sorted(missing),
        "tool_genes_found_with_zero_counts_in_all_samples": sorted(total_of_found[total_of_found == 0].index),
        "genes_of_each_model": models,
    }


def run_sencid_command(counts_by_symbol, out_dir):
    """Run SenCID's command line on the counts; return the folder with its output files, the command and the seconds."""
    input_dir = os.path.join(out_dir, "tool_input")
    # kopt and Keras, imported by SenCID, read and write settings in the home folder.
    home_dir = os.path.join(out_dir, "tool_home")
    tool_output_dir = os.path.join(out_dir, "sencid_output")
    os.makedirs(input_dir)
    os.makedirs(home_dir)
    counts_by_symbol.to_csv(os.path.join(input_dir, TOOL_INPUT_NAME + ".txt"), sep="\t")
    command = [
        os.path.join(os.path.dirname(sys.executable), "SenCID"),
        "--filepath", input_dir,
        "--fileclass", "txt",
        "--denoising", "f",
        "--binarize", "t",
        "--output_dir", tool_output_dir,
    ]
    start = time.time()
    with open(os.path.join(out_dir, "sencid_log.txt"), "w") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, cwd=input_dir,
                       env=dict(os.environ, HOME=home_dir), check=True)
    seconds = time.time() - start
    shutil.rmtree(input_dir)
    shutil.rmtree(home_dir)
    return tool_output_dir, command, seconds


def read_sencid_output(tool_output_dir, samples):
    """Return SenCID's prediction table of each model and its recommendation table, rows in the order of samples."""
    paths = {f"SID{n}": os.path.join(tool_output_dir, f"predictions_SID{n}_{TOOL_INPUT_NAME}.txt")
             for n in SID_NUMBERS}
    paths["recommendation"] = os.path.join(tool_output_dir, f"RecommendSID_Index_{TOOL_INPUT_NAME}.txt")
    tables = {}
    for name, path in paths.items():
        # SenCID's command line prints an error and exits with code 0 when a file fails.
        if not os.path.exists(path):
            raise RuntimeError(f"SenCID wrote no {os.path.basename(path)}; see sencid_log.txt")
        table = pd.read_csv(path, sep="\t", index_col=0)
        if list(table.index) != samples:
            raise RuntimeError(f"the rows of {path} are not the samples of the input")
        tables[name] = table
    recommendation = tables.pop("recommendation")
    return tables, recommendation


def build_scores(predictions, recommendation):
    """Return the long table sample, tool, score_name, score, call."""
    scores = pd.DataFrame(index=recommendation.index)
    calls = pd.DataFrame(index=recommendation.index)
    for name, prediction in predictions.items():
        scores[f"{name}_Score"] = prediction["SID_Score"]
        calls[f"{name}_Score"] = prediction["Binarization"].map(CALL_OF_BINARIZATION)
        scores[f"{name}_Decision"] = prediction["Decision"]
    for name in predictions:
        scores[f"rec_{name}"] = recommendation[f"rec_{name}"]
    scores["RecSID"] = recommendation["RecSID"].map({f"SID{n}": n for n in SID_NUMBERS})
    if scores["RecSID"].isna().any():
        raise RuntimeError("SenCID's RecSID holds a value other than SID1 to SID6")
    rows = np.arange(len(scores))
    recommended_column = scores["RecSID"].to_numpy(dtype=int) - 1
    score_names = [f"{name}_Score" for name in predictions]
    scores[MAIN_SCORE_NAME] = scores[score_names].to_numpy()[rows, recommended_column]
    calls[MAIN_SCORE_NAME] = calls[score_names].to_numpy()[rows, recommended_column]
    if not np.isfinite(scores.to_numpy(dtype=float)).all() or calls.isna().any().any():
        raise RuntimeError("SenCID returned a score that is not finite or a call that is not 0 or 1")

    long_table = scores.rename_axis(index="sample", columns="score_name").stack().rename("score").reset_index()
    call_of = calls.rename_axis(index="sample", columns="score_name").stack().rename("call").reset_index()
    long_table = long_table.merge(call_of, how="left", on=["sample", "score_name"], validate="one_to_one")
    long_table.insert(1, "tool", TOOL)
    return long_table.fillna({"call": ""})


def main():
    arguments = parse_arguments()
    out_dir = os.path.abspath(arguments.out_dir)
    if os.path.exists(out_dir):
        sys.exit(f"error: --out-dir {out_dir} exists")
    counts = read_counts(arguments.counts)
    counts_by_symbol = convert_gene_ids_to_symbols(counts, arguments.annotation)
    gene_counts = count_tool_genes(counts_by_symbol)
    os.makedirs(out_dir)

    tool_output_dir, command, seconds = run_sencid_command(counts_by_symbol, out_dir)
    predictions, recommendation = read_sencid_output(tool_output_dir, list(counts.columns))
    scores = build_scores(predictions, recommendation)
    scores.to_csv(os.path.join(out_dir, "scores.tsv"), sep="\t", index=False)

    with open(os.path.join(sys.prefix, "sencid_source.json")) as handle:
        source = json.load(handle)
    run_info = {
        "tool": TOOL,
        "tool_repository": source["repository"],
        "tool_commit": source["commit"],
        "tool_command": command,
        "main_score_name": MAIN_SCORE_NAME,
        "counts": os.path.abspath(arguments.counts),
        "annotation": os.path.abspath(arguments.annotation),
        "n_samples": counts.shape[1],
        "n_input_genes": counts.shape[0],
        "n_input_genes_without_symbol": counts.shape[0] - counts_by_symbol.shape[0],
        "n_symbols_on_more_than_one_input_gene": int(counts_by_symbol.index.duplicated().sum()),
        **gene_counts,
        "tool_seconds": round(seconds, 1),
        "python": platform.python_version(),
        "package_versions": {name: pkg_resources.get_distribution(name).version for name in REPORTED_PACKAGES},
    }
    with open(os.path.join(out_dir, "run_info.json"), "w") as handle:
        json.dump(run_info, handle, indent=2)
        handle.write("\n")
    print(f"{TOOL}: {counts.shape[1]} samples, {gene_counts['n_tool_genes_found']} of "
          f"{gene_counts['n_tool_genes']} genes found, {seconds:.0f} s, output in {out_dir}")


if __name__ == "__main__":
    main()
