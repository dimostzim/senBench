import gzip
import os
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd

TOOL = "SenCID"
SID_NUMBERS = [1, 2, 3, 4, 5, 6]
MAIN_SCORE_NAME = "RecSID_Score"
CALL_OF_BINARIZATION = {1: "senescent", 0: "not senescent"}
TOOL_INPUT_NAME = "counts_by_symbol"


def read_counts(counts_path):
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
    annotation = pd.read_csv(annotation_path, sep="\t", usecols=["GeneID", "Symbol"], dtype=str)
    symbol_of_gene_id = annotation.dropna(subset=["Symbol"]).set_index("GeneID")["Symbol"]
    symbols = counts.index.to_series().map(symbol_of_gene_id)
    counts_by_symbol = counts[symbols.notna().to_numpy()]
    counts_by_symbol.index = pd.Index(symbols.dropna(), name="Symbol")
    return counts_by_symbol


def run_sencid_command(counts_by_symbol, out_dir):
    input_dir = os.path.join(out_dir, "tool_input")
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
    with open(os.path.join(out_dir, "sencid_log.txt"), "w") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, cwd=input_dir,
                       env=dict(os.environ, HOME=home_dir), check=True)
    shutil.rmtree(input_dir)
    shutil.rmtree(home_dir)
    return tool_output_dir


def read_sencid_output(tool_output_dir, samples):
    paths = {f"SID{n}": os.path.join(tool_output_dir, f"predictions_SID{n}_{TOOL_INPUT_NAME}.txt")
             for n in SID_NUMBERS}
    paths["recommendation"] = os.path.join(tool_output_dir, f"RecommendSID_Index_{TOOL_INPUT_NAME}.txt")
    tables = {}
    for name, path in paths.items():
        if not os.path.exists(path):
            raise RuntimeError(f"SenCID wrote no {os.path.basename(path)}; see sencid_log.txt")
        table = pd.read_csv(path, sep="\t", index_col=0)
        if list(table.index) != samples:
            raise RuntimeError(f"the rows of {path} are not the samples of the input")
        tables[name] = table
    recommendation = tables.pop("recommendation")
    return tables, recommendation


def build_scores(predictions, recommendation):
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
    counts_path, annotation_path, out_dir = sys.argv[1:]
    out_dir = os.path.abspath(out_dir)
    if os.path.exists(out_dir):
        sys.exit(f"error: {out_dir} exists")
    counts = read_counts(counts_path)
    counts_by_symbol = convert_gene_ids_to_symbols(counts, annotation_path)
    os.makedirs(out_dir)
    tool_output_dir = run_sencid_command(counts_by_symbol, out_dir)
    predictions, recommendation = read_sencid_output(tool_output_dir, list(counts.columns))
    build_scores(predictions, recommendation).to_csv(os.path.join(out_dir, "scores.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
