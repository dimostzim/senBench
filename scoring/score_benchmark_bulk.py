#!/usr/bin/env python3
"""Collect the tools' scores on the bulk benchmark and compute a provisional AUROC per tool, score and dataset."""

import argparse
import csv
from pathlib import Path

import pandas as pd

SENESCENT = "senescent"
NOT_SENESCENT = "not senescent"
PROVISIONAL = ("provisional: no split exists yet, the overlap with the tools' training data is known only as far as "
               "the lists go, and nothing was tuned")


def _read(path, **options):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, **options)


def compute_auroc(scores, senescent):
    """Return the probability that a senescent sample has a higher score than a not senescent one (ties count half)."""
    ranks = scores.rank()
    positives = int(senescent.sum())
    negatives = len(scores) - positives
    return (ranks[senescent].sum() - positives * (positives + 1) / 2) / (positives * negatives)


def summarise(group):
    """Return the AUROC and the counts of one tool, score and dataset; the calls are counted where the tool makes them."""
    senescent = group.state == SENESCENT
    row = {"senescent": int(senescent.sum()), "not_senescent": int((~senescent).sum()),
           "auroc": round(compute_auroc(group.score.astype(float), senescent), 4)}
    called = group.call != ""
    if called.any() and not called.all():
        raise SystemExit(f"{group.name}: a call is missing for some samples")
    # Empty where the tool makes no call.
    for name, call, state in (("true_positives", SENESCENT, senescent), ("false_positives", SENESCENT, ~senescent),
                              ("false_negatives", NOT_SENESCENT, senescent),
                              ("true_negatives", NOT_SENESCENT, ~senescent)):
        row[name] = int(((group.call == call) & state).sum()) if called.any() else pd.NA
    return pd.Series(row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tool-runs-dir", type=Path, required=True, help="output folder of run_tools_on_benchmark.py")
    parser.add_argument("--units", type=Path, required=True, help="units.tsv of that run")
    parser.add_argument("--samples", type=Path, required=True, help="benchmark_bulk_samples.tsv")
    parser.add_argument("--main-scores", type=Path, required=True)
    parser.add_argument("--tool-overlap", type=Path, required=True, help="output of mark_tool_training_overlap.py")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.out_dir.exists():
        raise SystemExit(f"{args.out_dir} exists")

    samples = _read(args.samples)
    main_scores = _read(args.main_scores)
    overlap = _read(args.tool_overlap).set_index("dataset")
    with args.units.open(newline="") as handle:
        units = pd.DataFrame(list(csv.reader(handle, delimiter="\t")), columns=["unit", "status", "note"])
    units[["wrapper", "dataset"]] = units.unit.str.split("/", expand=True)

    frames = []
    for unit in units[units.status == "ok"].itertuples():
        scores = _read(args.tool_runs_dir / unit.wrapper / unit.dataset / "scores.tsv")
        frames.append(scores.assign(wrapper=unit.wrapper, dataset=unit.dataset))
    scores = pd.concat(frames).merge(samples[["dataset", "sample", "study", "state"]], on=["dataset", "sample"],
                                     validate="m:1")
    expected = units[units.status == "ok"].merge(samples.groupby("dataset").size().rename("samples").reset_index())
    per_run = scores.groupby(["wrapper", "dataset", "score_name"]).size().rename("rows").reset_index().merge(expected)
    if not (per_run.rows == per_run.samples).all():
        raise SystemExit(f"scores do not cover every sample:\n{per_run[per_run.rows != per_run.samples]}")

    aurocs = scores.groupby(["wrapper", "tool", "score_name", "dataset", "study"]).apply(
        summarise, include_groups=False).reset_index()
    counts = ["senescent", "not_senescent", "true_positives", "false_positives", "false_negatives", "true_negatives"]
    aurocs[counts] = aurocs[counts].astype("Int64")
    aurocs = aurocs.merge(main_scores[["wrapper", "tools_tsv_name", "main_score"]], on="wrapper", validate="m:1")
    aurocs["main_score"] = (aurocs.score_name == aurocs.main_score).map({True: "yes", False: "no"})

    def training_use(row):
        if row.dataset not in overlap.index:
            return "not looked up"
        if not row.tools_tsv_name:
            return "no list for this tool"
        for column, use in (("tools_training", "training or model selection"), ("tools_test", "test by the authors"),
                            ("tools_applied_only", "applied by the authors")):
            if row.tools_tsv_name in overlap.at[row.dataset, column].split(";"):
                return use
        return "in no list of this tool"

    aurocs["dataset_in_the_tool_lists"] = aurocs.apply(training_use, axis=1)
    aurocs["dataset_in_no_list_of_any_tool"] = aurocs.dataset.map(
        lambda dataset: overlap.at[dataset, "in_no_list"] if dataset in overlap.index else "not looked up")
    aurocs["status"] = PROVISIONAL
    failed = units[units.status != "ok"][["wrapper", "dataset", "status", "note"]]

    args.out_dir.mkdir(parents=True)
    scores.to_csv(args.out_dir / "tool_scores.tsv.gz", sep="\t", index=False)
    aurocs.drop(columns="tools_tsv_name").to_csv(args.out_dir / "provisional_auroc.tsv", sep="\t", index=False)
    failed.to_csv(args.out_dir / "failed_tool_runs.tsv", sep="\t", index=False)
    main_rows = aurocs[aurocs.main_score == "yes"]
    print(main_rows.groupby("wrapper").auroc.describe().round(3).to_string())
    print(f"{len(failed)} failed tool runs")


if __name__ == "__main__":
    main()
