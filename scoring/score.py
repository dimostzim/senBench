#!/usr/bin/env python3
"""Score the predictions of methods on the labelled samples of `test`, and compare every method with a reference.

A prediction file holds one score per sample (columns dataset, sample and score; a higher score means more senescent),
and a method is named after its file. Measures: the AUROC across studies, and the pooled AUROC and pooled AUPRC
(average precision) over all samples in one ranking. Every measure of every method has an interval from the same
resamples of laboratory groups, and the reference minus every other method is computed on those resamples, on the
samples both have scores for.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import across_studies_measures as measures
from benchmark_bulk_classes import NOT_SENESCENT, SENESCENT

AUROC_ACROSS = "AUROC across studies"
POOLED_AUROC = "pooled AUROC"
POOLED_AUPRC = "pooled AUPRC"
MEASURES = (AUROC_ACROSS, POOLED_AUROC, POOLED_AUPRC)
PREDICTION_COLUMNS = ["dataset", "sample", "score"]
TEST_SAMPLES = Path(__file__).resolve().parent / "test_samples.tsv"
RESULT_COLUMNS = ["method", "measure", "value", "interval_low", "interval_high", "random_score", "samples", "senescent",
                  "not_senescent", "datasets", "laboratory_groups", "reference_minus_method", "difference_interval_low",
                  "difference_interval_high", "difference_samples"]


def read_table(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def read_test_samples(path):
    """Return dataset, sample, study, laboratory group and state of every labelled sample of `test`."""
    test_set = read_table(path)
    test_set["laboratory_group"] = test_set.laboratory_group.astype(int)
    return test_set[["dataset", "sample", "study", "laboratory_group", "state"]]


def read_predictions(path):
    """Return dataset, sample and score of a prediction file, which names each sample once with a finite score."""
    predictions = pd.read_csv(path, sep="\t", dtype={"dataset": str, "sample": str, "score": float},
                              keep_default_na=False, float_precision="round_trip")
    if list(predictions.columns) != PREDICTION_COLUMNS:
        raise ValueError(f"{path}: the columns are not {', '.join(PREDICTION_COLUMNS)}")
    if predictions.duplicated(["dataset", "sample"]).any():
        raise ValueError(f"{path}: a sample occurs twice")
    if not np.isfinite(predictions.score).all():
        raise ValueError(f"{path}: a score is not a finite number")
    return predictions


def join_test_labels(predictions, test_set):
    """Return the predictions with the study, laboratory group and state of each sample."""
    rows = predictions.merge(test_set, on=["dataset", "sample"], how="left", validate="1:1", indicator=True)
    if (rows._merge != "both").any():
        raise ValueError("a predicted sample is not a labelled sample of `test`")
    return rows.drop(columns="_merge")


def count_entered(rows):
    return {"laboratory_groups": rows.laboratory_group.nunique(), "studies": rows.study.nunique(),
            "datasets": rows.dataset.nunique(), "samples": len(rows),
            "positives": int((rows.state == SENESCENT).sum()), "negatives": int((rows.state == NOT_SENESCENT).sum())}


def compute_average_precision(scores, positive, weights):
    """Return the average precision of the ranking by score for each row of sample weights; tied scores form one
    threshold, as in scikit-learn's average_precision_score. NaN for a row without positive weight."""
    order = np.argsort(-scores, kind="stable")
    scores, positive, weights = scores[order], positive[order], weights[:, order]
    last_of_tie = np.r_[scores[1:] != scores[:-1], True]
    true_positives = np.cumsum(weights * positive, axis=1)[:, last_of_tie]
    false_positives = np.cumsum(weights * ~positive, axis=1)[:, last_of_tie]
    total = true_positives[:, -1]
    with np.errstate(invalid="ignore", divide="ignore"):
        recall = true_positives / total[:, None]
        precision = true_positives / (true_positives + false_positives)
        steps = np.diff(recall, axis=1, prepend=0)
        value = np.where(steps > 0, steps * precision, 0).sum(axis=1)
    return np.where(total > 0, value, np.nan)


def resample_measures(rows, studies):
    """Return per measure the value on the samples and its value in every resample of laboratory groups."""
    table = measures.tabulate_pairs(*measures.split_states(rows), studies)
    across, across_draws, _ = measures.resample_auroc_across_studies(table)
    study_weights = measures.weigh_studies(table.groups)
    with np.errstate(invalid="ignore", divide="ignore"):
        pooled_draws = (((study_weights @ table.wins) * study_weights).sum(axis=1)
                        / ((study_weights @ table.comparisons) * study_weights).sum(axis=1))
    scores = rows.score.to_numpy(dtype=float)
    positive = (rows.state == SENESCENT).to_numpy()
    precision_draws = compute_average_precision(scores, positive,
                                                measures.weigh_studies(rows.laboratory_group.to_numpy()))
    return {AUROC_ACROSS: (across, across_draws),
            POOLED_AUROC: (float(table.wins.sum() / table.comparisons.sum()), pooled_draws),
            POOLED_AUPRC: (float(compute_average_precision(scores, positive, np.ones((1, len(rows))))[0]),
                           precision_draws)}


def describe_random_value(measure, rows):
    """Return the value of a random score: 0.5 for an AUROC, about the share of senescent samples for the AUPRC."""
    return float((rows.state == SENESCENT).mean()) if measure == POOLED_AUPRC else 0.5


def score_methods(reference_rows, rows_of_method):
    """Return one row per method and measure, with its interval and the number of the samples of `reference_rows`
    (the reference method's, or those of the `test` set) that the method has no score for."""
    results = []
    for method, rows in rows_of_method.items():
        for measure, (value, draws) in resample_measures(rows, measures.list_studies(rows)).items():
            estimate = measures.summarise_resamples(value, draws)
            results.append({"method": method, "measure": measure, "value": estimate.value,
                            "interval_low": estimate.interval_low, "interval_high": estimate.interval_high,
                            "resamples_without_value": estimate.resamples_without_value,
                            "random_score": describe_random_value(measure, rows), **count_entered(rows),
                            "samples_without_score": len(reference_rows) - len(rows)})
    return pd.DataFrame(results)


def compare_methods(reference, rows_of_method):
    """Return per other method and measure the reference's value minus the method's, on the samples both have
    scores for, the same studies and the same resamples."""
    results = []
    for method, rows in rows_of_method.items():
        if method == reference:
            continue
        common = rows_of_method[reference][["dataset", "sample"]].merge(rows[["dataset", "sample"]])
        kept = [of_method.merge(common) for of_method in (rows_of_method[reference], rows)]
        studies = measures.list_studies(kept[0])
        resampled = [resample_measures(of_method, studies) for of_method in kept]
        for measure in resampled[0]:
            (reference_value, reference_draws), (method_value, method_draws) = (
                of_method[measure] for of_method in resampled)
            estimate = measures.summarise_resamples(reference_value - method_value, reference_draws - method_draws)
            excludes_zero = estimate.interval_low > 0 or estimate.interval_high < 0
            results.append({"reference": reference, "method": method, "measure": measure,
                            "value_reference": reference_value, "value_method": method_value,
                            "difference": estimate.value, "interval_low": estimate.interval_low,
                            "interval_high": estimate.interval_high,
                            "resamples_without_value": estimate.resamples_without_value,
                            "interval_excludes_zero": "yes" if excludes_zero else "no", **count_entered(kept[0])})
    return pd.DataFrame(results)


def tabulate_results(scored, differences):
    """Return one row per method and measure, the methods ordered by their AUROC across studies, with the reference
    minus the method (empty for the reference itself)."""
    differences = differences.rename(columns={
        "difference": "reference_minus_method", "interval_low": "difference_interval_low",
        "interval_high": "difference_interval_high", "samples": "difference_samples"})
    table = scored.rename(columns={"positives": "senescent", "negatives": "not_senescent"}).merge(
        differences[["method", "measure", "reference_minus_method", "difference_interval_low",
                     "difference_interval_high", "difference_samples"]], on=["method", "measure"], how="left")
    order = scored[scored.measure == AUROC_ACROSS].sort_values("value", ascending=False, kind="stable").method
    table["method"] = pd.Categorical(table.method, categories=order)
    table["measure"] = pd.Categorical(table.measure, categories=MEASURES)
    table["difference_samples"] = table.difference_samples.astype("Int64")
    return table.sort_values(["method", "measure"])[RESULT_COLUMNS]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--predictions", type=Path, nargs="+", required=True,
                        help="one prediction file per method, named <method>.tsv")
    parser.add_argument("--reference", required=True, help="the method that every other method is compared with")
    parser.add_argument("--test-samples", type=Path, default=TEST_SAMPLES,
                        help="test_samples.tsv: study, dataset, sample, laboratory_group and state of every labelled "
                             "sample of `test`")
    parser.add_argument("--out", type=Path, required=True, help="the results table to write; it must not exist")
    args = parser.parse_args()
    names = [path.stem for path in args.predictions]
    if len(set(names)) != len(names):
        parser.error("two prediction files have the same name")
    if args.reference not in names:
        parser.error(f"no prediction file is named {args.reference}.tsv")
    if args.out.exists():
        raise SystemExit(f"{args.out} exists")
    test_set = read_test_samples(args.test_samples)
    rows_of_method = {name: join_test_labels(read_predictions(path), test_set)
                      for name, path in zip(names, args.predictions)}
    results = tabulate_results(score_methods(test_set, rows_of_method), compare_methods(args.reference, rows_of_method))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.out, sep="\t", index=False, float_format="%.4f")
    print(results.to_string(index=False))


if __name__ == "__main__":
    main()
