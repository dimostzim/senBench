#!/usr/bin/env python3
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

SENESCENT = "senescent"
NOT_SENESCENT = "not senescent"
AUROC_ACROSS = "AUROC across studies"
POOLED_AUROC = "pooled AUROC"
POOLED_AUPRC = "pooled AUPRC"
MEASURES = (AUROC_ACROSS, POOLED_AUROC, POOLED_AUPRC)
PREDICTION_COLUMNS = ["dataset", "sample", "score"]
ROOT = Path(__file__).resolve().parent.parent
TEST_SAMPLES = ROOT / "pipeline" / "test_samples.tsv"
PREDICTIONS = ROOT / "predictions"
RESULTS = ROOT / "results.tsv"
REFERENCE = "Agentomics"
RESULT_COLUMNS = ["method", "measure", "value", "interval_low", "interval_high", "random_score", "samples", "senescent",
                  "not_senescent", "datasets", "laboratory_groups", "reference_minus_method", "difference_interval_low",
                  "difference_interval_high", "difference_samples"]
SEED = 20261116
RESAMPLES = 2000
INTERVAL_PERCENTILES = (2.5, 97.5)


@dataclass
class PairTable:
    studies: np.ndarray
    groups: np.ndarray
    wins: np.ndarray
    comparisons: np.ndarray

    def mean_over(self, pairs):
        pairs = pairs & (self.comparisons > 0)
        return float((self.wins[pairs] / self.comparisons[pairs]).mean()) if pairs.any() else np.nan

    def across_groups(self):
        return self.groups[:, None] != self.groups[None, :]


@lru_cache(maxsize=256)
def draw_group_multiplicities(groups):
    generator = np.random.default_rng(SEED)
    draws = generator.integers(0, len(groups), size=(RESAMPLES, len(groups)))
    multiplicities = np.zeros((RESAMPLES, len(groups)))
    np.add.at(multiplicities, (np.arange(RESAMPLES)[:, None], draws), 1)
    return multiplicities


def weigh_studies(study_groups):
    groups, codes = np.unique(study_groups, return_inverse=True)
    return draw_group_multiplicities(tuple(groups.tolist()))[:, codes]


def summarise_resamples(value, resampled):
    defined = ~np.isnan(resampled)
    low, high = np.percentile(resampled[defined], INTERVAL_PERCENTILES) if defined.any() else (np.nan, np.nan)
    return value, float(low), float(high)


def tabulate_pairs(positives, negatives, studies):
    index = pd.Index(studies.study)
    positive_codes = index.get_indexer(positives.study)
    negative_codes = index.get_indexer(negatives.study)
    if (positive_codes < 0).any() or (negative_codes < 0).any():
        raise ValueError("a sample belongs to a study outside the analysis")
    return _tabulate(positives.score.to_numpy(dtype=float), positive_codes, negatives.score.to_numpy(dtype=float),
                     negative_codes, studies.study.to_numpy(), studies.laboratory_group.to_numpy())


def _tabulate(positive_scores, positive_codes, negative_scores, negative_codes, names, groups):
    comparison = ((positive_scores[:, None] > negative_scores[None, :])
                  + 0.5 * (positive_scores[:, None] == negative_scores[None, :]))
    positive_onehot = np.zeros((len(names), len(positive_scores)))
    positive_onehot[positive_codes, np.arange(len(positive_scores))] = 1
    negative_onehot = np.zeros((len(names), len(negative_scores)))
    negative_onehot[negative_codes, np.arange(len(negative_scores))] = 1
    return PairTable(names, groups, positive_onehot @ comparison @ negative_onehot.T,
                     np.outer(positive_onehot.sum(axis=1), negative_onehot.sum(axis=1)))


def resample_auroc_across_studies(table):
    valid = table.across_groups() & (table.comparisons > 0)
    probability = np.where(valid, table.wins / np.where(valid, table.comparisons, 1), 0)
    weights = weigh_studies(table.groups)
    numerator = ((weights @ probability) * weights).sum(axis=1)
    denominator = ((weights @ valid.astype(float)) * weights).sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        resampled = numerator / denominator
    return table.mean_over(valid), resampled


def split_states(samples):
    return samples[samples.state == SENESCENT], samples[samples.state == NOT_SENESCENT]


def list_studies(samples):
    return samples.drop_duplicates("study")[["study", "laboratory_group"]].sort_values("study").reset_index(drop=True)


def read_table(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def read_test_samples(path):
    test_set = read_table(path)
    test_set["laboratory_group"] = test_set.laboratory_group.astype(int)
    return test_set[["dataset", "sample", "study", "laboratory_group", "state"]]


def read_predictions(path):
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
    rows = predictions.merge(test_set, on=["dataset", "sample"], how="left", validate="1:1", indicator=True)
    if (rows._merge != "both").any():
        raise ValueError("a predicted sample is not a labelled sample of `test`")
    return rows.drop(columns="_merge")


def count_entered(rows):
    return {"laboratory_groups": rows.laboratory_group.nunique(), "datasets": rows.dataset.nunique(),
            "samples": len(rows), "senescent": int((rows.state == SENESCENT).sum()),
            "not_senescent": int((rows.state == NOT_SENESCENT).sum())}


def compute_average_precision(scores, positive, weights):
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
    table = tabulate_pairs(*split_states(rows), studies)
    across, across_draws = resample_auroc_across_studies(table)
    study_weights = weigh_studies(table.groups)
    with np.errstate(invalid="ignore", divide="ignore"):
        pooled_draws = (((study_weights @ table.wins) * study_weights).sum(axis=1)
                        / ((study_weights @ table.comparisons) * study_weights).sum(axis=1))
    scores = rows.score.to_numpy(dtype=float)
    positive = (rows.state == SENESCENT).to_numpy()
    precision_draws = compute_average_precision(scores, positive,
                                                weigh_studies(rows.laboratory_group.to_numpy()))
    return {AUROC_ACROSS: (across, across_draws),
            POOLED_AUROC: (float(table.wins.sum() / table.comparisons.sum()), pooled_draws),
            POOLED_AUPRC: (float(compute_average_precision(scores, positive, np.ones((1, len(rows))))[0]),
                           precision_draws)}


def describe_random_value(measure, rows):
    return float((rows.state == SENESCENT).mean()) if measure == POOLED_AUPRC else 0.5


def score_methods(rows_of_method):
    results = []
    for method, rows in rows_of_method.items():
        for measure, (value, draws) in resample_measures(rows, list_studies(rows)).items():
            value, low, high = summarise_resamples(value, draws)
            results.append({"method": method, "measure": measure, "value": value, "interval_low": low,
                            "interval_high": high, "random_score": describe_random_value(measure, rows),
                            **count_entered(rows)})
    return pd.DataFrame(results)


def compare_methods(rows_of_method):
    results = []
    for method, rows in rows_of_method.items():
        if method == REFERENCE:
            continue
        common = rows_of_method[REFERENCE][["dataset", "sample"]].merge(rows[["dataset", "sample"]])
        kept = [of_method.merge(common) for of_method in (rows_of_method[REFERENCE], rows)]
        studies = list_studies(kept[0])
        resampled = [resample_measures(of_method, studies) for of_method in kept]
        for measure in resampled[0]:
            (reference_value, reference_draws), (method_value, method_draws) = (
                of_method[measure] for of_method in resampled)
            difference, low, high = summarise_resamples(reference_value - method_value,
                                                        reference_draws - method_draws)
            results.append({"method": method, "measure": measure, "reference_minus_method": difference,
                            "difference_interval_low": low, "difference_interval_high": high,
                            "difference_samples": len(common)})
    return pd.DataFrame(results)


def tabulate_results(scored, differences):
    table = scored.merge(differences, on=["method", "measure"], how="left")
    order = scored[scored.measure == AUROC_ACROSS].sort_values("value", ascending=False, kind="stable").method
    table["method"] = pd.Categorical(table.method, categories=order)
    table["measure"] = pd.Categorical(table.measure, categories=MEASURES)
    table["difference_samples"] = table.difference_samples.astype("Int64")
    return table.sort_values(["method", "measure"])[RESULT_COLUMNS]


def main():
    test_set = read_test_samples(TEST_SAMPLES)
    rows_of_method = {path.stem: join_test_labels(read_predictions(path), test_set)
                      for path in sorted(PREDICTIONS.glob("*.tsv"))}
    results = tabulate_results(score_methods(rows_of_method), compare_methods(rows_of_method))
    results.to_csv(RESULTS, sep="\t", index=False, float_format="%.4f")
    print(results.to_string(index=False))


if __name__ == "__main__":
    main()
