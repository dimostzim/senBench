"""Measures of the bulk benchmark across studies, with intervals from resamples of laboratory groups.

A sample table here has one row per scored sample with the columns study, laboratory_group, dataset, sample, state,
score and call. Every study weighs the same; a pair of studies counts only when the studies are in different
laboratory groups.
"""

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import pandas as pd

from benchmark_bulk_classes import NOT_SENESCENT, SENESCENT
from score_benchmark_bulk import compute_auroc

SEED = 20261116
RESAMPLES = 2000
SHUFFLES = 20
COMPANION_DRAWS = 200
INTERVAL_PERCENTILES = (2.5, 97.5)


@dataclass
class Estimate:
    """A value with its percentile interval over the resamples; `resamples_without_value` had no study or pair."""
    value: float
    interval_low: float
    interval_high: float
    resamples_without_value: int


@dataclass
class PairTable:
    """For every ordered pair of studies (row S gives positives, column T gives negatives): the sum of the
    comparisons of a positive of S with a negative of T (a win is 1, a tie 0.5) and their number."""
    studies: np.ndarray
    groups: np.ndarray
    wins: np.ndarray
    comparisons: np.ndarray

    def mean_over(self, pairs):
        """Return the mean over the given pairs of the probability that the positive scores higher."""
        pairs = pairs & (self.comparisons > 0)
        return float((self.wins[pairs] / self.comparisons[pairs]).mean()) if pairs.any() else np.nan

    def across_groups(self):
        return self.groups[:, None] != self.groups[None, :]


@lru_cache(maxsize=256)
def draw_group_multiplicities(groups):
    """Return how often each laboratory group is drawn in each resample (resamples x groups, groups sorted).

    The draws depend only on the tuple of groups, so two measures on the same groups use the same resamples.
    """
    generator = np.random.default_rng(SEED)
    draws = generator.integers(0, len(groups), size=(RESAMPLES, len(groups)))
    multiplicities = np.zeros((RESAMPLES, len(groups)))
    np.add.at(multiplicities, (np.arange(RESAMPLES)[:, None], draws), 1)
    return multiplicities


def weigh_studies(study_groups):
    """Return the weight of every study in every resample: the number of times its laboratory group was drawn."""
    groups, codes = np.unique(study_groups, return_inverse=True)
    return draw_group_multiplicities(tuple(groups.tolist()))[:, codes]


def summarise_resamples(value, resampled):
    defined = ~np.isnan(resampled)
    low, high = np.percentile(resampled[defined], INTERVAL_PERCENTILES) if defined.any() else (np.nan, np.nan)
    return Estimate(value, float(low), float(high), int((~defined).sum()))


def estimate_mean_over_studies(values, study_groups):
    """Return the mean over studies of one value per study (NaN: the study has none), with its interval."""
    values = np.asarray(values, dtype=float)
    defined = ~np.isnan(values)
    weights = weigh_studies(np.asarray(study_groups))
    with np.errstate(invalid="ignore", divide="ignore"):
        resampled = (weights @ np.where(defined, values, 0)) / (weights @ defined)
    return summarise_resamples(float(values[defined].mean()) if defined.any() else np.nan, resampled)


def tabulate_pairs(positives, negatives, studies):
    """Return the PairTable of the positives and negatives over the given studies (a table with study and group)."""
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
    """Return the mean over ordered pairs of studies of different laboratory groups, its value in every resample,
    and the number of pairs.

    In a resample a pair counts with the product of the numbers of times its two groups were drawn; copies of the
    same group are never paired, because pairs within one group do not count.
    """
    valid = table.across_groups() & (table.comparisons > 0)
    probability = np.where(valid, table.wins / np.where(valid, table.comparisons, 1), 0)
    weights = weigh_studies(table.groups)
    numerator = ((weights @ probability) * weights).sum(axis=1)
    denominator = ((weights @ valid.astype(float)) * weights).sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        resampled = numerator / denominator
    return table.mean_over(valid), resampled, int(valid.sum())


def estimate_auroc_across_studies(table):
    """Return the AUROC across studies with its interval, and the number of pairs."""
    value, resampled, pairs = resample_auroc_across_studies(table)
    return summarise_resamples(value, resampled), pairs


def compute_reference_values(table):
    """Return the mean over pairs within one study, and over pairs of different studies of one laboratory group."""
    same_study = np.eye(len(table.studies), dtype=bool)
    same_group = ~table.across_groups() & ~same_study
    return table.mean_over(same_study), table.mean_over(same_group)


def split_states(samples):
    return samples[samples.state == SENESCENT], samples[samples.state == NOT_SENESCENT]


def list_studies(samples):
    return samples.drop_duplicates("study")[["study", "laboratory_group"]].sort_values("study").reset_index(drop=True)


def compute_shuffled_auroc(samples):
    """Return the mean AUROC across studies over shuffles of the labels within each study (expected 0.5)."""
    generator = np.random.default_rng(SEED)
    samples = samples.sort_values(["study", "dataset", "sample"]).reset_index(drop=True)
    studies = list_studies(samples)
    codes = pd.Index(studies.study).get_indexer(samples.study)
    scores = samples.score.to_numpy(dtype=float)
    positive = (samples.state == SENESCENT).to_numpy()
    values = []
    for _ in range(SHUFFLES):
        # The samples are sorted by study, so sorting by study and a random key permutes the labels within studies.
        shuffled = positive[np.lexsort((generator.random(len(samples)), codes))]
        table = _tabulate(scores[shuffled], codes[shuffled], scores[~shuffled], codes[~shuffled],
                          studies.study.to_numpy(), studies.laboratory_group.to_numpy())
        values.append(table.mean_over(table.across_groups()))
    return float(np.nanmean(values))


def compute_companion_random_groups(table, groups_in_class):
    """Return the mean AUROC across studies over random draws of as many laboratory groups as a class has."""
    generator = np.random.default_rng(SEED)
    groups = np.unique(table.groups)
    across = table.across_groups()
    values = []
    for _ in range(COMPANION_DRAWS):
        chosen = np.isin(table.groups, generator.choice(groups, size=groups_in_class, replace=False))
        values.append(table.mean_over(across & chosen[:, None] & chosen[None, :]))
    return float(np.nanmean(values))


def contribute_studies(table):
    """Return per study the mean over its valid pairs as S (its positives) and as T (its negatives)."""
    valid = table.across_groups() & (table.comparisons > 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        probability = np.where(valid, table.wins / np.where(valid, table.comparisons, 1), 0)
        return pd.DataFrame({
            "study": table.studies,
            "pairs_as_s": valid.sum(axis=1), "mean_as_s": probability.sum(axis=1) / valid.sum(axis=1),
            "pairs_as_t": valid.sum(axis=0), "mean_as_t": probability.sum(axis=0) / valid.sum(axis=0)})


def compute_calls_per_study(samples):
    """Return per study the sensitivity, the specificity and their mean (balanced accuracy)."""
    positives, negatives = split_states(samples)
    studies = list_studies(samples).set_index("study")
    studies["sensitivity"] = (positives.call == SENESCENT).groupby(positives.study).mean()
    studies["specificity"] = (negatives.call == NOT_SENESCENT).groupby(negatives.study).mean()
    studies["balanced accuracy"] = (studies.sensitivity + studies.specificity) / 2
    return studies


def compute_auroc_within_datasets(samples):
    """Return per study the mean over its datasets of the AUROC within the dataset (datasets with both states)."""
    per_dataset = samples.groupby(["study", "laboratory_group", "dataset"]).apply(
        lambda rows: compute_auroc(rows.score.astype(float), rows.state == SENESCENT)
        if rows.state.nunique() == 2 else np.nan, include_groups=False).rename("auroc").reset_index()
    return per_dataset.groupby(["study", "laboratory_group"]).auroc.mean().reset_index()


def compute_change_per_dataset(changed, labelled):
    """Return per dataset the median absolute change of the score divided by the standard deviation of the
    dataset's scores on the labelled samples (NaN where that is 0), and the share of samples whose call changes."""
    both = changed.merge(labelled[["dataset", "sample", "score", "call"]], on=["dataset", "sample"],
                         suffixes=("", "_labelled"), validate="1:1")
    spread = labelled.groupby("dataset").score.std()
    per_dataset = both.groupby(["study", "laboratory_group", "dataset"]).apply(lambda rows: pd.Series({
        "median_change": (rows.score - rows.score_labelled).abs().median(),
        "call_changed": (rows.call != rows.call_labelled).mean()}), include_groups=False).reset_index()
    with np.errstate(invalid="ignore", divide="ignore"):
        per_dataset["change"] = per_dataset.median_change / per_dataset.dataset.map(spread).replace(0, np.nan)
    return per_dataset
