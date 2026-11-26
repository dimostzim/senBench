"""Check the measures across studies on made-up scores whose right answer is known."""

import itertools

import numpy as np
import pandas as pd
import pytest

import across_studies_measures as measures

STATES = ["senescent"] * 4 + ["not senescent"] * 4


def make_samples(studies, score, groups=None, generator=None):
    """Return a sample table of `studies` studies in their own laboratory groups (or `groups`), 4 positives and 4
    negatives each; `score(study, positive, generator)` gives the score of a sample."""
    generator = generator or np.random.default_rng(1)
    rows = [{"study": f"S{study:03d}", "dataset": f"D{study:03d}", "sample": f"S{study:03d}_{position}",
             "laboratory_group": groups[study] if groups else study, "state": state,
             "score": score(study, state == "senescent", generator), "call": ""}
            for study in range(studies) for position, state in enumerate(STATES)]
    return pd.DataFrame(rows)


def auroc_across(samples):
    table = measures.tabulate_pairs(*measures.split_states(samples), measures.list_studies(samples))
    return measures.estimate_auroc_across_studies(table)[0].value, table


def offset(study):
    return np.sin(study * 7.3) * 50


def perfect_within_unrelated_across(study, positive, generator):
    return offset(study) + (1.0 if positive else 0.0) + generator.normal(0, 0.1)


def one_scale(study, positive, generator):
    return (2.0 if positive else 0.0) + generator.normal(0, 0.3)


def test_perfect_within_studies_and_unrelated_across_gives_one_half():
    value, table = auroc_across(make_samples(60, perfect_within_unrelated_across))
    within_study, _ = measures.compute_reference_values(table)
    assert within_study == 1.0
    assert value == pytest.approx(0.5, abs=0.05)


def test_one_scale_gives_near_one():
    value, _ = auroc_across(make_samples(60, one_scale))
    assert value > 0.98


def test_shuffled_labels_give_one_half():
    assert measures.compute_shuffled_auroc(make_samples(60, one_scale)) == pytest.approx(0.5, abs=0.03)


def rescale(values):
    return (values - values.min()) / (values.max() - values.min())


def test_score_rescaled_within_the_input_is_highest_on_labelled_samples_and_one_half_on_one_state():
    """A score that is perfect within a study and rescaled to 0 to 1 over the samples of its input: on the labelled
    samples it ranks well across studies only because of the composition of the input, on one state alone it does not,
    and with extra samples of the dataset it lies in between. A score that does not depend on the input stays the
    same on all three."""
    generator = np.random.default_rng(3)
    labelled = make_samples(60, perfect_within_unrelated_across, generator=generator)
    extra = pd.DataFrame({"study": labelled.study.unique().repeat(2), "laboratory_group": np.arange(60).repeat(2)})
    extra["score"] = extra.study.str[1:].astype(int).map(offset) + generator.uniform(-2, 3, len(extra))
    rescaled = {
        "labelled samples": labelled.assign(score=labelled.groupby("study").score.transform(rescale)),
        "one state": labelled.assign(score=labelled.groupby(["study", "state"]).score.transform(rescale)),
        "all samples": pd.concat([labelled, extra.assign(state="")]).assign(
            score=lambda table: table.groupby("study").score.transform(rescale)).query("state != ''"),
    }
    values = {name: auroc_across(samples)[0] for name, samples in rescaled.items()}
    assert values["labelled samples"] > 0.95
    assert values["one state"] == pytest.approx(0.5, abs=0.05)
    assert values["labelled samples"] > values["all samples"] > values["one state"]
    assert auroc_across(labelled)[0] == pytest.approx(0.5, abs=0.05)


def test_pairs_within_a_laboratory_group_do_not_count():
    # S0 and S1 share group 0 and rank each other perfectly; across groups every pair is 0.
    scores = {(0, True): 10, (0, False): 9, (1, True): 10, (1, False): 9, (2, True): 0, (2, False): 100}
    samples = make_samples(3, lambda study, positive, generator: scores[(study, positive)], groups=[0, 0, 1])
    value, table = auroc_across(samples)
    _, same_group = measures.compute_reference_values(table)
    assert value == 0.0
    assert same_group == 1.0


def test_resampled_value_counts_pairs_with_the_product_of_multiplicities():
    generator = np.random.default_rng(5)
    groups = [0, 0, 1, 2, 2, 3]
    samples = make_samples(6, lambda study, positive, generator: generator.normal(1 if positive else 0, 1),
                           groups=groups, generator=generator)
    table = measures.tabulate_pairs(*measures.split_states(samples), measures.list_studies(samples))
    _, resampled, _ = measures.resample_auroc_across_studies(table)
    multiplicities = measures.draw_group_multiplicities(tuple(sorted(set(groups))))
    probability = table.wins / table.comparisons
    for draw in range(20):
        # Every copy of a drawn group brings its studies; copies of the same group are never paired.
        copies = [(study, group) for study, group in enumerate(groups) for _ in range(int(multiplicities[draw, group]))]
        pairs = [probability[first, second] for (first, first_group), (second, second_group)
                 in itertools.product(copies, copies) if first_group != second_group]
        expected = np.mean(pairs) if pairs else np.nan
        assert resampled[draw] == pytest.approx(expected, nan_ok=True)


def test_mean_over_studies_is_weighted_by_the_multiplicity_of_the_group():
    values = np.array([1.0, 3.0, np.nan, 5.0])
    groups = np.array([7, 7, 8, 9])
    estimate = measures.estimate_mean_over_studies(values, groups)
    multiplicities = measures.draw_group_multiplicities((7, 8, 9))
    weights = multiplicities[:, [0, 0, 1, 2]]
    defined = ~np.isnan(values)
    with np.errstate(invalid="ignore"):
        expected = (weights[:, defined] * values[defined]).sum(axis=1) / weights[:, defined].sum(axis=1)
    assert estimate.value == 3.0
    low, high = np.nanpercentile(expected, measures.INTERVAL_PERCENTILES)
    assert (estimate.interval_low, estimate.interval_high) == pytest.approx((low, high))
