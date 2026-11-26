"""Check the pooled measures of score.py on made-up scores whose right answer is known, its reading of prediction
files, and that the shipped predictions give the shipped results table."""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import across_studies_measures as measures
from score import (POOLED_AUPRC, POOLED_AUROC, compute_average_precision, join_test_labels, read_predictions,
                   resample_measures)

REPOSITORY = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize("scores, positive, weights, expected", [
    # positives at ranks 1 and 3 of 4: precision 1 at recall 1/2, then 2/3 at recall 1
    ([0.9, 0.8, 0.7, 0.6], [True, False, True, False], [1, 1, 1, 1], (1 + 2 / 3) / 2),
    # one tie of a positive and a negative is one threshold: precision 1/2 at recall 1
    ([0.5, 0.5], [True, False], [1, 1], 0.5),
    # a weight of 2 counts like a copy of the sample
    ([0.9, 0.8, 0.7], [True, False, True], [1, 2, 1], (1 + 2 / 4) / 2),
    # a sample of weight 0 is left out
    ([0.9, 0.8, 0.7], [True, False, True], [1, 0, 1], 1.0),
])
def test_average_precision_matches_hand_values(scores, positive, weights, expected):
    value = compute_average_precision(np.array(scores), np.array(positive), np.array([weights], dtype=float))
    assert value[0] == pytest.approx(expected)


def test_average_precision_without_positive_weight_is_undefined():
    value = compute_average_precision(np.array([0.9, 0.8]), np.array([True, False]), np.array([[0.0, 1.0]]))
    assert np.isnan(value[0])


def test_pooled_auroc_counts_every_pair_and_resamples_like_copies():
    rows = pd.DataFrame({"study": ["A", "A", "B", "B"], "laboratory_group": [1, 1, 2, 2],
                         "state": ["senescent", "not senescent"] * 2, "score": [0.9, 0.2, 0.3, 0.6]})
    studies = measures.list_studies(rows)
    resampled = resample_measures(rows, studies)
    # pairs (positive, negative): (0.9, 0.2) (0.9, 0.6) (0.3, 0.2) (0.3, 0.6): 3 of 4 won
    assert resampled[POOLED_AUROC][0] == pytest.approx(0.75)
    weights = measures.weigh_studies(studies.laboratory_group.to_numpy())
    for draw in range(5):
        copies = rows.loc[rows.index.repeat(np.repeat(weights[draw].astype(int), 2))]
        positives = copies.score[copies.state == "senescent"].to_numpy()
        negatives = copies.score[copies.state == "not senescent"].to_numpy()
        expected = ((positives[:, None] > negatives[None, :]).sum()
                    + 0.5 * (positives[:, None] == negatives[None, :]).sum()) / (len(positives) * len(negatives))
        assert resampled[POOLED_AUROC][1][draw] == pytest.approx(expected)
        assert resampled[POOLED_AUPRC][1][draw] == pytest.approx(compute_average_precision(
            copies.score.to_numpy(), (copies.state == "senescent").to_numpy(), np.ones((1, len(copies))))[0])


@pytest.mark.parametrize("text, message", [
    ("dataset\tsample\tprobability\nD1\tS1\t0.5\n", "the columns are not"),
    ("dataset\tsample\tscore\nD1\tS1\t0.5\nD1\tS1\t0.6\n", "a sample occurs twice"),
    ("dataset\tsample\tscore\nD1\tS1\tinf\n", "not a finite number"),
])
def test_a_prediction_file_names_each_sample_once_with_a_finite_score(tmp_path, text, message):
    path = tmp_path / "method.tsv"
    path.write_text(text)
    with pytest.raises(ValueError, match=message):
        read_predictions(path)


def test_a_predicted_sample_outside_the_test_set_is_refused():
    test_set = pd.DataFrame({"dataset": ["D1"], "sample": ["S1"], "study": ["D1"], "laboratory_group": [1],
                             "state": ["senescent"]})
    predictions = pd.DataFrame({"dataset": ["D1", "D2"], "sample": ["S1", "S1"], "score": [0.5, 0.6]})
    with pytest.raises(ValueError, match="not a labelled sample"):
        join_test_labels(predictions, test_set)


def test_the_shipped_predictions_give_the_shipped_results_table(tmp_path):
    out = tmp_path / "results.tsv"
    predictions = sorted((REPOSITORY / "predictions").glob("*.tsv"))
    subprocess.run([sys.executable, str(REPOSITORY / "scoring" / "score.py"), "--predictions", *map(str, predictions),
                    "--reference", "Agentomics", "--out", str(out)], check=True, capture_output=True)
    assert out.read_bytes() == (REPOSITORY / "results.tsv").read_bytes()
