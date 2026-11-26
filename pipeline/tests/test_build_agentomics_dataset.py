import numpy as np
import pandas as pd
import pytest

from build_agentomics_dataset import (TEST, TRAIN, VALIDATION, assign_ids, choose_genes, draw_validation, list_samples,
                                      name_genes, validate_counts, validate_splits, write_split)

TRAINING = pd.DataFrame({
    "study": ["S1", "S1", "S2", "S2", "S3", "S3", "S4"],
    "dataset": ["GSE1", "GSE1", "GSE2", "GSE2", "E-MTAB-3", "E-MTAB-3", "GSE4"],
    "sample": ["GSM1", "GSM2", "GSM3", "GSM4", "Sample 1", "Sample 2", "GSM5"],
    "part": ["development", "development", "development", "development", "not in set 1", "not in set 1",
             "development"],
    "laboratory_group": ["1", "1", "2", "2", "3", "3", "4"],
    "state": ["senescent", "not senescent", "senescent", "not senescent", "senescent", "not senescent", "senescent"],
    "state_for_training": ["senescent", "not senescent", "senescent", "not senescent", "senescent", "not senescent",
                           "excluded"],
})
TEST_SAMPLES = pd.DataFrame({"study": ["S5", "S5"], "dataset": ["GSE5", "GSE5"], "sample": ["GSM6", "GSM7"],
                             "laboratory_group": ["5", "5"], "state": ["senescent", "not senescent"]})


def test_samples_keep_the_labelled_ones_with_their_keys_and_labels():
    samples = list_samples(TRAINING, TEST_SAMPLES)
    assert samples[["key", "split", "label"]].values.tolist() == [
        ["GSM1", TRAIN, "senescent"], ["GSM2", TRAIN, "not_senescent"], ["GSM3", TRAIN, "senescent"],
        ["GSM4", TRAIN, "not_senescent"], ["E-MTAB-3/Sample 1", TRAIN, "senescent"],
        ["E-MTAB-3/Sample 2", TRAIN, "not_senescent"], ["GSM6", TEST, "senescent"], ["GSM7", TEST, "not_senescent"]]


def test_validation_takes_whole_groups_of_development_only():
    samples = draw_validation(list_samples(TRAINING, TEST_SAMPLES), share=0.5, seed=1)
    in_validation = samples[samples.split == VALIDATION]
    assert set(in_validation.part) == {"development"}
    assert in_validation.laboratory_group.nunique() == 1
    assert set(samples.loc[samples.part == "not in set 1", "split"]) == {TRAIN}
    validate_splits(samples)


def test_a_group_in_two_splits_is_refused():
    samples = list_samples(TRAINING, TEST_SAMPLES)
    samples.loc[samples.key == "GSM2", "split"] = VALIDATION
    with pytest.raises(ValueError, match="laboratory_group is in two splits"):
        validate_splits(samples)


def test_genes_are_chosen_on_train_only():
    counts = pd.DataFrame({"a": [1.0, np.nan, np.nan], "b": [np.nan, 2.0, np.nan], "t": [np.nan, np.nan, 5.0]},
                          index=[10, 20, 30])
    assert choose_genes(counts, ["a", "b"]).tolist() == [10, 20]


@pytest.mark.parametrize("count", [1.5, -1.0])
def test_a_count_that_is_not_a_whole_number_is_refused(count):
    counts = pd.DataFrame({"a": [0.0, 3.0, np.nan]})
    validate_counts(counts)
    with pytest.raises(ValueError, match="negative or not a whole number"):
        validate_counts(counts.assign(a=[count, 3.0, np.nan]))


def test_a_split_holds_the_counts_as_whole_numbers_with_na(tmp_path):
    samples = pd.DataFrame({"key": ["GSM1", "GSM2"], "id": ["sample_0001", "sample_0002"],
                            "label": ["senescent", "not_senescent"]})
    values = pd.DataFrame({"GSM1": [3.0, np.nan], "GSM2": [0.0, 120000.0]}, index=["CDKN1A", "MKI67"])
    write_split(tmp_path / "train", samples, values)
    data = (tmp_path / "train/input/data.csv").read_text()
    assert data == "id,CDKN1A,MKI67\nsample_0001,3,NA\nsample_0002,0,120000\n"
    assert (tmp_path / "train/labels.csv").read_text() == "id,label\nsample_0001,senescent\nsample_0002,not_senescent\n"


def test_genes_are_named_by_symbol_unless_it_is_empty_or_repeated():
    annotation = pd.Series({1: "CDKN1A", 2: "DUP", 3: "DUP"})
    assert name_genes([1, 2, 3, 4], annotation) == ["CDKN1A", "DUP_2", "DUP_3", "GeneID_4"]


def test_ids_are_anonymous_and_unique():
    samples = assign_ids(list_samples(TRAINING, TEST_SAMPLES), seed=1)
    assert samples.id.is_unique and samples.id.str.fullmatch(r"sample_\d{4}").all()
