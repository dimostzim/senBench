"""Test of write_final_sets_all980_v10.py (task M.6; after test_write_final_sets_all980_v9.py).
Run: PYTHONPATH=<repository>/scripts:<this folder> python3 -m pytest -q <this file>"""

import pandas as pd
import pytest

from write_final_sets import build_final_sets

# GSE1 is `development`, E-MTAB-1 is `not in set 1`, GSE2 and GSE4 are `test`, GSE3 is left out of both parts. The
# round's labels of GSE2 hold the M.4 row of GSM4 already; M.7 excludes GSM8, the only senescent sample of GSE4.
SAMPLES = pd.DataFrame({
    "study": ["GSE1", "GSE1", "E-MTAB-1", "E-MTAB-1", "GSE2", "GSE2", "GSE3", "GSE3", "GSE4", "GSE4"],
    "dataset": ["GSE1", "GSE1", "E-MTAB-1", "E-MTAB-1", "GSE2", "GSE2", "GSE3", "GSE3", "GSE4", "GSE4"],
    "sample": ["GSM1", "GSM2", "Sample 1", "Sample 2", "GSM3", "GSM4", "GSM5", "GSM6", "GSM8", "GSM10"],
    "state": ["senescent", "not senescent", "senescent", "not senescent", "senescent", "not senescent", "senescent",
              "not senescent", "senescent", "not senescent"],
})
STUDIES = pd.DataFrame({"study": ["GSE1", "E-MTAB-1", "GSE2", "GSE3", "GSE4"], "laboratory_group": ["1", "2", "3", "3", "4"],
                        "part": ["development", "not in set 1", "test", "left out of both parts", "test"]})
CORRECTIONS = pd.DataFrame({
    "dataset": ["GSE1", "GSE1", "E-MTAB-1", "GSE2", "GSE9", "GSE2", "GSE4", "GSE2"],
    "sample": ["GSM1", "GSM2", "Sample 2", "GSM4", "GSM9", "GSM4", "GSM8", "GSM3"],
    "state": ["excluded", "unclear", "unclear", "excluded", "unclear", "not senescent", "excluded", "senescent"],
    "reason": ["r"] * 8,
    "selection": ["label check", "test reading", "label check", "label check", "trial", "test reading", "test reading",
                  "label check"],
    "source": ["L.4", "M.7", "L.4", "L.4", "M.2", "M.4", "M.7", "L.4"],
})


def test_training_data_take_every_correction_and_test_only_the_readings():
    training, test, left = build_final_sets(SAMPLES, STUDIES, CORRECTIONS)
    assert training[["study", "sample", "part", "state", "state_for_training"]].values.tolist() == [
        ["E-MTAB-1", "Sample 1", "not in set 1", "senescent", "senescent"],
        ["E-MTAB-1", "Sample 2", "not in set 1", "not senescent", "excluded"],
        ["GSE1", "GSM1", "development", "senescent", "excluded"],
        ["GSE1", "GSM2", "development", "not senescent", "excluded"]]
    assert test.columns.tolist() == ["study", "dataset", "sample", "laboratory_group", "state"]
    assert test[["sample", "state"]].values.tolist() == [["GSM3", "senescent"], ["GSM4", "not senescent"]]
    assert left == ["GSE4"]


def test_two_corrections_with_another_effect_are_refused_and_two_that_unlabel_agree():
    agreeing = pd.concat([CORRECTIONS, CORRECTIONS.iloc[[0]].assign(state="unclear", selection="test reading",
                                                                    source="M.7")], ignore_index=True)
    build_final_sets(SAMPLES, STUDIES, agreeing)
    other = pd.concat([CORRECTIONS, CORRECTIONS.iloc[[0]].assign(state="senescent", selection="test reading",
                                                                 source="M.7")], ignore_index=True)
    with pytest.raises(ValueError, match="different states"):
        build_final_sets(SAMPLES, STUDIES, other)


def test_an_m4_row_with_another_state_than_the_round_is_refused():
    other = CORRECTIONS.assign(state=CORRECTIONS.state.where(CORRECTIONS.source != "M.4", "senescent"))
    with pytest.raises(ValueError, match="M.4"):
        build_final_sets(SAMPLES, STUDIES, other)


def test_a_trial_correction_of_a_sample_of_the_round_is_refused():
    trial = pd.concat([CORRECTIONS, CORRECTIONS.iloc[[0]].assign(sample="Sample 1", dataset="E-MTAB-1",
                                                                 selection="trial")], ignore_index=True)
    with pytest.raises(ValueError, match="trial"):
        build_final_sets(SAMPLES, STUDIES, trial)
