#!/usr/bin/env python3
"""Write the training data and `test` set of round all980_v10 (task M.6, phase B; after write_final_sets_all980_v9.py),
from the round's split and sample table and kb/data/label_corrections.tsv (fourth version); reads no score and writes
none.

- training_data_samples.tsv: every sample of `development` and of the studies `not in set 1`, with `state` (the
  round's label) and `state_for_training`: the round's label with every correction that names the sample applied
  (`label check`, `trial` and `test reading` rows). A correction to a state other than `senescent` or `not senescent`
  gives `excluded`; two corrections of one sample must have the same effect.
- test_samples.tsv: the labelled samples of `test`, with `state`: the round's label with only the `test reading` rows
  applied (one label column). A `test` dataset left with labelled samples of one state or none leaves `test`, by the
  rule of the earlier rounds.
- final_sets_sha256.txt: the sha256 of both tables.
The `trial` rows are in the round's labels already, so none may name a sample of the round; the `test reading` rows of
M.4 are in the round's labels too, so each of them that names a sample of the round must give its state. Studies left
out of both parts are in neither table. Refuses existing outputs."""

import argparse
import hashlib
from pathlib import Path

import pandas as pd

from accessions import sample_keys
from build_benchmark_bulk import STATES
from build_study_split import DEVELOPMENT, NOT_IN_SET_1, TEST

EXCLUDED = "excluded"
TRIAL = "trial"
TEST_READING = "test reading"
READING_IN_THE_ROUND = "M.4"
TRAINING_NAME = "training_data_samples.tsv"
TEST_NAME = "test_samples.tsv"
CHECKSUMS_NAME = "final_sets_sha256.txt"


def _read(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def correct_states(samples, corrections):
    """Return the state of every sample with the corrections that name it applied; a correction to a state other than
    STATES gives `excluded`, and two corrections of one sample must give the same state then."""
    corrections = corrections[["dataset", "sample", "state"]].merge(samples[["dataset", "sample"]]).assign(
        state=lambda table: table.state.where(table.state.isin(STATES), EXCLUDED)).drop_duplicates()
    if corrections.duplicated(["dataset", "sample"]).any():
        raise ValueError("two corrections of one sample give different states")
    corrected = samples.merge(corrections.rename(columns={"state": "corrected"}), on=["dataset", "sample"],
                              how="left", validate="1:1").corrected
    return corrected.fillna(samples.state.reset_index(drop=True)).values


def validate_corrections(samples, corrections):
    """Refuse a `trial` row that names a sample of the round, an M.4 row with another state than the round's, and a
    correction that names a sample of the round only under another dataset."""
    keys = set(sample_keys(samples.dataset, samples["sample"]))
    trial = corrections[corrections.selection == TRIAL]
    if set(sample_keys(trial.dataset, trial["sample"])) & keys:
        raise ValueError("a `trial` correction names a sample of the round: it should be applied to the labels already")
    in_the_round = corrections[(corrections.selection == TEST_READING) & (corrections.source == READING_IN_THE_ROUND)]
    named = samples.merge(in_the_round[["dataset", "sample", "state"]], on=["dataset", "sample"], suffixes=("", "_read"))
    if (named.state != named.state_read).any():
        raise ValueError(f"a `{TEST_READING}` row of {READING_IN_THE_ROUND} names a sample of the round that has "
                         "another state")
    matched = samples.merge(corrections[["dataset", "sample"]].drop_duplicates(), on=["dataset", "sample"])
    if keys & set(sample_keys(corrections.dataset, corrections["sample"])) != set(
            sample_keys(matched.dataset, matched["sample"])):
        raise ValueError("a correction names a sample of the round only under another dataset")


def keep_datasets_with_both_states(test):
    """Return the labelled samples of the datasets that keep both states, and the datasets that leave."""
    labelled = test[test.state.isin(STATES)]
    states = labelled.groupby("dataset").state.nunique()
    both = set(states.index[states == len(STATES)])
    return labelled[labelled.dataset.isin(both)].reset_index(drop=True), sorted(set(test.dataset) - both)


def build_final_sets(samples, studies, corrections):
    """Return the table of the training data, the table of `test`, and the `test` datasets that leave `test`."""
    validate_corrections(samples, corrections)
    samples = samples.merge(studies[["study", "laboratory_group", "part"]], on="study", how="left", validate="m:1")
    order = ["study", "dataset", "sample"]
    training = samples[samples.part.isin([DEVELOPMENT, NOT_IN_SET_1])].sort_values(order, ignore_index=True)
    training = training.assign(state_for_training=correct_states(training, corrections))
    test = samples[samples.part == TEST].sort_values(order, ignore_index=True)
    test = test.assign(state=correct_states(test, corrections[corrections.selection == TEST_READING]))
    test, left = keep_datasets_with_both_states(test)
    return (training[["study", "dataset", "sample", "part", "laboratory_group", "state", "state_for_training"]],
            test[["study", "dataset", "sample", "laboratory_group", "state"]], left)


def summarise(table, state_column, by):
    """Return per value of `by` (and in all) the studies, laboratory groups and samples per state of a column."""
    rows = []
    for name, rows_of in [*table.groupby(by), ("together", table)]:
        labelled = rows_of[rows_of[state_column].isin(STATES)]
        rows.append({by: name, "studies": labelled.study.nunique(), "laboratory_groups": labelled.laboratory_group.nunique(),
                     "datasets": labelled.dataset.nunique(), "samples": len(labelled),
                     **{state: int((labelled[state_column] == state).sum()) for state in STATES},
                     EXCLUDED: int((rows_of[state_column] == EXCLUDED).sum())})
    return pd.DataFrame(rows)


def compute_sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--split-dir", type=Path, required=True, help="split folder of the round (study_split.tsv)")
    parser.add_argument("--samples", type=Path, required=True, help="benchmark_bulk_samples.tsv of the round")
    parser.add_argument("--label-corrections", type=Path, required=True, help="kb/data/label_corrections.tsv")
    arguments = parser.parse_args()
    out = arguments.split_dir
    for name in (TRAINING_NAME, TEST_NAME, CHECKSUMS_NAME):
        if (out / name).exists():
            raise SystemExit(f"{out / name} exists")
    training, test, left = build_final_sets(_read(arguments.samples), _read(out / "study_split.tsv"),
                                            _read(arguments.label_corrections))
    training.to_csv(out / TRAINING_NAME, sep="\t", index=False)
    test.to_csv(out / TEST_NAME, sep="\t", index=False)
    (out / CHECKSUMS_NAME).write_text("".join(f"{compute_sha256(out / name)}  {name}\n"
                                              for name in (TRAINING_NAME, TEST_NAME)))
    pd.set_option("display.width", 200)
    print(f"=== {TRAINING_NAME}: {len(training)} rows; as labelled (`state`)")
    print(summarise(training, "state", "part").to_string(index=False))
    print("=== with every correction (`state_for_training`)")
    print(summarise(training, "state_for_training", "part").to_string(index=False))
    print(f"=== {TEST_NAME}: {len(test)} rows (`state`); `test` datasets that leave `test` (one state or none): {left}")
    print(summarise(test.assign(part=TEST), "state", "part").to_string(index=False))
    print((out / CHECKSUMS_NAME).read_text(), end="")


if __name__ == "__main__":
    main()
