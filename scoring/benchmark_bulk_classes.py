#!/usr/bin/env python3
"""Reduce the inducer, the cell type and the time after induction of the sample labels to the classes of the design.

Run as a script, it writes the class of every distinct pair of `cell_line` and `cell_type` and of every distinct
`time_after_induction` of a sample table, for a person to check.
"""

import argparse
import re
from pathlib import Path

import pandas as pd

SENESCENT = "senescent"
NOT_SENESCENT = "not senescent"
QUIESCENT = "quiescent"

OTHER_INDUCER = "other"
# The order also breaks ties between classes with the same number of positives.
INDUCER_CLASSES = ("replicative", "irradiation", "oncogene", "DNA-damaging drug", "CDK4/6 inhibitor", "oxidative stress",
                   OTHER_INDUCER)
INDUCERS_IN_OTHER = ("therapy", "other", "not stated", "none")

OTHER_CELL_TYPE = "other or not stated"
# Class -> keys; the order of the classes decides when several match, and breaks ties between classes.
CELL_TYPE_KEYS = {
    "fibroblast": ("fibroblast", "fibroblasts", "imr-90", "wi-38", "bj", "hff", "mrc-5", "hdf", "lf1", "tig-1", "tig-3"),
    "endothelial": ("endothelial", "huvec", "haec", "hmvec", "hrmec", "hmec-1"),
    "mesenchymal stromal cell": ("mesenchymal", "msc", "mscs", "stromal", "adsc", "adipose-derived"),
    "epithelial or keratinocyte": ("epithelial", "keratinocyte", "keratinocytes", "rpe", "arpe-19", "hmec", "rptec",
                                   "saec", "hepatocyte", "hepatocytes", "mcf-10a", "rpe1", "cholangiocyte",
                                   "cholangiocytes"),
    "cancer cell line": ("cancer", "carcinoma", "adenocarcinoma", "melanoma", "tumor", "tumour", "leukemia", "lymphoma",
                         "glioma", "glioblastoma", "sarcoma", "osteosarcoma", "malignant", "a549", "mcf-7", "hct-116",
                         "hela", "u2os", "hepg2", "mda-mb-231", "a375", "lovo", "h1299", "sk-mel-28", "huh-7"),
}
# Tumour words that also count at the end of a longer word (liposarcoma, cholangiocarcinoma, xanthoastrocytoma).
TUMOUR_WORD_ENDINGS = ("carcinoma", "sarcoma", "blastoma", "cytoma", "glioma", "melanoma", "lymphoma", "leukemia",
                       "leukaemia", "mesothelioma")
CANCER_CELL_LINE = "cancer cell line"
CELL_TYPE_CLASSES = (*CELL_TYPE_KEYS, OTHER_CELL_TYPE)

NO_TIME = "not stated or not a time"
# Upper bound in days (included) -> class; a longer time is the last class.
TIME_CLASS_BOUNDS = ((4, "4 days or less"), (9, "more than 4 to 9 days"), (20, "more than 9 to 20 days"))
LONGEST_TIME_CLASS = "more than 20 days"
TIME_CLASSES = (*(name for _, name in TIME_CLASS_BOUNDS), LONGEST_TIME_CLASS, NO_TIME)
DAYS_PER_UNIT = {"h": 1 / 24, "d": 1, "w": 7, "m": 30}
NUMBER_WORDS = {"one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6, "seven": 7, "eight": 8, "nine": 9,
                "ten": 10, "eleven": 11, "twelve": 12}
_NUMBER = r"\d+(?:\.\d+)?|" + "|".join(NUMBER_WORDS)
# "7 days", "7d", "72 h", "2 weeks", "2 wk", "1 month", "five days"; or the unit first: "day 7", "d7", "week 1", "month 2".
TIME = re.compile(rf"(?<![\w.])(?P<number>{_NUMBER})\s*(?P<unit>hours?|hrs?|h|days?|d|weeks?|wks?|months?)(?![a-z])"
                  rf"|(?<![a-z])(?P<unit_first>day|d|week|wk|month)\s*(?P<number_after>{_NUMBER})(?![\w.])")


def _whole_word(key):
    """Return a pattern for a key as a whole word; a hyphen or a space inside the key is optional."""
    return re.compile(r"(?<![a-z0-9])" + re.sub(r"[- ]", "[- ]?", key) + r"(?![a-z0-9])")


def names_key(text, key):
    """Return whether a text holds a key as a whole word, in any case."""
    return bool(_whole_word(key).search(text.lower()))


_CELL_TYPE_PATTERNS = {name: [_whole_word(key) for key in keys] for name, keys in CELL_TYPE_KEYS.items()}
_CELL_TYPE_PATTERNS[CANCER_CELL_LINE].append(re.compile(r"[a-z]*(?:" + "|".join(TUMOUR_WORD_ENDINGS) + r")(?![a-z0-9])"))


def classify_inducer(inducer):
    if inducer in INDUCERS_IN_OTHER:
        return OTHER_INDUCER
    if inducer not in INDUCER_CLASSES:
        raise ValueError(f"inducer outside the allowed list: {inducer!r}")
    return inducer


def classify_cell_type(cell_line, cell_type):
    """Return the first class with a key in `cell_line`; if none has one, the first class with a key in `cell_type`."""
    for text in (cell_line.lower(), cell_type.lower()):
        for name, patterns in _CELL_TYPE_PATTERNS.items():
            if any(pattern.search(text) for pattern in patterns):
                return name
    return OTHER_CELL_TYPE


def read_days(time_after_induction):
    """Return the first time named in the text, in days, or None if the text names none."""
    found = TIME.search(time_after_induction.lower())
    if not found:
        return None
    number, unit = (found["number"], found["unit"]) if found["number"] else (found["number_after"], found["unit_first"])
    return float(NUMBER_WORDS.get(number, number)) * DAYS_PER_UNIT[unit[0]]


def classify_time(time_after_induction):
    days = read_days(time_after_induction)
    if days is None:
        return NO_TIME
    return next((name for bound, name in TIME_CLASS_BOUNDS if days <= bound), LONGEST_TIME_CLASS)


def add_classes(samples):
    """Return the sample table with the columns inducer_class, cell_type_class and time_class.

    The inducer class and the time class describe a positive; they are empty for a negative.
    """
    positive = samples.state == SENESCENT
    cell_type_class = [classify_cell_type(line, kind) for line, kind in zip(samples.cell_line, samples.cell_type)]
    return samples.assign(
        inducer_class=samples.inducer.map(classify_inducer).where(positive, ""),
        cell_type_class=cell_type_class,
        time_class=samples.time_after_induction.map(classify_time).where(positive, ""))


def tabulate_cell_type_classes(samples):
    """Return one row per distinct pair of cell_line and cell_type with its class and its number of samples."""
    table = add_classes(samples).groupby(["cell_line", "cell_type", "cell_type_class"]).size().rename("samples")
    return table.reset_index().sort_values(["cell_type_class", "samples"], ascending=[True, False])


def tabulate_time_classes(samples):
    """Return one row per distinct time_after_induction of the positives with its days, its class and its positives."""
    positives = add_classes(samples[samples.state == SENESCENT])
    table = positives.groupby(["time_after_induction", "time_class"]).size().rename("positives").reset_index()
    table.insert(1, "days", table.time_after_induction.map(read_days))
    return table.sort_values(["time_class", "positives"], ascending=[True, False])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=Path, required=True, help="benchmark_bulk_samples.tsv")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    cell_type_path = args.out_dir / "cell_type_classes.tsv"
    time_path = args.out_dir / "time_classes.tsv"
    for path in (cell_type_path, time_path):
        if path.exists():
            raise SystemExit(f"{path} exists")
    samples = pd.read_csv(args.samples, sep="\t", dtype=str, keep_default_na=False)
    cell_types = tabulate_cell_type_classes(samples)
    times = tabulate_time_classes(samples)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    cell_types.to_csv(cell_type_path, sep="\t", index=False)
    times.to_csv(time_path, sep="\t", index=False)
    other = cell_types[cell_types.cell_type_class == OTHER_CELL_TYPE].samples.sum()
    no_time = times[times.time_class == NO_TIME].positives.sum()
    print(f"{other} of {len(samples)} samples in `{OTHER_CELL_TYPE}`; "
          f"{no_time} of {times.positives.sum()} positives in `{NO_TIME}`")


if __name__ == "__main__":
    main()
