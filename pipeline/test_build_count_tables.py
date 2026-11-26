#!/usr/bin/env python3
"""Check on hand-made frames which values build_count_tables.py takes as counts, which it rounds as expected counts
and which fail, and how it gives the columns of a GEO uniform table to samples by the table of reassignments."""

import sys
import tempfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import build_count_tables as b  # noqa: E402

GENES = ["ENSG00000000003", "ENSG00000000005", "ENSG00000000419", "ENSG00000000457", "ENSG00000000460"]
TRANSCRIPTS = ["ENST00000373020.9", "ENST00000494424.1", "NM_000014", "XR_001737835.2", "ENSG00000000460"]


def check(identifiers, **samples):
    """Return the checked values per sample and the number of values rounded, or the problem found."""
    frame = pd.DataFrame({b.GENE_COLUMN: identifiers, **samples})
    try:
        counts, _, rounded = b._check_counts(frame, "file")
    except b._CheckProblem as problem:
        return str(problem)
    return {sample: list(counts[sample]) for sample in samples}, rounded


CASES = {
    "whole numbers stay": (check(GENES, s1=["0", "1", "7", "1", "12"]), ({"s1": [0, 1, 7, 1, 12]}, 0)),
    "whole numbers written with decimals stay": (check(GENES, s1=["0.00", "1.00", "7.0", "1", "12"]),
                                                 ({"s1": [0, 1, 7, 1, 12]}, 0)),
    "expected counts are rounded half to even": (check(GENES, s1=["0.5", "1", "1.5", "1", "2.5"],
                                                       s2=["3", "3", "10.49", "0", "10.51"]),
                                                 ({"s1": [0, 1, 2, 1, 2], "s2": [3, 3, 10, 0, 11]}, 5)),
    "a value next to a half is rounded by its text": (
        check(GENES, s1=["3653.4999999999995", "1964.5000000000002", "487.49999999999994", "1", "1"]),
        ({"s1": [3653, 1965, 487, 1, 1]}, 3)),
    "a sample without a value above zero does not fail expected counts": (
        check(GENES, s1=["0", "0", "0", "0", "0"], s2=["2", "2", "0.4", "0", "7.7"]),
        ({"s1": [0, 0, 0, 0, 0], "s2": [2, 2, 0, 0, 8]}, 2)),
    "counts divided by a factor fail": (
        check(GENES, s1=["0", "0.9553", "0.9553", "1.9106", "23.88"]),
        "file: values are not whole numbers, and the most common value above zero of s1 is 0.9553, not a whole "
        "number"),
    "one scaled sample fails the table": (
        check(GENES, s1=["1", "1", "2.5", "0", "4"], s2=["0.01", "0.01", "3.2", "0", "1"]),
        "file: values are not whole numbers, and the most common value above zero of s2 is 0.01, not a whole number"),
    "expected counts per transcript fail": (
        check(TRANSCRIPTS, s1=["1", "1", "2.5", "0", "4"]),
        "file: values are not whole numbers, and 4 of 5 rows are named by a transcript accession"),
    "negative values fail": (check(GENES, s1=["1", "1", "-2.5", "0", "4"]),
                             "file: 1 values are negative, the first '-2.5' of s1"),
    "values that are not numbers fail": (check(GENES, s1=["1", "1", "NA", "0", "4.5"]),
                                         "file: 1 values are not numbers, the first 'NA' of s1"),
}



def read_uniform(source_columns):
    """Return the columns per sample, the match per sample and the columns that are not samples, from a made-up GEO
    uniform table with the columns GSM1, GSM2, GSM3 and GSM9 (GSM9 is not a sample of the dataset)."""
    with tempfile.TemporaryDirectory() as folder:
        path = Path(folder) / "uniform.tsv.gz"
        pd.DataFrame({"GeneID": ["1", "2"], "GSM1": ["1", "10"], "GSM2": ["2", "20"], "GSM3": ["3", "30"],
                      "GSM9": ["9", "90"]}).to_csv(path, sep="\t", index=False)
        frame, _, matches, others = b._read_uniform_table(path, dict.fromkeys(["GSM1", "GSM2", "GSM3"], ""),
                                                          source_columns)
    return ({sample: list(frame[sample]) for sample in frame.columns[1:]},
            {match["sample"]: (match["column_name"], match["match"]) for match in matches}, others)


def read_reassignments(rows):
    with tempfile.TemporaryDirectory() as folder:
        path = Path(folder) / "reassignments.tsv"
        pd.DataFrame(rows, columns=["dataset", "sample", "outcome", "count_source_column"]).to_csv(path, sep="\t",
                                                                                                  index=False)
        try:
            return b.read_reassignments(path)
        except ValueError as problem:
            return str(problem).split(": ")[1]


CASES.update({
    "every sample takes its own GEO column": (
        read_uniform({"GSM1": "GSM1", "GSM2": "GSM2", "GSM3": "GSM3"}),
        ({"GSM1": [1, 10], "GSM2": [2, 20], "GSM3": [3, 30]},
         {"GSM1": ("GSM1", "accession"), "GSM2": ("GSM2", "accession"), "GSM3": ("GSM3", "accession")}, 1)),
    "exchanged GEO columns go to their samples, and a sample left out has no column": (
        read_uniform({"GSM1": "GSM2", "GSM2": "GSM1"}),
        ({"GSM1": [2, 20], "GSM2": [1, 10]}, {"GSM1": ("GSM2", "reassigned"), "GSM2": ("GSM1", "reassigned")}, 1)),
    "the table of reassignments gives the columns, and the datasets of submitters' raw counts apart": (
        read_reassignments([("D1", "GSM1", "reassigned", "GSM2"), ("D1", "GSM2", "reassigned", "GSM1"),
                            ("D1", "GSM3", "left out", ""), ("D2", "GSM5", "kept", "GSM5"),
                            ("D3", "GSM7", "left out", ""), ("D4", "GSM8", "submitters' raw counts", "c1")]),
        ({"D4"}, {"D1": {"GSM1": "GSM2", "GSM2": "GSM1"}, "D2": {"GSM5": "GSM5"}, "D3": {}})),
    "a GEO column given to two samples fails": (
        read_reassignments([("D1", "GSM1", "kept", "GSM1"), ("D1", "GSM2", "reassigned", "GSM1")]),
        "a GEO column is given to several samples"),
})


def build_counted():
    """Return the units row and the match per sample of dataset D1 (samples GSM1, GSM2, GSM3), built from a made-up
    table counted from raw reads that holds GSM1 and GSM2."""
    with tempfile.TemporaryDirectory() as folder:
        folder = Path(folder)
        (folder / "counted").mkdir()
        (folder / "out" / b.RECORDS_FOLDER).mkdir(parents=True)
        pd.DataFrame({"GeneID": ["1", "2"], "GSM1": ["1", "10"], "GSM2": ["2", "20"]}).to_csv(
            folder / "counted" / "D1.tsv.gz", sep="\t", index=False)
        unit = b._build_unit({"D1": dict.fromkeys(["GSM1", "GSM2", "GSM3"], "")}, {}, folder, folder, folder, set(),
                             folder / "out", False, set(), {}, folder / "counted", "D1")
        record = (folder / "out" / b.RECORDS_FOLDER / "D1.json").read_text()
    return unit, {column["sample"]: column["match"] for column in b.json.loads(record)["columns"]}


CASES["a table counted from raw reads gives its columns to its samples, and a sample without one has no column"] = (
    build_counted(),
    (("D1", "ok", "counted from raw reads; 2 of 3 samples"),
     {"GSM1": "accession", "GSM2": "accession", "GSM3": "no column"}))

failed = 0
for name, (got, expected) in CASES.items():
    ok = got == expected
    failed += not ok
    print(("ok   " if ok else "FAIL ") + name + ("" if ok else f"\n      got      {got}\n      expected {expected}"))
sys.exit(1 if failed else 0)
