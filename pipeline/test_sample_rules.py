#!/usr/bin/env python3
"""Check the senescence wording rules and the field mapping of build_sample_table.py on hand-made datasets."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import build_sample_table as b  # noqa: E402


def label(*samples):
    rows = [{"sample": str(i), "title": title, "fields": list(fields), "label_fields": list(fields)}
            for i, (title, fields) in enumerate(samples)]
    b._label_dataset(rows)
    return [(row.get("senescence_state"), row.get("control_kind", ""), row.get("senescence_inducer", "")) for row in rows]


CASES = {
    "stated in titles": (label(("IMR90 senescent rep1", []), ("IMR90 proliferating rep1", [])),
                         [("senescent", "", ""), ("not senescent", "proliferating", "")]),
    "experiment name in every title": (label(("Senescence_timecourse_day0_rep1", []), ("Senescence_timecourse_day10_rep1", [])),
                                       [("unclear", "", ""), ("unclear", "", "")]),
    "negated form": (label(("HDF non-senescent 1", []), ("HDF senescent 1", []), ("HDF pre-senescent 1", [])),
                     [("not senescent", "non-senescent as stated", ""), ("senescent", "", ""),
                      ("not senescent", "non-senescent as stated", "")]),
    "inducer without a stated state": (label(("WI38 IR 10Gy day 10", []), ("WI38 non-irradiated", [])),
                                       [("inducer applied", "", "irradiation"), ("unclear", "", "")]),
    "OIS and vector control": (label(("IMR90 OIS d6 r1", []), ("IMR90 vector control d6 r1", [])),
                               [("senescent", "", "oncogene"), ("not senescent", "empty vector", "")]),
    "field named after senescence": (label(("s1", [("senescence", "yes")]), ("s2", [("senescence", "no")])),
                                     [("senescent", "", ""), ("not senescent", "non-senescent as stated", "")]),
    "both states in one title": (label(("senescent vs proliferating rep1", []), ("other sample", [])),
                                 [("unclear", "", ""), ("unclear", "", "")]),
    "state shared by all samples": (label(("drug A", [("cell state", "senescent")]), ("drug B", [("cell state", "senescent")])),
                                    [("senescent", "", ""), ("senescent", "", "")]),
    "single sample": (label(("replicative senescent HDF", []),), [("senescent", "", "replicative")]),
    "senolytic": (label(("senescent + senolytic ABT263", []), ("senescent + DMSO", [])),
                  [("unclear", "", ""), ("senescent", "", "")]),
    "control word after a stated sample": (label(("etoposide senescent", []), ("DMSO", [])),
                                           [("senescent", "", "DNA-damaging drug"), ("not senescent", "untreated or vehicle", "")]),
    "quiescent": (label(("quiescent 1", []), ("senescent 1", [])), [("not senescent", "quiescent", ""), ("senescent", "", "")]),
    "reagent control is not a control of the state": (label(("etoposide senescent siRELA", []), ("Sen siControl", []), ("DMSO", [])),
                                                      [("senescent", "", "DNA-damaging drug"), ("unclear", "", ""),
                                                       ("not senescent", "untreated or vehicle", "")]),
    "both states are not relabelled as control": (label(("senescent rep1", []), ("early passage proliferating vs senescence", [])),
                                                  [("senescent", "", ""), ("unclear", "", "")]),
    "long description is not a sample field": (label(("a", [("note", "x" * 60 + " cells were followed until replicative senescence")]), ("b", [("note", "other")])),
                                               [("unclear", "", ""), ("unclear", "", "")]),
}
failed = 0
for name, (got, expected) in CASES.items():
    ok = got == expected
    failed += not ok
    print(("ok   " if ok else "FAIL ") + name + ("" if ok else f"\n      got      {got}\n      expected {expected}"))

derived = b._derived([("Characteristics [OrganismPart]", "liver"), ("cell type", "hepatocyte"), ("gender", "Female"),
                      ("age (weeks)", "12"), ("Factor Value[compound]", "doxorubicin"), ("genotype/variation", "WT"),
                      ("characteristics[organism]", "Mus musculus"), ("population doubling", "PD32")])
expected = {"tissue": "liver", "cell_type": "hepatocyte", "sex": "female", "age": "age (weeks): 12",
            "treatment": "doxorubicin", "genotype": "WT", "organism": "Mus musculus", "passage": "PD32"}
ok = derived == expected
failed += not ok
print(("ok   " if ok else "FAIL ") + "field mapping" + ("" if ok else f"\n      got {derived}"))
sys.exit(1 if failed else 0)
