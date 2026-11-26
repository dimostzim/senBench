#!/usr/bin/env python3
"""Build the laboratory groups, the sets and the split of the bulk benchmark across studies into `test` and `development`."""

import argparse
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from benchmark_bulk_classes import (CELL_TYPE_CLASSES, INDUCER_CLASSES, NOT_SENESCENT, QUIESCENT, SENESCENT, add_classes,
                                    names_key)
from build_method_source_datasets import APPLIED_ONLY, COUNTED_ROLES, validate_roles
from laboratory_groups import build_laboratory_groups, rank_accessions

SEED = 20261116
TEST = "test"
DEVELOPMENT = "development"
NOT_IN_SET_1 = "not in set 1"
LEFT_OUT_OF_BOTH_PARTS = "left out of both parts"
PARTS = (TEST, DEVELOPMENT)
AMENDED_PARTS = (DEVELOPMENT, LEFT_OUT_OF_BOTH_PARTS)
DATASET_SEPARATOR = ";"
SPLIT_NAME = "study_split.tsv"
RECORD_NAME = "study_split_record.json"
SERIES_NAME = "laboratory_group_series.tsv"
MATCHES_NAME = "list_matches.tsv"
SUMMARY_NAME = "split_summary.tsv"
STRATA_NAME = "split_strata.tsv"
SENCAT_CANDIDATES_NAME = "sencat_unnamed_dataset_candidates.tsv"
SENCAT_CELL_LINE = "imr-90"
SET_2 = "set 2"
SET_3 = "set 3"
SET_COLUMN_PREFIX = "set for "


def _read(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def compute_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _name_most_frequent(groups, classes, order):
    """Return per group the class with the most samples; the earlier class of `order` wins a tie."""
    return pd.crosstab(groups, classes).reindex(columns=list(order), fill_value=0).idxmax(axis=1)


def name_strata(samples, group_of_study):
    """Return the stratum key of every laboratory group: `<inducer class>|<cell type class>|<quiescent yes or no>`."""
    group = samples.study.map(group_of_study)
    positive = samples.state == SENESCENT
    inducer = _name_most_frequent(group[positive], samples.inducer_class[positive], INDUCER_CLASSES)
    cell_type = _name_most_frequent(group, samples.cell_type_class, CELL_TYPE_CLASSES)
    quiescent = ((samples.state == NOT_SENESCENT) & (samples.control_kind == QUIESCENT)).groupby(group).any()
    return (inducer + "|" + cell_type + "|" + quiescent.map({True: "quiescent yes", False: "quiescent no"})).rename("stratum")


def assign_parts(stratum_of_group, first_part, seed):
    """Return the part of every laboratory group, and the part that the next group would get.

    The strata are taken in the order of their key. Within a stratum the groups are sorted by number and reordered
    with one call of the generator's `permutation`; they get the parts in turn, and the turn carries over.
    """
    generator = np.random.default_rng(seed)
    turn = PARTS.index(first_part)
    parts = {}
    for stratum in sorted(stratum_of_group.unique()):
        groups = np.sort(stratum_of_group.index[stratum_of_group == stratum].to_numpy())
        for group in generator.permutation(groups):
            parts[group] = PARTS[turn]
            turn = 1 - turn
    return pd.Series(parts, dtype=str), PARTS[turn]


def read_earlier_parts(studies, earlier_split, amended_parts):
    """Return the part of every laboratory group that holds a study of the earlier split, and the earlier studies' parts.

    A study is found in the earlier split through its datasets, and keeps its part, unless `amended_parts` (study ->
    part) gives it another. A study left out of both parts gives its laboratory group no part.
    """
    earlier = earlier_split.assign(dataset=earlier_split.datasets.str.split(DATASET_SEPARATOR)).explode("dataset")
    now = studies.assign(dataset=studies.datasets.str.split(DATASET_SEPARATOR)).explode("dataset")
    now = now.merge(earlier[["dataset", "part"]], on="dataset")
    now["part"] = now.study.map(amended_parts).fillna(now.part)
    in_parts = now[now.part != LEFT_OUT_OF_BOTH_PARTS]
    for unit, rows in (("study", now), ("laboratory_group", in_parts)):
        parts = rows.groupby(unit).part.nunique()
        if (parts > 1).any():
            raise ValueError(f"the earlier split gives more than one part to: {unit} {list(parts.index[parts > 1])}")
    return (in_parts.drop_duplicates("laboratory_group").set_index("laboratory_group").part,
            now.drop_duplicates("study").set_index("study").part)


def validate_amendments(amendments, studies):
    """Refuse an amendment that moves a study into `test`, gives a part other than AMENDED_PARTS, names a study twice,
    or names a study that is not in set 1."""
    into_test = list(amendments.study[amendments.part == TEST])
    if into_test:
        raise ValueError(f"an amendment may not move a study into {TEST}: {into_test}")
    if not amendments.part.isin(AMENDED_PARTS).all() or amendments.study.duplicated().any():
        raise ValueError(f"an amendment gives a part other than {AMENDED_PARTS}, or names a study twice")
    in_set_1 = set(studies.study[studies.in_set_1 == "yes"])
    if not set(amendments.study) <= in_set_1:
        raise ValueError(f"amended studies not in set 1: {sorted(set(amendments.study) - in_set_1)}")


def validate_parts(studies, parts):
    """Refuse a split in which one laboratory group has studies in both parts."""
    in_parts = studies.assign(part=parts)[parts.isin(PARTS)]
    both = in_parts.groupby("laboratory_group").part.nunique()
    if (both > 1).any():
        raise ValueError(f"laboratory groups with studies in both parts: {list(both.index[both > 1])}")


def describe_studies(samples, series_table, matches, methods):
    """Return one row per study: datasets, laboratory group, institutes, name hashes, stratum, and the set per method."""
    in_table = series_table[series_table.study != ""]
    if (in_table.groupby("study").laboratory_group.nunique() > 1).any():
        raise ValueError("a study has datasets in more than one laboratory group")
    counted = matches[matches.role.isin(COUNTED_ROLES)]
    by_study = in_table.groupby("study")

    def join_distinct(column):
        return by_study[column].agg(lambda values: "; ".join(sorted(set(values) - {""})))

    studies = pd.DataFrame({
        "datasets": by_study.series.agg(lambda values: DATASET_SEPARATOR.join(sorted(values))),
        "laboratory_group": by_study.laboratory_group.first(),
        "institutes": join_distinct("contact_institute"),
        "contact_name_sha256": join_distinct("contact_name_sha256"),
        "last_contributor_sha256": join_distinct("last_contributor_sha256"),
    })
    number = series_table.series.map(rank_accessions(series_table.series))
    smallest = series_table.assign(number=number).sort_values("number").drop_duplicates(
        "laboratory_group").set_index("laboratory_group").series
    studies.insert(2, "laboratory_group_smallest_series", studies.laboratory_group.map(smallest))
    studies["stratum"] = studies.laboratory_group.map(name_strata(samples, studies.laboratory_group))
    studies["in_set_1"] = (~studies.laboratory_group.isin(counted.laboratory_group)).map({True: "yes", False: "no"})
    listed_studies = set(in_table.study[in_table.series.isin(counted.series)])
    studies["study_itself_in_a_list"] = studies.index.isin(listed_studies)
    studies["study_itself_in_a_list"] = studies.study_itself_in_a_list.map({True: "yes", False: "no"})
    for method, sets in describe_sets(matches, methods).items():
        studies[SET_COLUMN_PREFIX + method] = studies.laboratory_group.map(sets).fillna(SET_2)
    return studies.reset_index()


def describe_sets(matches, methods):
    """Return per method the set of every laboratory group that a list names, with the roles: `set 3: training, test`.

    A group in a method's list with a counted role is in the method's set 3; any other group is in its set 2, with
    `: applied only` where the method's authors only applied their method to a series of the group.
    """
    order = {role: position for position, role in enumerate((*COUNTED_ROLES, APPLIED_ONLY))}
    sets = {}
    for method in methods:
        of_method = matches[matches.method == method].drop_duplicates(["laboratory_group", "role"])
        of_method = of_method.sort_values("role", key=lambda role: role.map(order))
        roles = of_method.groupby("laboratory_group").role.agg(list)
        sets[method] = roles.map(lambda names: (SET_3 if set(names) & set(COUNTED_ROLES) else SET_2) + ": "
                                 + ", ".join(names))
    return sets


def read_set(value):
    """Return `set 2` or `set 3` from a value of a column `set for <method>`."""
    return value.split(":")[0]


def split_studies(studies, earlier_split, first_part, seed, amendments):
    """Return the part of every study, and the part that the next new laboratory group of set 1 would get.

    `amendments` (study, part) give studies a part of AMENDED_PARTS on top of the split continued from the earlier one.
    """
    validate_amendments(amendments, studies)
    amended_parts = amendments.set_index("study").part
    if earlier_split is None:
        part_of_earlier_group = part_of_earlier_study = pd.Series(dtype=str)
    else:
        part_of_earlier_group, part_of_earlier_study = read_earlier_parts(studies, earlier_split, amended_parts)
    groups = studies.drop_duplicates("laboratory_group").set_index("laboratory_group")
    new_in_set_1 = groups[(groups.in_set_1 == "yes") & ~groups.index.isin(part_of_earlier_group.index)]
    part_of_new_group, next_part = assign_parts(new_in_set_1.stratum, first_part, seed)
    part_of_group = pd.concat([part_of_earlier_group, part_of_new_group])
    part = studies.laboratory_group.map(part_of_group).fillna(NOT_IN_SET_1)
    # No study changes its part, also when its laboratory group has changed since the earlier split.
    earlier_part = studies.study.map(part_of_earlier_study)
    part = studies.study.map(amended_parts).fillna(earlier_part.fillna(part))
    validate_parts(studies, part)
    return part, next_part


def summarise_parts(studies, samples):
    """Return per part the numbers of laboratory groups, studies, datasets and samples, and the largest group's share."""
    samples = samples.merge(studies[["study", "part"]], on="study", validate="m:1")
    rows = []
    for part, rows_of_part in [*studies.groupby("part"), ("set 1", studies[studies.in_set_1 == "yes"])]:
        of_part = samples[samples.study.isin(rows_of_part.study)]
        largest = rows_of_part.laboratory_group.value_counts().iloc[0] if len(rows_of_part) else 0
        rows.append({"part": part, "laboratory_groups": rows_of_part.laboratory_group.nunique(),
                     "studies": len(rows_of_part), "datasets": of_part.dataset.nunique(), "samples": len(of_part),
                     "positives": int((of_part.state == SENESCENT).sum()),
                     "negatives": int((of_part.state == NOT_SENESCENT).sum()),
                     "studies_in_the_largest_laboratory_group": largest,
                     "share_of_studies_in_the_largest_laboratory_group":
                         round(largest / len(rows_of_part), 4) if len(rows_of_part) else ""})
    return pd.DataFrame(rows)


def list_sencat_candidates(studies, samples):
    """Return the studies of set 1 with IMR-90 cells and quiescent controls, candidates for SenCat's unnamed dataset."""
    imr90 = [names_key(line, SENCAT_CELL_LINE) or names_key(kind, SENCAT_CELL_LINE)
             for line, kind in zip(samples.cell_line, samples.cell_type)]
    quiescent = (samples.state == NOT_SENESCENT) & (samples.control_kind == QUIESCENT)
    candidates = set(samples.study[imr90]) & set(samples.study[quiescent])
    return studies.loc[studies.study.isin(candidates) & (studies.in_set_1 == "yes"),
                       ["study", "datasets", "laboratory_group", "part"]]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=Path, required=True, help="benchmark_bulk_samples.tsv of the round")
    parser.add_argument("--contacts", type=Path, required=True, help="geo_series_contacts.tsv (hashes of names)")
    parser.add_argument("--lists", type=Path, required=True, help="method_source_datasets.tsv")
    parser.add_argument("--main-scores", type=Path, required=True,
                        help="main_scores.tsv; every `list_method` in it gets a column `set for <method>`")
    parser.add_argument("--accession-pubmed-ids", type=Path, required=True,
                        help="table with the columns accession and pmid_or_doi, for listed accessions outside GEO")
    parser.add_argument("--laboratory-group-links", type=Path, required=True,
                        help="laboratory_group_links.tsv: links made by hand (series_a, series_b, evidence, date)")
    parser.add_argument("--split-amendments", type=Path, required=True,
                        help="split_amendments.tsv: amendments of parts (study, part, reason, date)")
    parser.add_argument("--earlier-split-dir", type=Path,
                        help="folder of the split to continue: its studies keep their parts and the turn continues")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.out_dir.exists():
        raise SystemExit(f"{args.out_dir} exists")

    samples = add_classes(_read(args.samples))
    lists = _read(args.lists)
    validate_roles(lists)
    outside_geo = _read(args.accession_pubmed_ids).rename(columns={"pmid_or_doi": "pubmed_id"})
    outside_geo = outside_geo.loc[outside_geo.pubmed_id.str.fullmatch(r"\d+"), ["accession", "pubmed_id"]]
    study_of_dataset = samples.drop_duplicates("dataset").set_index("dataset").study
    series_table, matches = build_laboratory_groups(study_of_dataset, args.contacts, lists, outside_geo,
                                                    _read(args.laboratory_group_links))
    methods = sorted(set(lists.method) | set(_read(args.main_scores).list_method))
    studies = describe_studies(samples, series_table, matches, methods)

    earlier_split, first_part, earlier_record = None, TEST, None
    if args.earlier_split_dir:
        earlier_split = _read(args.earlier_split_dir / SPLIT_NAME)
        earlier_record = json.loads((args.earlier_split_dir / RECORD_NAME).read_text())
        first_part = earlier_record["next_part"]
    part, next_part = split_studies(studies, earlier_split, first_part, SEED, _read(args.split_amendments))
    studies.insert(studies.columns.get_loc("stratum") + 1, "part", part)
    studies["method_source_datasets_sha256"] = compute_sha256(args.lists)
    studies["sample_table_sha256"] = compute_sha256(args.samples)
    lost = int(((studies.in_set_1 == "no") & (studies.study_itself_in_a_list == "no")).sum())
    record = {
        "written": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "seed": SEED, "generator": "numpy.random.default_rng", "numpy": np.__version__, "pandas": pd.__version__,
        "python": platform.python_version(), "command": sys.argv,
        "first_part": first_part, "next_part": next_part,
        "earlier_split": str(args.earlier_split_dir) if args.earlier_split_dir else None,
        "earlier_split_sha256": compute_sha256(args.earlier_split_dir / SPLIT_NAME) if args.earlier_split_dir else None,
        "samples": str(args.samples), "sample_table_sha256": compute_sha256(args.samples),
        "lists": str(args.lists), "method_source_datasets_sha256": compute_sha256(args.lists),
        "contacts": str(args.contacts), "contacts_sha256": compute_sha256(args.contacts),
        "laboratory_group_links": str(args.laboratory_group_links),
        "laboratory_group_links_sha256": compute_sha256(args.laboratory_group_links),
        "split_amendments": str(args.split_amendments),
        "split_amendments_sha256": compute_sha256(args.split_amendments),
        "studies_out_of_set_1_through_the_laboratory_rule_only": lost,
    }
    summary = summarise_parts(studies, samples)
    strata = studies[studies.in_set_1 == "yes"].groupby(["stratum", "part"]).agg(
        laboratory_groups=("laboratory_group", "nunique"), studies=("study", "size")).reset_index()

    args.out_dir.mkdir(parents=True)
    studies.to_csv(args.out_dir / SPLIT_NAME, sep="\t", index=False)
    series_table.to_csv(args.out_dir / SERIES_NAME, sep="\t", index=False)
    matches.sort_values(["laboratory_group", "series", "method"]).to_csv(args.out_dir / MATCHES_NAME, sep="\t", index=False)
    summary.to_csv(args.out_dir / SUMMARY_NAME, sep="\t", index=False)
    strata.to_csv(args.out_dir / STRATA_NAME, sep="\t", index=False)
    list_sencat_candidates(studies, samples).to_csv(args.out_dir / SENCAT_CANDIDATES_NAME, sep="\t", index=False)
    record["study_split_sha256"] = compute_sha256(args.out_dir / SPLIT_NAME)
    (args.out_dir / RECORD_NAME).write_text(json.dumps(record, indent=2) + "\n")
    print(summary.to_string(index=False))
    print(f"{lost} studies are out of set 1 through the laboratory rule only; numpy {np.__version__}")


if __name__ == "__main__":
    main()
