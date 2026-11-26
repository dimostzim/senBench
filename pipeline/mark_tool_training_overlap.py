#!/usr/bin/env python3
"""Mark which datasets with both states are in the dataset lists of the published tools, by role."""

import argparse
import re
from pathlib import Path

import pandas as pd

UNKNOWN_ROLE = "unknown role"
ROLE_COLUMNS = {
    "training": "tools_training",
    "validation or model selection": "tools_training",
    "test": "tools_test",
    "applied only": "tools_applied_only",
}
RELATION = re.compile(r"(SubSeries of|SuperSeries of|BioProject|SRA): ([^\s;]+)")
RELATED_KINDS = {"SubSeries of": "SuperSeries", "SuperSeries of": "SubSeries", "BioProject": "BioProject",
                 "SRA": "SRA study"}


def read_related_accessions(series_path):
    """Return one row per GEO series and related accession: series, kind (of the related accession), related."""
    series = pd.read_csv(series_path, sep="\t", dtype=str, keep_default_na=False, usecols=["series", "relations"])
    found = series.relations.str.extractall(RELATION)
    return pd.DataFrame({
        "series": series.series.loc[found.index.get_level_values(0)].values,
        "kind": found[0].map(RELATED_KINDS).values,
        "related": found[1].str.split(r"[/=]").str[-1].values,  # BioProject and SRA are given as addresses
    }).drop_duplicates()


def list_matching_accessions(datasets, related):
    """Return one row per dataset and accession that names it: dataset, accession, relation, counted."""
    own = pd.DataFrame({"dataset": datasets, "accession": datasets, "relation": "", "counted": True})
    direct = related[related.series.isin(datasets)]
    direct = pd.DataFrame({"dataset": direct.series, "accession": direct.related,
                           "relation": "its " + direct.kind + " " + direct.related, "counted": True})
    # The table states a SuperSeries from both sides; SubSeries of one SuperSeries are reported and not counted.
    below = related[related.kind == "SuperSeries"].rename(columns={"series": "sub", "related": "super"})
    above = related[related.kind == "SubSeries"].rename(columns={"series": "super", "related": "sub"})
    pairs = pd.concat([below, above])[["super", "sub"]].drop_duplicates()
    siblings = pairs[pairs["sub"].isin(datasets)].merge(pairs, on="super", suffixes=("", "_other"))
    siblings = siblings[siblings["sub"] != siblings.sub_other]
    siblings = pd.DataFrame({
        "dataset": siblings["sub"], "accession": siblings.sub_other,
        "relation": siblings.sub_other + ", a SubSeries of the same SuperSeries " + siblings["super"] + ",",
        "counted": False})
    return pd.concat([own, direct, siblings], ignore_index=True)


def join_unique(values, separator=";"):
    return separator.join(dict.fromkeys(values))


def write_match_notes(matches):
    """Return one note per dataset: the matches that are not a counted match by the dataset's own accession."""
    noted = matches[~matches.counted | (matches.relation != "") | (matches.role == UNKNOWN_ROLE)].copy()
    relation = noted.relation.where(noted.relation != "", "its own accession")
    noted["start"] = ("listed as " + relation + " by").where(
        noted.counted, "not counted: " + relation + " is listed by")
    named_role = ~noted.counted | (noted.role == UNKNOWN_ROLE)
    noted["tool_named"] = noted.tool.where(~named_role, noted.tool + " (" + noted.role + ")")
    by_relation = noted.groupby(["dataset", "start"], sort=False).tool_named.agg(join_unique, separator=", ")
    by_relation = (by_relation.index.get_level_values("start") + " " + by_relation).rename("note")
    return by_relation.groupby("dataset").agg(join_unique, separator="; ")


def mark_overlap(datasets, tool_datasets, related):
    """Return one row per dataset with the tools that list it, by role, and whether no tool lists it."""
    unknown = set(tool_datasets.role) - set(ROLE_COLUMNS) - {UNKNOWN_ROLE}
    if unknown:
        raise ValueError(f"roles outside the allowed list: {sorted(unknown)}")
    matches = list_matching_accessions(datasets, related).merge(tool_datasets, on="accession")
    table = pd.DataFrame({"dataset": datasets}).set_index("dataset")
    counted = matches[matches.counted]
    for column in dict.fromkeys(ROLE_COLUMNS.values()):
        roles = [role for role, role_column in ROLE_COLUMNS.items() if role_column == column]
        table[column] = counted[counted.role.isin(roles)].groupby("dataset").tool.agg(join_unique)
    table["in_no_list"] = ["no" if dataset in set(counted.dataset) else "yes" for dataset in table.index]
    table["match_note"] = write_match_notes(matches)
    return table.fillna("").reset_index()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--datasets", type=Path, required=True)
    parser.add_argument("--yes-column", default="both_states_by_rules",
                        help="column of --datasets; the datasets with `yes` in it are marked")
    parser.add_argument("--tool-datasets", type=Path, required=True)
    parser.add_argument("--geo-series", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    datasets = pd.read_csv(args.datasets, sep="\t", dtype=str, keep_default_na=False)
    datasets = datasets[datasets[args.yes_column] == "yes"].dataset.reset_index(drop=True)
    tool_datasets = pd.read_csv(args.tool_datasets, sep="\t", dtype=str, keep_default_na=False)
    table = mark_overlap(datasets, tool_datasets, read_related_accessions(args.geo_series))
    table.to_csv(args.out, sep="\t", index=False)
    print(f"{len(table)} datasets, {(table.in_no_list == 'yes').sum()} in no list, "
          f"{(table.tools_training != '').sum()} listed for training or validation by at least one tool")


if __name__ == "__main__":
    main()
