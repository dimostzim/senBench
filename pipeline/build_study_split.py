#!/usr/bin/env python3
import argparse
import hashlib
import re
import string
import unicodedata
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 20261116
TEST = "test"
DEVELOPMENT = "development"
NOT_IN_SET_1 = "not in set 1"
PARTS = (TEST, DEVELOPMENT)
DATASET_SEPARATOR = ";"
SPLIT_NAME = "study_split.tsv"
SERIES_NAME = "laboratory_group_series.tsv"
MATCHES_NAME = "list_matches.tsv"
SUMMARY_NAME = "split_summary.tsv"
STRATA_NAME = "split_strata.tsv"
SET_2 = "set 2"
SET_3 = "set 3"
SET_COLUMN_PREFIX = "set for "


DERIVED_FROM = "derived from"
APPLIED_ONLY = "applied only"
COUNTED_ROLES = ("training", "validation or model selection", "test", "unknown role", DERIVED_FROM)


def validate_roles(lists):
    unknown = sorted(set(lists.role) - set(COUNTED_ROLES) - {APPLIED_ONLY})
    if unknown:
        raise ValueError(f"roles outside the allowed list: {unknown}")


SENESCENT = "senescent"
NOT_SENESCENT = "not senescent"
QUIESCENT = "quiescent"

OTHER_INDUCER = "other"
INDUCER_CLASSES = ("replicative", "irradiation", "oncogene", "DNA-damaging drug", "CDK4/6 inhibitor", "oxidative stress",
                   OTHER_INDUCER)
INDUCERS_IN_OTHER = ("therapy", "other", "not stated", "none")

OTHER_CELL_TYPE = "other or not stated"
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
TUMOUR_WORD_ENDINGS = ("carcinoma", "sarcoma", "blastoma", "cytoma", "glioma", "melanoma", "lymphoma", "leukemia",
                       "leukaemia", "mesothelioma")
CANCER_CELL_LINE = "cancer cell line"
CELL_TYPE_CLASSES = (*CELL_TYPE_KEYS, OTHER_CELL_TYPE)

def _whole_word(key):
    return re.compile(r"(?<![a-z0-9])" + re.sub(r"[- ]", "[- ]?", key) + r"(?![a-z0-9])")


_CELL_TYPE_PATTERNS = {name: [_whole_word(key) for key in keys] for name, keys in CELL_TYPE_KEYS.items()}
_CELL_TYPE_PATTERNS[CANCER_CELL_LINE].append(re.compile(r"[a-z]*(?:" + "|".join(TUMOUR_WORD_ENDINGS) + r")(?![a-z0-9])"))


def classify_inducer(inducer):
    if inducer in INDUCERS_IN_OTHER:
        return OTHER_INDUCER
    if inducer not in INDUCER_CLASSES:
        raise ValueError(f"inducer outside the allowed list: {inducer!r}")
    return inducer


def classify_cell_type(cell_line, cell_type):
    for text in (cell_line.lower(), cell_type.lower()):
        for name, patterns in _CELL_TYPE_PATTERNS.items():
            if any(pattern.search(text) for pattern in patterns):
                return name
    return OTHER_CELL_TYPE


def add_classes(samples):
    positive = samples.state == SENESCENT
    cell_type_class = [classify_cell_type(line, kind) for line, kind in zip(samples.cell_line, samples.cell_type)]
    return samples.assign(
        inducer_class=samples.inducer.map(classify_inducer).where(positive, ""),
        cell_type_class=cell_type_class)


GEO_SERIES = re.compile(r"GSE\d+")
SEPARATOR = "; "
_PUNCTUATION = str.maketrans(string.punctuation, " " * len(string.punctuation))
_UMLAUTS = str.maketrans({"ä": "ae", "ö": "oe", "ü": "ue", "ß": "ss"})


def normalise_name(name):
    letters = unicodedata.normalize("NFKD", unicodedata.normalize("NFC", name.lower()).translate(_UMLAUTS))
    unaccented = "".join(letter for letter in letters if not unicodedata.combining(letter))
    return " ".join(unaccented.translate(_PUNCTUATION).split())


RELATION = re.compile(r"(SubSeries of|SuperSeries of|BioProject|SRA): ([^\s;]+)")
RELATED_KINDS = {"SubSeries of": "SuperSeries", "SuperSeries of": "SubSeries", "BioProject": "BioProject",
                 "SRA": "SRA study"}


def read_related_accessions(series_path):
    series = pd.read_csv(series_path, sep="\t", dtype=str, keep_default_na=False, usecols=["series", "relations"])
    found = series.relations.str.extractall(RELATION)
    return pd.DataFrame({
        "series": series.series.loc[found.index.get_level_values(0)].values,
        "kind": found[0].map(RELATED_KINDS).values,
        "related": found[1].str.split(r"[/=]").str[-1].values,
    }).drop_duplicates()


SUPERSERIES = "SuperSeries"
SUBSERIES = "SubSeries"
BIOPROJECT = "BioProject"
PUBMED = "PubMed id"
CONTACT = "contact person with the institute"
LAST_CONTRIBUTOR = "last contributor"
TWO_CONTRIBUTORS = "two contributors"
LINK_BY_HAND = "link made by hand"
KEYS_FROM_RELATIVES = (CONTACT, LAST_CONTRIBUTOR, TWO_CONTRIBUTORS, PUBMED, BIOPROJECT)


def _split_values(column):
    values = column.str.split(SEPARATOR).explode()
    return values[values != ""]


def list_pubmed_ids(contacts):
    return _split_values(contacts.set_index("series").pubmed_ids).rename("pubmed_id").reset_index()


def list_contributor_pairs(contacts):
    names = _split_values(contacts.set_index("series").contributors_sha256).rename("name").reset_index()
    names = names.drop_duplicates()
    pairs = names.merge(names, on="series")
    pairs = pairs[pairs.name_x < pairs.name_y]
    return pd.DataFrame({"series": pairs.series, "kind": TWO_CONTRIBUTORS, "key": pairs.name_x + " " + pairs.name_y})


def list_link_keys(contacts, related, study_of_dataset):
    has_contact = (contacts.contact_name_sha256 != "") & (contacts.contact_institute != "")
    contact = contacts[has_contact]
    last = contacts[contacts.last_contributor_sha256 != ""]
    pubmed = list_pubmed_ids(contacts)
    linking = related[related.kind.isin([SUPERSERIES, BIOPROJECT])]
    above = related[related.kind == SUBSERIES]
    keys = [
        pd.DataFrame({"series": contact.series, "kind": CONTACT,
                      "key": contact.contact_name_sha256 + " " + contact.contact_institute.map(normalise_name)}),
        pd.DataFrame({"series": last.series, "kind": LAST_CONTRIBUTOR, "key": last.last_contributor_sha256}),
        list_contributor_pairs(contacts),
        pd.DataFrame({"series": pubmed.series, "kind": PUBMED, "key": pubmed.pubmed_id}),
        pd.DataFrame({"series": linking.series, "kind": linking.kind, "key": linking.related}),
        pd.DataFrame({"series": above.series, "kind": SUPERSERIES, "key": above.series}),
        pd.DataFrame({"series": study_of_dataset.index, "kind": "study", "key": study_of_dataset.values}),
    ]
    return pd.concat(keys, ignore_index=True).drop_duplicates()


def list_keys_from_relatives(table_series, contacts, related):
    pairs = related[related.kind.isin([SUPERSERIES, SUBSERIES])][["series", "related"]]
    pairs = pd.concat([pairs, pairs.rename(columns={"series": "related", "related": "series"})]).drop_duplicates()
    pairs = pairs[pairs.series.isin(table_series) & pairs.related.isin(contacts.series)]
    relatives = contacts[contacts.series.isin(pairs.related)]
    keys = list_link_keys(relatives, related[related.series.isin(relatives.series)], pd.Series(dtype=str))
    keys = keys[keys.kind.isin(KEYS_FROM_RELATIVES)].rename(columns={"series": "related"})
    return pairs.merge(keys, on="related")[["series", "kind", "key"]].drop_duplicates()


def list_links_by_hand(series, links):
    unknown = sorted((set(links.series_a) | set(links.series_b)) - set(series))
    if unknown:
        raise ValueError(f"links made by hand name series outside the grouping: {unknown}")
    key = links.series_a + " " + links.series_b
    return pd.DataFrame({"series": pd.concat([links.series_a, links.series_b], ignore_index=True),
                         "kind": LINK_BY_HAND, "key": pd.concat([key, key], ignore_index=True)})


def rank_accessions(accessions):
    def order(accession):
        return (0, int(accession[3:]), "") if GEO_SERIES.fullmatch(accession) else (1, 0, accession)

    ordered = sorted(set(accessions), key=order)
    return pd.Series(range(len(ordered)), index=ordered)


def find_groups(series, link_keys):
    parent = {accession: accession for accession in series}

    def find(accession):
        while parent[accession] != accession:
            parent[accession] = parent[parent[accession]]
            accession = parent[accession]
        return accession

    for members in link_keys[link_keys.series.isin(parent)].groupby(["kind", "key"]).series.agg(list):
        for other in members[1:]:
            parent[find(other)] = find(members[0])
    root = pd.Series({accession: find(accession) for accession in series})
    number = root.index.to_series().map(rank_accessions(root.index))
    smallest = number.groupby(root).transform("min")
    return smallest.rank(method="dense").astype(int).rename("laboratory_group")


def list_accessions_of_series(series, related):
    own = pd.DataFrame({"series": series, "accession": series})
    return pd.concat([own, related.rename(columns={"related": "accession"})[["series", "accession"]]]).drop_duplicates()


def match_series_to_lists(series, related, pubmed_ids, lists, pubmed_ids_outside_geo):
    listed = lists[["method", "role", "accession"]].drop_duplicates()
    names = list_accessions_of_series(pd.concat([pd.Series(series), listed.accession]).drop_duplicates(), related)
    listed_names = listed.merge(names.rename(columns={"series": "accession", "accession": "through"}), on="accession")
    by_accession = names[names.series.isin(series)].rename(columns={"accession": "through"}).merge(
        listed_names, on="through")
    listed_pubmed = pd.concat([pubmed_ids.rename(columns={"series": "accession"}), pubmed_ids_outside_geo])
    listed_pubmed = listed.merge(listed_pubmed, on="accession")
    by_pubmed = pubmed_ids[pubmed_ids.series.isin(series)].merge(listed_pubmed, on="pubmed_id")
    by_pubmed = by_pubmed.assign(through=PUBMED + " " + by_pubmed.pubmed_id)
    columns = ["series", "method", "role", "accession", "through"]
    matches = pd.concat([by_accession[columns], by_pubmed[columns]], ignore_index=True)
    return matches.rename(columns={"through": "matched_through"}).drop_duplicates()


def describe_lists(memberships, roles):
    kept = memberships.loc[memberships.role.isin(roles), ["method", "role"]].rename_axis("key").reset_index()
    order = {role: position for position, role in enumerate(roles)}
    kept = kept.drop_duplicates().sort_values("role", key=lambda role: role.map(order))
    per_method = kept.groupby(["key", "method"]).role.agg(", ".join).reset_index()
    return (per_method.method + " (" + per_method.role + ")").groupby(per_method.key).agg("; ".join)


def build_laboratory_groups(study_of_dataset, contacts_path, lists, pubmed_ids_outside_geo, links_by_hand):
    all_contacts = pd.read_csv(contacts_path, sep="\t", dtype=str, keep_default_na=False)
    listed_series = lists.accession[lists.accession.str.fullmatch(GEO_SERIES)]
    series = pd.concat([pd.Series(study_of_dataset.index), listed_series]).drop_duplicates()
    not_harvested = sorted(set(series) - set(all_contacts.series))
    if not_harvested:
        raise ValueError(f"no contact fields for {len(not_harvested)} series, the first: {not_harvested[:5]}")
    contacts = all_contacts[all_contacts.series.isin(series)]
    all_related = read_related_accessions(contacts_path)
    related = all_related[all_related.series.isin(series)]
    link_keys = pd.concat([list_link_keys(contacts, related, study_of_dataset),
                           list_keys_from_relatives(study_of_dataset.index, all_contacts, all_related),
                           list_links_by_hand(series, links_by_hand)], ignore_index=True)
    groups = find_groups(list(series), link_keys)
    pubmed_ids = list_pubmed_ids(contacts)
    matches = match_series_to_lists(list(series), related, pubmed_ids, lists, pubmed_ids_outside_geo)
    matches.insert(0, "laboratory_group", matches.series.map(groups))

    table = contacts.drop(columns="relations").assign(laboratory_group=contacts.series.map(groups))
    table["study"] = table.series.map(study_of_dataset).fillna("")
    for kind, column in ((BIOPROJECT, "bioprojects"), (SUPERSERIES, "superseries"), (SUBSERIES, "subseries")):
        values = related[related.kind == kind].groupby("series").related.agg(SEPARATOR.join)
        table[column] = table.series.map(values).fillna("")
    listed_by = describe_lists(matches.set_index("series"), COUNTED_ROLES)
    table["matches_the_lists_of"] = table.series.map(listed_by).fillna("")
    return table.sort_values(["laboratory_group", "series"]), matches


def _read(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def compute_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _name_most_frequent(groups, classes, order):
    return pd.crosstab(groups, classes).reindex(columns=list(order), fill_value=0).idxmax(axis=1)


def name_strata(samples, group_of_study):
    group = samples.study.map(group_of_study)
    positive = samples.state == SENESCENT
    inducer = _name_most_frequent(group[positive], samples.inducer_class[positive], INDUCER_CLASSES)
    cell_type = _name_most_frequent(group, samples.cell_type_class, CELL_TYPE_CLASSES)
    quiescent = ((samples.state == NOT_SENESCENT) & (samples.control_kind == QUIESCENT)).groupby(group).any()
    return (inducer + "|" + cell_type + "|" + quiescent.map({True: "quiescent yes", False: "quiescent no"})).rename("stratum")


def assign_parts(stratum_of_group, first_part, seed):
    generator = np.random.default_rng(seed)
    turn = PARTS.index(first_part)
    parts = {}
    for stratum in sorted(stratum_of_group.unique()):
        groups = np.sort(stratum_of_group.index[stratum_of_group == stratum].to_numpy())
        for group in generator.permutation(groups):
            parts[group] = PARTS[turn]
            turn = 1 - turn
    return pd.Series(parts, dtype=str)


def validate_parts(studies, parts):
    in_parts = studies.assign(part=parts)[parts.isin(PARTS)]
    both = in_parts.groupby("laboratory_group").part.nunique()
    if (both > 1).any():
        raise ValueError(f"laboratory groups with studies in both parts: {list(both.index[both > 1])}")


def describe_studies(samples, series_table, matches, methods):
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
    order = {role: position for position, role in enumerate((*COUNTED_ROLES, APPLIED_ONLY))}
    sets = {}
    for method in methods:
        of_method = matches[matches.method == method].drop_duplicates(["laboratory_group", "role"])
        of_method = of_method.sort_values("role", key=lambda role: role.map(order))
        roles = of_method.groupby("laboratory_group").role.agg(list)
        sets[method] = roles.map(lambda names: (SET_3 if set(names) & set(COUNTED_ROLES) else SET_2) + ": "
                                 + ", ".join(names))
    return sets


def split_studies(studies, seed):
    groups = studies.drop_duplicates("laboratory_group").set_index("laboratory_group")
    in_set_1 = groups[groups.in_set_1 == "yes"]
    part_of_group = assign_parts(in_set_1.stratum, TEST, seed)
    part = studies.laboratory_group.map(part_of_group).fillna(NOT_IN_SET_1)
    validate_parts(studies, part)
    return part


def summarise_parts(studies, samples):
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=Path, required=True, help="benchmark_bulk_samples.tsv of the round")
    parser.add_argument("--contacts", type=Path, required=True, help="geo_series_contacts.tsv (hashes of names)")
    parser.add_argument("--lists", type=Path, required=True, help="method_source_datasets.tsv")
    parser.add_argument("--list-methods", type=Path, required=True,
                        help="list_methods.tsv; every `list_method` in it gets a column `set for <method>`")
    parser.add_argument("--accession-pubmed-ids", type=Path, required=True,
                        help="table with the columns accession and pmid_or_doi, for listed accessions outside GEO")
    parser.add_argument("--laboratory-group-links", type=Path, required=True,
                        help="laboratory_group_links.tsv: links made by hand (series_a, series_b)")
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
    methods = sorted(set(lists.method) | set(_read(args.list_methods).list_method))
    studies = describe_studies(samples, series_table, matches, methods)

    studies.insert(studies.columns.get_loc("stratum") + 1, "part", split_studies(studies, SEED))
    studies["method_source_datasets_sha256"] = compute_sha256(args.lists)
    studies["sample_table_sha256"] = compute_sha256(args.samples)
    summary = summarise_parts(studies, samples)
    strata = studies[studies.in_set_1 == "yes"].groupby(["stratum", "part"]).agg(
        laboratory_groups=("laboratory_group", "nunique"), studies=("study", "size")).reset_index()

    args.out_dir.mkdir(parents=True)
    studies.to_csv(args.out_dir / SPLIT_NAME, sep="\t", index=False)
    series_table.to_csv(args.out_dir / SERIES_NAME, sep="\t", index=False)
    matches.sort_values(["laboratory_group", "series", "method"]).to_csv(args.out_dir / MATCHES_NAME, sep="\t", index=False)
    summary.to_csv(args.out_dir / SUMMARY_NAME, sep="\t", index=False)
    strata.to_csv(args.out_dir / STRATA_NAME, sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
