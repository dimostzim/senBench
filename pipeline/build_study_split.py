#!/usr/bin/env python3
"""Build the laboratory groups, the sets and the split of the benchmark across studies into `test` and `development`,
and the tables of the classes of cell type and of time after induction, for a person to check.

Two series are connected when they share the contact person with the institute, the last contributor name, two
contributor names, a PubMed id, a BioProject or a SuperSeries, when they are datasets of one study, or when a link made
by hand names them. A series of the benchmark table also takes these keys from the rows of its SuperSeries and SubSeries
in the contact table. A laboratory group is a connected set. Set 1 holds the laboratory groups of which no series
matches the list of a method; they are divided into `test` and `development` within strata of the class of the
inducer, the class of the cell type and the presence of quiescent controls.
"""

import argparse
import hashlib
import json
import platform
import re
import string
import sys
import unicodedata
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

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

CELL_TYPE_CLASSES_NAME = "cell_type_classes.tsv"
TIME_CLASSES_NAME = "time_classes.tsv"


# The roles of the lists of the methods.

DERIVED_FROM = "derived from"
APPLIED_ONLY = "applied only"
# A laboratory group with a listed accession in one of these roles is in the list; `applied only` is reported only.
COUNTED_ROLES = ("training", "validation or model selection", "test", "unknown role", DERIVED_FROM)


def validate_roles(lists):
    unknown = sorted(set(lists.role) - set(COUNTED_ROLES) - {APPLIED_ONLY})
    if unknown:
        raise ValueError(f"roles outside the allowed list: {unknown}")


# The classes of the inducer, the cell type and the time after induction.

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


# The laboratory groups.

GEO_SERIES = re.compile(r"GSE\d+")
# Separates several values in a cell of the contact table.
SEPARATOR = "; "
_PUNCTUATION = str.maketrans(string.punctuation, " " * len(string.punctuation))
# Records write a German umlaut either as the letter or as two letters (Göttingen, Goettingen).
_UMLAUTS = str.maketrans({"ä": "ae", "ö": "oe", "ü": "ue", "ß": "ss"})


def normalise_name(name):
    """Return a name in lower case without punctuation and accents, with German umlauts written as two letters."""
    letters = unicodedata.normalize("NFKD", unicodedata.normalize("NFC", name.lower()).translate(_UMLAUTS))
    unaccented = "".join(letter for letter in letters if not unicodedata.combining(letter))
    return " ".join(unaccented.translate(_PUNCTUATION).split())


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


SUPERSERIES = "SuperSeries"
SUBSERIES = "SubSeries"
BIOPROJECT = "BioProject"
PUBMED = "PubMed id"
CONTACT = "contact person with the institute"
LAST_CONTRIBUTOR = "last contributor"
TWO_CONTRIBUTORS = "two contributors"
LINK_BY_HAND = "link made by hand"
# The keys that a series of the table takes from the rows of its SuperSeries and SubSeries.
KEYS_FROM_RELATIVES = (CONTACT, LAST_CONTRIBUTOR, TWO_CONTRIBUTORS, PUBMED, BIOPROJECT)


def _split_values(column):
    """Return one row per value of a column whose cells hold several values; empty cells give no row."""
    values = column.str.split(SEPARATOR).explode()
    return values[values != ""]


def list_pubmed_ids(contacts):
    """Return one row per series and PubMed id: series, pubmed_id."""
    return _split_values(contacts.set_index("series").pubmed_ids).rename("pubmed_id").reset_index()


def list_contributor_pairs(contacts):
    """Return one row per series and pair of its contributors: series, kind, key (the two name hashes in sorted
    order)."""
    names = _split_values(contacts.set_index("series").contributors_sha256).rename("name").reset_index()
    names = names.drop_duplicates()
    pairs = names.merge(names, on="series")
    pairs = pairs[pairs.name_x < pairs.name_y]
    return pd.DataFrame({"series": pairs.series, "kind": TWO_CONTRIBUTORS, "key": pairs.name_x + " " + pairs.name_y})


def list_link_keys(contacts, related, study_of_dataset):
    """Return one row per series and key that connects it to other series: series, kind, key."""
    has_contact = (contacts.contact_name_sha256 != "") & (contacts.contact_institute != "")
    contact = contacts[has_contact]
    last = contacts[contacts.last_contributor_sha256 != ""]
    pubmed = list_pubmed_ids(contacts)
    linking = related[related.kind.isin([SUPERSERIES, BIOPROJECT])]
    # A SuperSeries shares the key of its SubSeries: its own accession.
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
    """Return the keys of KEYS_FROM_RELATIVES of the SuperSeries and SubSeries rows of the contact table, as keys of
    the series of the table that they are related to (in either row's relations): series, kind, key."""
    pairs = related[related.kind.isin([SUPERSERIES, SUBSERIES])][["series", "related"]]
    pairs = pd.concat([pairs, pairs.rename(columns={"series": "related", "related": "series"})]).drop_duplicates()
    pairs = pairs[pairs.series.isin(table_series) & pairs.related.isin(contacts.series)]
    relatives = contacts[contacts.series.isin(pairs.related)]
    keys = list_link_keys(relatives, related[related.series.isin(relatives.series)], pd.Series(dtype=str))
    keys = keys[keys.kind.isin(KEYS_FROM_RELATIVES)].rename(columns={"series": "related"})
    return pairs.merge(keys, on="related")[["series", "kind", "key"]].drop_duplicates()


def list_links_by_hand(series, links):
    """Return the keys of the links made by hand (series_a, series_b): series, kind, key. Both series of a link must
    take part in the grouping."""
    unknown = sorted((set(links.series_a) | set(links.series_b)) - set(series))
    if unknown:
        raise ValueError(f"links made by hand name series outside the grouping: {unknown}")
    key = links.series_a + " " + links.series_b
    return pd.DataFrame({"series": pd.concat([links.series_a, links.series_b], ignore_index=True),
                         "kind": LINK_BY_HAND, "key": pd.concat([key, key], ignore_index=True)})


def rank_accessions(accessions):
    """Return the rank of every accession: GEO series in the order of their number, then other accessions by text."""
    def order(accession):
        return (0, int(accession[3:]), "") if GEO_SERIES.fullmatch(accession) else (1, 0, accession)

    ordered = sorted(set(accessions), key=order)
    return pd.Series(range(len(ordered)), index=ordered)


def find_groups(series, link_keys):
    """Return the laboratory group number of every series; groups are numbered by their smallest accession in the
    order of rank_accessions()."""
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
    """Return one row per series and accession that names it or is related to it: series, accession."""
    own = pd.DataFrame({"series": series, "accession": series})
    return pd.concat([own, related.rename(columns={"related": "accession"})[["series", "accession"]]]).drop_duplicates()


def match_series_to_lists(series, related, pubmed_ids, lists, pubmed_ids_outside_geo):
    """Return one row per series, list row and way of matching: series, method, role, accession, matched_through.

    A series matches a listed accession when it is that accession, or shares a SuperSeries, SubSeries, BioProject, SRA
    study or PubMed id with it. `pubmed_ids` has the PubMed ids of GEO series (series, pubmed_id), and
    `pubmed_ids_outside_geo` those of listed accessions outside GEO (accession, pubmed_id).
    """
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
    """Return per key of the index the text `method (role, role); method (role)` for the memberships with these roles."""
    kept = memberships.loc[memberships.role.isin(roles), ["method", "role"]].rename_axis("key").reset_index()
    order = {role: position for position, role in enumerate(roles)}
    kept = kept.drop_duplicates().sort_values("role", key=lambda role: role.map(order))
    per_method = kept.groupby(["key", "method"]).role.agg(", ".join).reset_index()
    return (per_method.method + " (" + per_method.role + ")").groupby(per_method.key).agg("; ".join)


def build_laboratory_groups(study_of_dataset, contacts_path, lists, pubmed_ids_outside_geo, links_by_hand):
    """Return the table of series with their laboratory group, and the matches of the series to the lists.

    `study_of_dataset` maps every dataset of the benchmark table to its study. Listed GEO series take part too.
    `links_by_hand` has the columns series_a and series_b.
    """
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


# The split.

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


def split_studies(studies, seed, amendments):
    """Return the part of every study, and the part that the next laboratory group of set 1 would get.

    `amendments` (study, part) give studies a part of AMENDED_PARTS on top of the split.
    """
    validate_amendments(amendments, studies)
    groups = studies.drop_duplicates("laboratory_group").set_index("laboratory_group")
    in_set_1 = groups[groups.in_set_1 == "yes"]
    part_of_group, next_part = assign_parts(in_set_1.stratum, TEST, seed)
    part = studies.laboratory_group.map(part_of_group).fillna(NOT_IN_SET_1)
    part = studies.study.map(amendments.set_index("study").part).fillna(part)
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
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--classes-dir", type=Path, required=True,
                        help=f"folder for {CELL_TYPE_CLASSES_NAME} and {TIME_CLASSES_NAME}: the class of every pair of "
                             "cell line and cell type and of every time after induction of the samples")
    args = parser.parse_args()
    class_tables = (args.classes_dir / CELL_TYPE_CLASSES_NAME, args.classes_dir / TIME_CLASSES_NAME)
    for path in (args.out_dir, *class_tables):
        if path.exists():
            raise SystemExit(f"{path} exists")

    labelled = _read(args.samples)
    samples = add_classes(labelled)
    lists = _read(args.lists)
    validate_roles(lists)
    outside_geo = _read(args.accession_pubmed_ids).rename(columns={"pmid_or_doi": "pubmed_id"})
    outside_geo = outside_geo.loc[outside_geo.pubmed_id.str.fullmatch(r"\d+"), ["accession", "pubmed_id"]]
    study_of_dataset = samples.drop_duplicates("dataset").set_index("dataset").study
    series_table, matches = build_laboratory_groups(study_of_dataset, args.contacts, lists, outside_geo,
                                                    _read(args.laboratory_group_links))
    methods = sorted(set(lists.method) | set(_read(args.main_scores).list_method))
    studies = describe_studies(samples, series_table, matches, methods)

    part, next_part = split_studies(studies, SEED, _read(args.split_amendments))
    studies.insert(studies.columns.get_loc("stratum") + 1, "part", part)
    studies["method_source_datasets_sha256"] = compute_sha256(args.lists)
    studies["sample_table_sha256"] = compute_sha256(args.samples)
    lost = int(((studies.in_set_1 == "no") & (studies.study_itself_in_a_list == "no")).sum())
    record = {
        "written": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "seed": SEED, "generator": "numpy.random.default_rng", "numpy": np.__version__, "pandas": pd.__version__,
        "python": platform.python_version(), "command": sys.argv,
        "first_part": TEST, "next_part": next_part,
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
    cell_types = tabulate_cell_type_classes(labelled)
    times = tabulate_time_classes(labelled)
    args.classes_dir.mkdir(parents=True, exist_ok=True)
    cell_types.to_csv(class_tables[0], sep="\t", index=False)
    times.to_csv(class_tables[1], sep="\t", index=False)
    print(summary.to_string(index=False))
    print(f"{lost} studies are out of set 1 through the laboratory rule only; numpy {np.__version__}")
    other = cell_types[cell_types.cell_type_class == OTHER_CELL_TYPE].samples.sum()
    no_time = times[times.time_class == NO_TIME].positives.sum()
    print(f"{other} of {len(labelled)} samples in `{OTHER_CELL_TYPE}`; "
          f"{no_time} of {times.positives.sum()} positives in `{NO_TIME}`")


if __name__ == "__main__":
    main()
