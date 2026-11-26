"""Laboratory groups of GEO series, and which laboratory groups are in the list of a method.

Two series are connected when they share the contact person with the institute, the last contributor name, two
contributor names, a PubMed id, a BioProject or a SuperSeries, when they are datasets of one study, or when a link made
by hand names them. A series of the benchmark table also takes these keys from the rows of its SuperSeries and
SubSeries in the contact table. A laboratory group is a connected set.
"""

import pandas as pd

from build_method_source_datasets import COUNTED_ROLES
from harvest_geo_contacts import GEO_SERIES, normalise_name
from harvest_geo_samples import SEPARATOR
from mark_tool_training_overlap import read_related_accessions

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


def list_groups_in_lists(groups, matches):
    """Return one row per laboratory group, method and role in which a series of the group matches the method's list."""
    in_lists = matches.assign(laboratory_group=matches.series.map(groups))
    return in_lists[["laboratory_group", "method", "role"]].drop_duplicates()


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
