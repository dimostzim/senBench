"""Check the assignment of laboratory groups to the parts `test` and `development`."""

import pandas as pd
import pytest

from build_study_split import DEVELOPMENT, LEFT_OUT_OF_BOTH_PARTS, SEED, TEST, assign_parts, split_studies
from laboratory_groups import build_laboratory_groups, find_groups

NO_AMENDMENTS = pd.DataFrame(columns=["study", "part"])
NO_LINKS = pd.DataFrame(columns=["series_a", "series_b"])


def make_strata(groups=40):
    return pd.Series([f"stratum {group % 4}" for group in range(1, groups + 1)], index=range(1, groups + 1))


def test_split_is_the_same_on_a_second_run_and_changes_with_the_seed():
    first, _ = assign_parts(make_strata(), TEST, SEED)
    again, _ = assign_parts(make_strata(), TEST, SEED)
    other, _ = assign_parts(make_strata(), TEST, SEED + 1)
    pd.testing.assert_series_equal(first, again)
    assert not first.equals(other)


def test_turn_carries_over_between_strata():
    parts, next_part = assign_parts(make_strata(41), TEST, SEED)
    counts = parts.value_counts()
    assert counts[TEST] == 21 and counts[DEVELOPMENT] == 20
    assert next_part == DEVELOPMENT


def test_later_study_joins_the_part_of_its_laboratory_group():
    earlier = pd.DataFrame({"study": ["GSE1", "GSE2"], "datasets": ["GSE1", "GSE2"], "part": [TEST, DEVELOPMENT]})
    now = pd.DataFrame({"study": ["GSE1", "GSE2", "GSE3", "GSE4"], "datasets": ["GSE1", "GSE2", "GSE3", "GSE4"],
                        "laboratory_group": [1, 2, 1, 3], "in_set_1": ["yes"] * 4,
                        "stratum": ["a|b|quiescent no"] * 4})
    parts, next_part = split_studies(now, earlier, TEST, SEED, NO_AMENDMENTS)
    assert list(parts) == [TEST, DEVELOPMENT, TEST, TEST]
    assert next_part == DEVELOPMENT


# GSE1 (`test`) and GSE2 (`development`) of the earlier split are now one laboratory group, as are GSE3 (`test`) and
# GSE4 (`development`); GSE5 is a new group of its own.
EARLIER = pd.DataFrame({"study": ["GSE1", "GSE2", "GSE3", "GSE4"], "datasets": ["GSE1", "GSE2", "GSE3", "GSE4"],
                        "part": [TEST, DEVELOPMENT, TEST, DEVELOPMENT]})
NOW = pd.DataFrame({"study": ["GSE1", "GSE2", "GSE3", "GSE4", "GSE5"], "laboratory_group": [1, 1, 2, 2, 3],
                    "in_set_1": ["yes"] * 5, "stratum": ["a|b|quiescent no"] * 5}).assign(datasets=lambda t: t.study)


def amend(rows):
    return pd.DataFrame(rows, columns=["study", "part"])


def test_amendments_move_a_study_into_development_or_out_of_both_parts():
    amendments = amend([("GSE1", DEVELOPMENT), ("GSE4", LEFT_OUT_OF_BOTH_PARTS)])
    parts, next_part = split_studies(NOW, EARLIER, TEST, SEED, amendments)
    assert list(parts) == [DEVELOPMENT, DEVELOPMENT, TEST, LEFT_OUT_OF_BOTH_PARTS, TEST]
    assert next_part == DEVELOPMENT


# GSE6 is a new study in the new laboratory group of GSE5, which takes the part `test` by its turn.
WITH_GSE6 = pd.concat([NOW, NOW.iloc[[-1]].assign(study="GSE6", datasets="GSE6")], ignore_index=True)


@pytest.mark.parametrize(("now", "amendments", "message"), [
    (NOW, amend([("GSE2", TEST), ("GSE4", LEFT_OUT_OF_BOTH_PARTS)]), "may not move a study into test"),
    (NOW, amend([("GSE1", DEVELOPMENT)]), "more than one part"),
    (WITH_GSE6, amend([("GSE1", DEVELOPMENT), ("GSE4", LEFT_OUT_OF_BOTH_PARTS), ("GSE6", DEVELOPMENT)]),
     "studies in both parts"),
])
def test_amendments_into_test_and_laboratory_groups_in_both_parts_are_refused(now, amendments, message):
    with pytest.raises(ValueError, match=message):
        split_studies(now, EARLIER, TEST, SEED, amendments)


def write_contacts(path, rows):
    columns = ["series", "contact_institute", "contact_name_sha256", "last_contributor_sha256", "contributors_sha256",
               "pubmed_ids", "relations"]
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)
    return path


@pytest.fixture
def contacts(tmp_path):
    """GSE10 is a SuperSeries without contributors in the table; its SubSeries GSE11 is not in the table and shares its
    last contributor with GSE20. GSE30 and E-MTAB-1 share nothing with anyone."""
    return write_contacts(tmp_path / "contacts.tsv", [
        ("GSE10", "Institute A", "a", "", "", "", "SuperSeries of: GSE11"),
        ("GSE11", "Institute A", "a", "x", "", "", "SubSeries of: GSE10"),
        ("GSE20", "Institute B", "b", "x", "", "", ""),
        ("GSE30", "Institute C", "c", "y", "", "", ""),
        ("E-MTAB-1", "Institute D", "d", "z", "", "", ""),
    ])


def group(contacts, links):
    study_of_dataset = pd.Series({series: series for series in ("GSE10", "GSE20", "GSE30", "E-MTAB-1")})
    lists = pd.DataFrame(columns=["method", "role", "accession"])
    outside_geo = pd.DataFrame(columns=["accession", "pubmed_id"])
    table, _ = build_laboratory_groups(study_of_dataset, contacts, lists, outside_geo, links)
    return table.set_index("series").laboratory_group.to_dict()


def test_a_series_takes_the_keys_of_its_subseries_and_links_made_by_hand_join_series(contacts):
    assert group(contacts, NO_LINKS) == {"GSE10": 1, "GSE20": 1, "GSE30": 2, "E-MTAB-1": 3}
    linked = pd.DataFrame({"series_a": ["E-MTAB-1"], "series_b": ["GSE30"]})
    assert group(contacts, linked) == {"GSE10": 1, "GSE20": 1, "GSE30": 2, "E-MTAB-1": 2}


def test_series_that_share_two_contributors_are_one_group_and_one_shared_contributor_is_not_enough(tmp_path):
    # GSE10 and GSE20 share the contributors p and q; GSE30 shares only q with them.
    contacts = write_contacts(tmp_path / "contacts.tsv", [
        ("GSE10", "Institute A", "a", "x", "p; q; r", "", ""),
        ("GSE20", "Institute B", "b", "y", "q; s; p", "", ""),
        ("GSE30", "Institute C", "c", "z", "q; t", "", ""),
        ("E-MTAB-1", "Institute D", "d", "w", "", "", ""),
    ])
    assert group(contacts, NO_LINKS) == {"GSE10": 1, "GSE20": 1, "GSE30": 2, "E-MTAB-1": 3}


def test_a_link_made_by_hand_to_a_series_outside_the_grouping_is_refused(contacts):
    with pytest.raises(ValueError, match="outside the grouping"):
        group(contacts, pd.DataFrame({"series_a": ["E-MTAB-1"], "series_b": ["GSE99"]}))


def test_groups_are_numbered_by_geo_series_first_then_other_accessions():
    links = pd.DataFrame({"series": ["GSE20", "E-MTAB-5"], "kind": ["contact", "contact"], "key": ["x", "x"]})
    groups = find_groups(["GSE300", "E-MTAB-1", "GSE20", "E-MTAB-5", "GSE9"], links)
    assert groups.to_dict() == {"GSE300": 3, "E-MTAB-1": 4, "GSE20": 2, "E-MTAB-5": 2, "GSE9": 1}
