"""Check the fields that link series of one laboratory, read from an ArrayExpress investigation description."""

import pytest

from harvest_geo_contacts import describe_experiment, hash_name, normalise_name
from mark_tool_training_overlap import RELATION


def test_experiment_takes_the_submitter_as_contact_and_the_last_person_as_last_contributor():
    idf = {"Person First Name": ["Ann", "Bo"], "Person Mid Initials": ["", "K"], "Person Last Name": ["Lee", "Wu"],
           "Person Affiliation": ["Institute A", "Institute B"], "Person Roles": ["biomaterial provider;submitter",
                                                                                  "investigator"],
           "PubMed ID": [""], "Comment[SecondaryAccession]": ["ERP172172"]}
    row = describe_experiment("E-MTAB-1", idf)
    assert (row["contact_name"], row["contact_institute"], row["last_contributor"]) == ("Ann,,Lee", "Institute A",
                                                                                         "Bo,K,Wu")
    assert row["pubmed_ids"] == ""
    assert RELATION.findall(row["relations"]) == [("SRA", "ERP172172")]


@pytest.mark.parametrize("tag", ["Comment[SecondaryAccession]", "Comment [SecondaryAccession]"])
def test_experiment_reads_the_secondary_accession_with_and_without_a_space(tag):
    idf = {"Person First Name": ["Ann"], "Person Last Name": ["Lee"], "Person Affiliation": ["Institute A"],
           "Person Roles": ["submitter"], tag: ["ERP117860"]}
    assert RELATION.findall(describe_experiment("E-MTAB-1", idf)["relations"]) == [("SRA", "ERP117860")]


@pytest.mark.parametrize("one, other", [
    ("University Medical Center Göttingen", "University Medical Center Goettingen"),
    ("Université de Montréal", "Universite de Montreal"),
    ("Karolinska Institutet, Stockholm", "karolinska institutet stockholm"),
])
def test_names_written_with_and_without_umlauts_accents_or_punctuation_are_equal(one, other):
    assert normalise_name(one) == normalise_name(other)


@pytest.mark.parametrize("one, other", [
    ("Ann,,Lee", "Lee,,Ann"),
    ("Ann,K,Lee", "Ann,,Lee"),
    ("Ann,,Lee", "Ann Lee"),
    ("Jürgen,,Groß", "Juergen,,Gross"),
])
def test_a_person_has_one_hash_with_first_and_last_name_exchanged_or_without_the_middle_name(one, other):
    assert hash_name(one) == hash_name(other)


def test_different_people_and_an_empty_name_have_different_hashes():
    assert hash_name("Ann,,Lee") != hash_name("Ann,,Li")
    assert hash_name("") == hash_name(",,") == ""
