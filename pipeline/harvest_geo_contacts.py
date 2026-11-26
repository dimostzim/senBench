#!/usr/bin/env python3
"""Read the fields that link GEO series of one laboratory from their SOFT files: contact, contributors, PubMed ids, relations.
ArrayExpress experiments are read from their investigation description (IDF) into the same fields.

Two tables are written. The one with names holds personal data of the submitters and stays where it is written.
The other holds a sha256 of every normalised name in place of the name and is the input of the laboratory groups.
"""

import argparse
import hashlib
import os
import re
import string
import unicodedata
from itertools import zip_longest
from pathlib import Path

import pandas as pd

from harvest_geo_samples import SEPARATOR, _records
from prepare_label_inputs import find_idf, read_idf

WITH_NAMES_NAME = "geo_series_contacts_with_names.tsv"
HASHED_NAME = "geo_series_contacts.tsv"
NOT_READ_NAME = "geo_series_contacts_not_read.tsv"
GEO_SERIES = re.compile(r"GSE\d+")
ARRAYEXPRESS_EXPERIMENT = re.compile(r"E-[A-Z]+-\d+")
SUBMITTER = "submitter"
# Secondary accession of an IDF, and the relation that a GEO series record gives an accession of that kind.
SECONDARY_RELATIONS = ((re.compile(r"[EDS]RP\d+"), "SRA"), (re.compile(r"PRJ[A-Z]+\d+"), "BioProject"))
# IDF files write the tag of the secondary accessions with and without a space.
SECONDARY_ACCESSION_TAGS = ("Comment[SecondaryAccession]", "Comment [SecondaryAccession]")
_PUNCTUATION = str.maketrans(string.punctuation, " " * len(string.punctuation))
# Records write a German umlaut either as the letter or as two letters (Göttingen, Goettingen).
_UMLAUTS = str.maketrans({"ä": "ae", "ö": "oe", "ü": "ue", "ß": "ss"})


def normalise_name(name):
    """Return a name in lower case without punctuation and accents, with German umlauts written as two letters."""
    letters = unicodedata.normalize("NFKD", unicodedata.normalize("NFC", name.lower()).translate(_UMLAUTS))
    unaccented = "".join(letter for letter in letters if not unicodedata.combining(letter))
    return " ".join(unaccented.translate(_PUNCTUATION).split())


def hash_name(name):
    """Return the sha256 of the words of a person's normalised first and last name in sorted order, so that a record
    with the two exchanged or without the middle name gives the same hash; an empty text for an empty name. GEO writes
    a person as `First,Middle,Last`."""
    fields = name.split(",")
    words = sorted(word for field in fields[:1] + fields[1:][-1:] for word in normalise_name(field).split())
    return hashlib.sha256(" ".join(words).encode()).hexdigest() if words else ""


def read_series_fields(soft):
    """Return the fields of the series record of a SOFT file, which stands before the platform and sample records."""
    for kind, _, fields in _records(soft):
        if kind == "SERIES":
            return fields
    raise ValueError(f"no series record in {soft}")


def describe_series(series, fields):
    """Return the row of one series for the table with names."""
    contributors = fields.get("Series_contributor", [])
    return {
        "series": series,
        "contact_name": SEPARATOR.join(fields.get("Series_contact_name", [])),
        "contact_institute": SEPARATOR.join(fields.get("Series_contact_institute", [])),
        "contributors": SEPARATOR.join(contributors),
        "last_contributor": contributors[-1] if contributors else "",
        "pubmed_ids": SEPARATOR.join(fields.get("Series_pubmed_id", [])),
        "relations": SEPARATOR.join(fields.get("Series_relation", [])),
    }


def _relate_secondary_accession(accession):
    for pattern, relation in SECONDARY_RELATIONS:
        if pattern.fullmatch(accession):
            return f"{relation}: {accession}"
    raise ValueError(f"secondary accession of an unknown kind: {accession}")


def describe_experiment(experiment, idf):
    """Return the row of one ArrayExpress experiment for the table with names, in the fields of a GEO series: the person
    with the role `submitter` is the contact, and the last person listed is the last contributor."""
    people = list(zip_longest(*(idf.get(tag, []) for tag in (
        "Person First Name", "Person Mid Initials", "Person Last Name", "Person Affiliation", "Person Roles")),
        fillvalue=""))
    names = [",".join(person[:3]) for person in people]
    submitters = [person for person in people if SUBMITTER in person[4].split(";")]
    if len(submitters) != 1:
        raise ValueError(f"{experiment}: {len(submitters)} persons with the role {SUBMITTER}")
    return {
        "series": experiment,
        "contact_name": ",".join(submitters[0][:3]),
        "contact_institute": submitters[0][3],
        "contributors": SEPARATOR.join(names),
        "last_contributor": names[-1],
        "pubmed_ids": SEPARATOR.join(value for value in idf.get("PubMed ID", []) if value),
        "relations": SEPARATOR.join(_relate_secondary_accession(value) for tag in SECONDARY_ACCESSION_TAGS
                                    for value in idf.get(tag, []) if value),
    }


def hash_names(with_names):
    """Return the table without names: a sha256 of the contact name, of the last contributor name and of every
    contributor name instead."""
    return pd.DataFrame({
        "series": with_names.series,
        "contact_institute": with_names.contact_institute,
        "contact_name_sha256": with_names.contact_name.map(hash_name),
        "last_contributor_sha256": with_names.last_contributor.map(hash_name),
        "contributors_sha256": with_names.contributors.map(
            lambda names: SEPARATOR.join(hash_name(name) for name in names.split(SEPARATOR) if name)),
        "pubmed_ids": with_names.pubmed_ids,
        "relations": with_names.relations,
    })


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--geo-roots", type=Path, nargs="+", required=True,
                        help="folders with <series>/soft/<series>_family.soft.gz; the first that has the file is read")
    parser.add_argument("--accessions", nargs=2, action="append", required=True,
                        metavar=("TABLE", "COLUMN"), help="a table and its column of accessions; can be repeated")
    parser.add_argument("--arrayexpress-dir", type=Path,
                        help="folder with one folder per ArrayExpress experiment, holding its IDF; needed when the "
                             "accessions include ArrayExpress experiments")
    parser.add_argument("--with-names-dir", type=Path, required=True,
                        help="folder for the table with names; personal data, not for a repository")
    parser.add_argument("--hashed-dir", type=Path, required=True, help="folder for the table with hashes of names")
    args = parser.parse_args()
    for path in (args.with_names_dir / WITH_NAMES_NAME, args.hashed_dir / HASHED_NAME):
        if path.exists():
            raise SystemExit(f"{path} exists")

    accessions = pd.concat(pd.read_csv(table, sep="\t", dtype=str, keep_default_na=False)[column]
                           for table, column in args.accessions).drop_duplicates()
    series = sorted(accessions[accessions.str.fullmatch(GEO_SERIES)], key=lambda accession: int(accession[3:]))
    rows, not_read = [], []
    for accession in series:
        soft_files = [root / accession / "soft" / f"{accession}_family.soft.gz" for root in args.geo_roots]
        soft_files = [soft for soft in soft_files if soft.is_file()]
        if soft_files:
            rows.append(describe_series(accession, read_series_fields(soft_files[0])))
        else:
            not_read.append({"series": accession, "reason": "no SOFT file in the GEO folders"})
    experiments = sorted(accessions[accessions.str.fullmatch(ARRAYEXPRESS_EXPERIMENT)])
    if experiments and args.arrayexpress_dir is None:
        raise SystemExit(f"{len(experiments)} ArrayExpress experiments among the accessions: give --arrayexpress-dir")
    rows += [describe_experiment(experiment, read_idf(find_idf(args.arrayexpress_dir / experiment)))
             for experiment in experiments]
    with_names = pd.DataFrame(rows)

    args.with_names_dir.mkdir(parents=True, exist_ok=True)
    args.hashed_dir.mkdir(parents=True, exist_ok=True)
    # Readable by the owner only, from the moment it exists.
    descriptor = os.open(args.with_names_dir / WITH_NAMES_NAME, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    with open(descriptor, "w", newline="") as handle:
        with_names.to_csv(handle, sep="\t", index=False)
    hash_names(with_names).to_csv(args.hashed_dir / HASHED_NAME, sep="\t", index=False)
    pd.DataFrame(not_read, columns=["series", "reason"]).to_csv(args.hashed_dir / NOT_READ_NAME, sep="\t", index=False)
    print(f"{len(accessions)} accessions, {len(series)} of them GEO series: {len(series) - len(not_read)} read, "
          f"{len(not_read)} without a SOFT file; {len(experiments)} ArrayExpress experiments read from their IDF")


if __name__ == "__main__":
    main()
