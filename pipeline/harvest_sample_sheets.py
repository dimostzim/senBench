#!/usr/bin/env python3
"""Copy the sample sheets of downloaded datasets outside GEO into one long table, unchanged.

Read are the tab-separated sample sheets of ArrayExpress and PRIDE (SDRF), the metadata tables of SenNet, and the
dataset records of SenNet (metadata.json), from which the fields that describe the dataset, its donor and its organ
are taken.
"""

import argparse
import csv
import gzip
import json
import re
import sys
from pathlib import Path

UNITS_NAME = "units.tsv"
COLUMNS = ("dataset", "data_folder", "sheet", "row", "column", "value")
# Data folder -> name patterns of its tab-separated sample sheets, searched in every dataset folder.
SHEET_PATTERNS = {
    "biostudies_arrayexpress": ("*.sdrf.txt",),
    "proteomics_metabolomics": ("*sdrf*.tsv", "*sdrf*.txt", "*SDRF*.tsv", "*SDRF*.txt"),
    "sennet": ("*metadata.tsv",),
}
RECORD_FOLDER = "sennet"
RECORD_NAME = "metadata.json"
# Branches of a SenNet dataset record that describe the dataset, its donors and its samples.
RECORD_BRANCHES = ("metadata", "sources", "origin_samples", "direct_ancestors")
# Fields that name people, internal identifiers or processing details, and branches that repeat the others.
RECORD_SKIP = re.compile(r"user|email|contact|contributor|timestamp|uuid|history|pipeline|ingest|provider|upload|cedar|"
                         r"^files$|creation_action|local_directory|run_id|schema|_path$|^source$|^origin_samples$")
DONOR_TERM = "grouping_concept_preferred_term"
PROBE_BYTES = 4096

csv.field_size_limit(sys.maxsize)


def _one_line(text):
    return " ".join(str(text).split())


class _NotATable(Exception):
    """A file that carries the name of a sample sheet and is not a text table."""


def _sheet_rows(path):
    """Yield (row number, column, value) for every filled cell of a tab-separated sheet."""
    with path.open("rb") as handle:
        if b"\x00" in handle.read(PROBE_BYTES):
            raise _NotATable
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = [_one_line(name) for name in next(reader, [])]
        for number, row in enumerate(reader, start=1):
            for column, value in zip(header, row):
                if value.strip():
                    yield number, column, _one_line(value)


def _flatten(value, name):
    """Yield (column, value) for every filled field below a branch of a SenNet record."""
    if isinstance(value, dict):
        for key, item in value.items():
            if not RECORD_SKIP.search(key):
                yield from _flatten(item, f"{name}.{key}")
    elif isinstance(value, list) and value and all(isinstance(item, dict) for item in value):
        if all(DONOR_TERM in item for item in value):
            # Donor data are a list of terms such as Age or Sex, each with a value and a unit.
            for item in value:
                yield f"{name}.{item[DONOR_TERM]}", _one_line(f"{item.get('data_value', '')} {item.get('units', '')}")
        else:
            for index, item in enumerate(value):
                yield from _flatten(item, name if len(value) == 1 else f"{name}.{index + 1}")
    elif isinstance(value, list):
        if value:
            yield name, "; ".join(_one_line(item) for item in value)
    elif value not in (None, ""):
        yield name, _one_line(value)


def _record_rows(path):
    """Yield (row number, column, value) for the fields of a SenNet dataset record."""
    with path.open(encoding="utf-8") as handle:
        record = json.load(handle)
    for key, value in record.items():
        if key in RECORD_BRANCHES or not isinstance(value, (dict, list)):
            if not RECORD_SKIP.search(key):
                for column, text in _flatten(value, key):
                    yield 1, column, text


def _sheets(data_root):
    """Yield (data folder, dataset folder, sheet path, reader) for every sample sheet under the data root."""
    for folder, patterns in SHEET_PATTERNS.items():
        for dataset in sorted(path for path in (data_root / folder).glob("*") if path.is_dir()):
            for sheet in sorted({path for pattern in patterns for path in dataset.rglob(pattern)}):
                yield folder, dataset, sheet, _sheet_rows
            if folder == RECORD_FOLDER:
                for sheet in sorted(dataset.rglob(RECORD_NAME)):
                    yield folder, dataset, sheet, _record_rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--units-dir", type=Path, required=True)
    args = parser.parse_args()

    failed = 0
    sheets = 0
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.out, "wt", newline="") as out_file, \
            (args.units_dir / UNITS_NAME).open("w", newline="") as units_file:
        writer = csv.writer(out_file, delimiter="\t", lineterminator="\n")
        units = csv.writer(units_file, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        for folder, dataset, sheet, read in _sheets(args.data_root):
            name = str(sheet.relative_to(dataset))
            sheets += 1
            try:
                rows = list(read(sheet))
            except _NotATable:
                units.writerow((f"{dataset.name}/{name}", "skipped", "not a text table"))
                continue
            except (OSError, json.JSONDecodeError, csv.Error) as problem:
                units.writerow((f"{dataset.name}/{name}", "failed", f"cannot be read: {problem!r}"))
                failed += 1
                continue
            writer.writerows((dataset.name, folder, name, number, column, value) for number, column, value in rows)
            units.writerow((f"{dataset.name}/{name}", "ok", f"{len({number for number, _, _ in rows})} rows"))
            if sheets % 500 == 0:
                print(f"PROGRESS {sheets} sheets read, the last: {dataset.name}/{name}", flush=True)
    print(f"{sheets} sheets read, {failed} failed", flush=True)
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
