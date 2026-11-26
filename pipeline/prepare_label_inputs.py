#!/usr/bin/env python3
"""Collect, per dataset, the series and sample fields that a label reading reads, unchanged: from GEO, or for an
ArrayExpress experiment from its sample sheet (SDRF) and its investigation description (IDF)."""

import argparse
import csv
import gzip
import itertools
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

from list_bulk_human_datasets import BULK_RNA_SEQ_ASSAYS

SERIES_FIELDS = ("title", "summary", "overall_design")
# Name in the output, and the column of the GEO sample table it is copied from.
SAMPLE_FIELDS = (
    ("title", "title"),
    ("source_name", "source_name_ch1"),
    ("organism", "organism_ch1"),
    ("molecule", "molecule_ch1"),
    ("treatment_protocol", "treatment_protocol_ch1"),
    ("growth_protocol", "growth_protocol_ch1"),
    ("description", "description"),
    ("library_source", "library_source"),
    ("library_selection", "library_selection"),
)
CHANNEL = "ch1"
ARRAYEXPRESS = "ArrayExpress"
# The data folder of the ArrayExpress experiments in the table of sample sheets.
ARRAYEXPRESS_FOLDER = "biostudies_arrayexpress"
IDF_PATTERN = "*.idf.txt"
# Series field, and the IDF row whose values it holds.
IDF_SERIES_ROWS = (("title", "Investigation Title"), ("summary", "Experiment Description"),
                   ("overall_design", "Experimental Design"))
# Sample field, and the sample sheet column it is copied from.
SHEET_FIELDS = (("organism", "Characteristics[organism]"), ("description", "Description"),
                ("library_source", "Comment[LIBRARY_SOURCE]"), ("library_selection", "Comment[LIBRARY_SELECTION]"))
# Sample field, and the IDF protocol type whose descriptions it holds. GEO has no field for the sample collection
# protocol, in which ArrayExpress submitters often describe the culture and the treatment.
PROTOCOL_FIELDS = (("treatment_protocol", "treatment protocol"), ("growth_protocol", "growth protocol"),
                   ("sample_collection_protocol", "sample collection protocol"))
ARRAYEXPRESS_SAMPLE_FIELDS = ("sample", "characteristics", "title", "source_name", "organism", "molecule",
                              "treatment_protocol", "growth_protocol", "sample_collection_protocol", "description",
                              "library_source", "library_selection")
# Sample sheet columns that describe the sample; they are given as characteristics, with their column names.
_DESCRIBING_COLUMN = re.compile(r"(?:characteristics|factor ?value|unit) ?\[", re.I)
SOURCE_NAME = "sourcename"
EXTRACT_NAME = "extractname"
MATERIAL_TYPE = "materialtype"
PROTOCOL_REF = "protocolref"
VALUE_SEPARATOR = "; "

csv.field_size_limit(sys.maxsize)


def _read_rows(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def _column_key(column):
    """Return a sample sheet column name in one form: lower case, without spaces."""
    return column.replace(" ", "").lower()


def _join(values):
    return VALUE_SEPARATOR.join(dict.fromkeys(value for value in values if value))


def read_idf(path):
    """Return the rows of an investigation description (IDF): its tag mapped to its values, empty values kept."""
    with path.open(newline="", encoding="utf-8") as handle:
        return {row[0].strip(): [value.strip() for value in row[1:]]
                for row in csv.reader(handle, delimiter="\t") if row and row[0].strip()}


def arrayexpress_series(idf):
    """Return the series fields of an ArrayExpress experiment from its IDF."""
    return {field: _join(idf.get(tag, [])) for field, tag in IDF_SERIES_ROWS}


def _protocols(idf):
    """Return, per protocol accession of the IDF, its type and its description."""
    rows = (idf.get(tag, []) for tag in ("Protocol Name", "Protocol Type", "Protocol Description"))
    return {name: (kind, description) for name, kind, description in itertools.zip_longest(*rows, fillvalue="")
            if name}


def _extract_material_types(row):
    """Return the material types that a sample sheet row gives after its extract name, that is of the extract."""
    keys = [_column_key(column) for column, _ in row]
    if EXTRACT_NAME not in keys:
        return []
    start = keys.index(EXTRACT_NAME)
    return [value for key, (_, value) in zip(keys[start:], row[start:]) if key == MATERIAL_TYPE]


def arrayexpress_sample(sample, rows, idf):
    """Return the label fields of one sample of an ArrayExpress experiment.

    `rows` are the sample sheet rows of the sample (one per data file), each a list of (column, value) in the order of
    the sheet; `idf` is the experiment's IDF as read_idf() returns it.
    """
    cells = [cell for row in rows for cell in row]
    fields = {"sample": sample, "title": sample, "source_name": sample,
              "characteristics": [list(cell) for cell in dict.fromkeys(cells) if _DESCRIBING_COLUMN.match(cell[0])],
              "molecule": _join(value for row in rows for value in _extract_material_types(row))}
    for field, column in SHEET_FIELDS:
        fields[field] = _join(value for name, value in cells if _column_key(name) == _column_key(column))
    protocols = _protocols(idf)
    references = [value for name, value in cells if _column_key(name) == PROTOCOL_REF]
    for field, kind in PROTOCOL_FIELDS:
        fields[field] = _join(protocols[reference][1] for reference in references if protocols[reference][0] == kind)
    return {field: fields[field] for field in ARRAYEXPRESS_SAMPLE_FIELDS}


def read_sheet_rows(path, datasets):
    """Return, per (ArrayExpress dataset, source name), the sample sheet rows of the sample from the table of sample
    sheets; each row is a list of (column, value) in the order of the sheet."""
    rows = defaultdict(list)
    for cell in _read_rows(path):
        if cell["dataset"] in datasets and cell["data_folder"] == ARRAYEXPRESS_FOLDER:
            rows[(cell["dataset"], cell["sheet"], cell["row"])].append((cell["column"], cell["value"]))
    samples = defaultdict(list)
    for (dataset, sheet, number), row in rows.items():
        names = [value for column, value in row if _column_key(column) == SOURCE_NAME]
        if len(names) != 1:
            raise SystemExit(f"{dataset}: row {number} of {sheet} has {len(names)} source names")
        samples[(dataset, names[0])].append(row)
    return samples


def find_idf(folder):
    idfs = sorted(folder.glob(IDF_PATTERN))
    if len(idfs) != 1:
        raise SystemExit(f"{folder}: {len(idfs)} IDF files")
    return idfs[0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--datasets", type=Path, required=True, help="table with a column `dataset`")
    parser.add_argument("--samples", type=Path, required=True, help="the sample table")
    parser.add_argument("--geo-series", type=Path, required=True)
    parser.add_argument("--geo-samples", type=Path, required=True)
    parser.add_argument("--geo-characteristics", type=Path, required=True)
    parser.add_argument("--sample-sheets", type=Path, required=True,
                        help="the table of sample sheets outside GEO, for the ArrayExpress experiments")
    parser.add_argument("--arrayexpress-dir", type=Path, required=True,
                        help="folder with one folder per ArrayExpress experiment, holding its IDF")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    datasets = [row["dataset"] for row in _read_rows(args.datasets)]
    wanted = set(datasets)
    samples = {dataset: {} for dataset in datasets}
    arrayexpress = set()
    for row in _read_rows(args.samples):
        if row["dataset"] in wanted and row["level"] == "sample" and row["assay"] in BULK_RNA_SEQ_ASSAYS:
            if row["sample"] in samples[row["dataset"]]:
                raise SystemExit(f"{row['dataset']}: sample {row['sample']} has two rows in the sample table")
            samples[row["dataset"]][row["sample"]] = {"sample": row["sample"], "characteristics": []}
            if row["repository"] == ARRAYEXPRESS:
                arrayexpress.add(row["dataset"])
    print(f"PROGRESS 1/5 sample table read: {sum(map(len, samples.values()))} samples of {len(datasets)} datasets",
          flush=True)

    series = {row["series"]: {field: row[field] for field in SERIES_FIELDS}
              for row in _read_rows(args.geo_series) if row["series"] in wanted}
    print(f"PROGRESS 2/5 series fields read: {len(series)} series", flush=True)

    seen = set()
    for row in _read_rows(args.geo_samples):
        sample = samples.get(row["series"], {}).get(row["sample"])
        if sample is not None:
            sample.update({name: row[column] for name, column in SAMPLE_FIELDS})
            seen.add((row["series"], row["sample"]))
    print(f"PROGRESS 3/5 sample fields read: {len(seen)} samples", flush=True)

    for row in _read_rows(args.geo_characteristics):
        sample = samples.get(row["series"], {}).get(row["sample"])
        if sample is not None and row["channel"] == CHANNEL:
            sample["characteristics"].append([row["key"], row["value"]])
    print("PROGRESS 4/5 characteristics read", flush=True)

    sheet_rows = read_sheet_rows(args.sample_sheets, arrayexpress)
    for dataset in sorted(arrayexpress):
        idf = read_idf(find_idf(args.arrayexpress_dir / dataset))
        series[dataset] = arrayexpress_series(idf)
        for sample in samples[dataset]:
            if (dataset, sample) in sheet_rows:
                samples[dataset][sample] = arrayexpress_sample(sample, sheet_rows[(dataset, sample)], idf)
                seen.add((dataset, sample))
    print(f"PROGRESS 5/5 sample sheets and IDF read: {len(arrayexpress)} ArrayExpress experiments", flush=True)

    missing_series = wanted - set(series)
    missing_samples = {(dataset, sample) for dataset in datasets for sample in samples[dataset]} - seen
    if missing_series or missing_samples:
        raise SystemExit(f"without fields: {len(missing_series)} series, {len(missing_samples)} samples; "
                         f"the first: {sorted(missing_series)[:3]} {sorted(missing_samples)[:3]}")

    with gzip.open(args.out, "wt") as handle:
        for dataset in datasets:
            record = {"dataset": dataset, "series": series[dataset], "samples": list(samples[dataset].values())}
            handle.write(json.dumps(record, ensure_ascii=False) + "\n")
    print(f"{len(datasets)} datasets written to {args.out}")


if __name__ == "__main__":
    main()
