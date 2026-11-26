#!/usr/bin/env python3
"""Copy the sample and series descriptions of downloaded GEO series from their SOFT files into tables, unchanged."""

import argparse
import csv
import gzip
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

UNITS_NAME = "units.tsv"
SAMPLES_NAME = "geo_samples.tsv.gz"
CHARACTERISTICS_NAME = "geo_sample_characteristics.tsv.gz"
SERIES_NAME = "geo_series.tsv"
SEPARATOR = "; "
CHANNELS = ("ch1", "ch2")

# Table column -> SOFT field. A field that occurs several times in a record is joined with SEPARATOR.
SAMPLE_FIELDS = {
    "title": "Sample_title",
    "type": "Sample_type",
    "channel_count": "Sample_channel_count",
    "source_name_ch1": "Sample_source_name_ch1",
    "organism_ch1": "Sample_organism_ch1",
    "taxid_ch1": "Sample_taxid_ch1",
    "characteristics_ch1": "Sample_characteristics_ch1",
    "molecule_ch1": "Sample_molecule_ch1",
    "treatment_protocol_ch1": "Sample_treatment_protocol_ch1",
    "growth_protocol_ch1": "Sample_growth_protocol_ch1",
    "source_name_ch2": "Sample_source_name_ch2",
    "organism_ch2": "Sample_organism_ch2",
    "characteristics_ch2": "Sample_characteristics_ch2",
    "description": "Sample_description",
    "platform": "Sample_platform_id",
    "instrument_model": "Sample_instrument_model",
    "library_strategy": "Sample_library_strategy",
    "library_source": "Sample_library_source",
    "library_selection": "Sample_library_selection",
    "relations": "Sample_relation",
    "supplementary_files": "Sample_supplementary_file",
    "series_of_sample": "Sample_series_id",
    "submission_date": "Sample_submission_date",
}
SERIES_FIELDS = {
    "title": "Series_title",
    "types": "Series_type",
    "summary": "Series_summary",
    "overall_design": "Series_overall_design",
    "pubmed_ids": "Series_pubmed_id",
    "relations": "Series_relation",
    "platforms": "Series_platform_id",
    "submission_date": "Series_submission_date",
    "last_update_date": "Series_last_update_date",
}
SAMPLE_COLUMNS = ("series", "sample", *SAMPLE_FIELDS)
SERIES_COLUMNS = ("series", *SERIES_FIELDS, "samples")
CHARACTERISTICS_COLUMNS = ("series", "sample", "channel", "key", "value")

csv.field_size_limit(sys.maxsize)


def _one_line(text):
    return " ".join(text.split())


def _records(path):
    """Yield (kind, accession, fields) for every record of a SOFT file; fields maps a field name to its values."""
    kind = accession = None
    fields = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line[0] == "^":
                if kind:
                    yield kind, accession, fields
                kind, _, accession = line[1:].partition("=")
                kind, accession, fields = kind.strip(), accession.strip(), {}
            # Lines without "!" are rows of the data tables, which are not read.
            elif line[0] == "!" and " = " in line:
                name, _, value = line[1:].partition(" = ")
                # A supplementary file field carries a number, as in Sample_supplementary_file_1.
                if name.startswith("Sample_supplementary_file"):
                    name = "Sample_supplementary_file"
                fields.setdefault(name, []).append(_one_line(value))
    if kind:
        yield kind, accession, fields


def _row(accession_columns, columns, fields):
    return {**accession_columns, **{column: SEPARATOR.join(fields.get(name, [])) for column, name in columns.items()}}


def _harvest(series_dir):
    """Return (series, status, note, series row, sample rows, characteristics rows) for one series folder."""
    series = series_dir.name
    soft = series_dir / "soft" / f"{series}_family.soft.gz"
    if not soft.is_file():
        return series, "failed", "no SOFT file", None, [], []
    series_row = None
    samples = []
    characteristics = []
    try:
        for kind, accession, fields in _records(soft):
            if kind == "SERIES":
                series_row = _row({"series": series}, SERIES_FIELDS, fields)
            elif kind == "SAMPLE":
                samples.append(_row({"series": series, "sample": accession}, SAMPLE_FIELDS, fields))
                for channel in CHANNELS:
                    for entry in fields.get(f"Sample_characteristics_{channel}", []):
                        key, separator, value = entry.partition(": ")
                        characteristics.append((series, accession, channel, key if separator else "",
                                                value if separator else entry))
    except (OSError, EOFError) as problem:
        return series, "failed", f"SOFT file cannot be read: {problem}", None, [], []
    if series_row is None:
        return series, "failed", "no series record in the SOFT file", None, [], []
    series_row["samples"] = len(samples)
    return series, "ok", f"{len(samples)} samples", series_row, samples, characteristics


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--geo-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--units-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()

    series_dirs = sorted(path for path in args.geo_root.iterdir() if path.is_dir())[:args.limit]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    failed = []
    with gzip.open(args.out_dir / SAMPLES_NAME, "wt", newline="") as samples_file, \
            gzip.open(args.out_dir / CHARACTERISTICS_NAME, "wt", newline="") as characteristics_file, \
            (args.out_dir / SERIES_NAME).open("w", newline="") as series_file, \
            (args.units_dir / UNITS_NAME).open("w", newline="") as units_file, \
            ProcessPoolExecutor(max_workers=args.workers) as pool:
        samples_writer = csv.DictWriter(samples_file, SAMPLE_COLUMNS, delimiter="\t", lineterminator="\n")
        series_writer = csv.DictWriter(series_file, SERIES_COLUMNS, delimiter="\t", lineterminator="\n")
        characteristics_writer = csv.writer(characteristics_file, delimiter="\t", lineterminator="\n")
        units_writer = csv.writer(units_file, delimiter="\t", lineterminator="\n")
        samples_writer.writeheader()
        series_writer.writeheader()
        characteristics_writer.writerow(CHARACTERISTICS_COLUMNS)
        results = pool.map(_harvest, series_dirs, chunksize=8)
        for done, (series, status, note, series_row, samples, characteristics) in enumerate(results, start=1):
            units_writer.writerow((series, status, note))
            if status == "ok":
                series_writer.writerow(series_row)
                samples_writer.writerows(samples)
                characteristics_writer.writerows(characteristics)
            else:
                failed.append(series)
            if done % 500 == 0 or done == len(series_dirs):
                print(f"PROGRESS {done}/{len(series_dirs)} {series} {status}: {note}", flush=True)
    if failed:
        print(f"failed series ({len(failed)}): {' '.join(failed)}", flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
