#!/usr/bin/env python3
import argparse
import csv
import gzip
import itertools
import json
import re
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

GEO_FOLDER = "geo"
ARRAYEXPRESS_FOLDER = "biostudies_arrayexpress"
SEPARATOR = "; "
GEO_SAMPLES = "geo_samples.tsv.gz"
GEO_CHARACTERISTICS = "geo_sample_characteristics.tsv.gz"
GEO_SERIES = "geo_series.tsv"
SHEETS_NAME = "sample_sheets_downloads.tsv.gz"
SAMPLES_NAME = "samples.tsv.gz"
LABEL_INPUTS_NAME = "label_inputs.jsonl.gz"
BULK_RNA_SEQ_ASSAYS = ("RNA-Seq", "RNA-Seq; sequencing assay")

csv.field_size_limit(sys.maxsize)


def _one_line(text):
    return " ".join(str(text).split())


CHANNELS = ("ch1", "ch2")
SOFT_SAMPLE_FIELDS = {
    "title": "Sample_title",
    "type": "Sample_type",
    "source_name_ch1": "Sample_source_name_ch1",
    "organism_ch1": "Sample_organism_ch1",
    "molecule_ch1": "Sample_molecule_ch1",
    "treatment_protocol_ch1": "Sample_treatment_protocol_ch1",
    "growth_protocol_ch1": "Sample_growth_protocol_ch1",
    "description": "Sample_description",
    "library_strategy": "Sample_library_strategy",
    "library_source": "Sample_library_source",
    "library_selection": "Sample_library_selection",
}
SOFT_SERIES_FIELDS = {
    "title": "Series_title",
    "summary": "Series_summary",
    "overall_design": "Series_overall_design",
}
GEO_SAMPLE_COLUMNS = ("series", "sample", *SOFT_SAMPLE_FIELDS)
GEO_SERIES_COLUMNS = ("series", *SOFT_SERIES_FIELDS, "samples")
GEO_CHARACTERISTICS_COLUMNS = ("series", "sample", "channel", "key", "value")


def _records(path):
    kind = accession = None
    fields = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line[0] == "^":
                if kind:
                    yield kind, accession, fields
                kind, _, accession = line[1:].partition("=")
                kind, accession, fields = kind.strip(), accession.strip(), {}
            elif line[0] == "!" and " = " in line:
                name, _, value = line[1:].partition(" = ")
                fields.setdefault(name, []).append(_one_line(value))
    if kind:
        yield kind, accession, fields


def _row(accession_columns, columns, fields):
    return {**accession_columns, **{column: SEPARATOR.join(fields.get(name, [])) for column, name in columns.items()}}


def _harvest(series_dir):
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
                series_row = _row({"series": series}, SOFT_SERIES_FIELDS, fields)
            elif kind == "SAMPLE":
                samples.append(_row({"series": series, "sample": accession}, SOFT_SAMPLE_FIELDS, fields))
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


def harvest_geo_records(geo_dir, out_dir, workers):
    series_dirs = sorted(path for path in geo_dir.iterdir() if path.is_dir())
    failed = []
    with gzip.open(out_dir / GEO_SAMPLES, "wt", newline="") as samples_file, \
            gzip.open(out_dir / GEO_CHARACTERISTICS, "wt", newline="") as characteristics_file, \
            (out_dir / GEO_SERIES).open("w", newline="") as series_file, \
            ProcessPoolExecutor(max_workers=workers) as pool:
        samples_writer = csv.DictWriter(samples_file, GEO_SAMPLE_COLUMNS, delimiter="\t", lineterminator="\n")
        series_writer = csv.DictWriter(series_file, GEO_SERIES_COLUMNS, delimiter="\t", lineterminator="\n")
        characteristics_writer = csv.writer(characteristics_file, delimiter="\t", lineterminator="\n")
        samples_writer.writeheader()
        series_writer.writeheader()
        characteristics_writer.writerow(GEO_CHARACTERISTICS_COLUMNS)
        results = pool.map(_harvest, series_dirs, chunksize=8)
        for done, (series, status, note, series_row, samples, characteristics) in enumerate(results, start=1):
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


SHEET_COLUMNS = ("dataset", "data_folder", "sheet", "row", "column", "value")
SHEET_PATTERN = "*.sdrf.txt"
PROBE_BYTES = 4096


class _NotATable(Exception):
    pass


def _sheet_rows(path):
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


def harvest_sample_sheets(arrayexpress_dir, out_dir):
    failed = 0
    sheets = 0
    with gzip.open(out_dir / SHEETS_NAME, "wt", newline="") as out_file:
        writer = csv.writer(out_file, delimiter="\t", lineterminator="\n")
        writer.writerow(SHEET_COLUMNS)
        for dataset in sorted(path for path in arrayexpress_dir.glob("*") if path.is_dir()):
            for sheet in sorted(dataset.rglob(SHEET_PATTERN)):
                name = str(sheet.relative_to(dataset))
                sheets += 1
                try:
                    rows = list(_sheet_rows(sheet))
                except _NotATable:
                    continue
                except (OSError, csv.Error) as problem:
                    print(f"{dataset.name}/{name} cannot be read: {problem!r}", flush=True)
                    failed += 1
                    continue
                writer.writerows((dataset.name, ARRAYEXPRESS_FOLDER, name, number, column, value)
                                 for number, column, value in rows)
                if sheets % 500 == 0:
                    print(f"PROGRESS {sheets} sheets read, the last: {dataset.name}/{name}", flush=True)
    print(f"{sheets} sheets read, {failed} failed", flush=True)
    if failed:
        sys.exit(1)


SAMPLE_COLUMNS = ("dataset", "repository", "sample", "organism", "assay")
GEO_REPOSITORY = "GEO"
ARRAYEXPRESS_REPOSITORY = "ArrayExpress"
GEO_ARRAY_ASSAYS = {"RNA": "array, RNA", "genomic": "array, genomic DNA", "protein": "array, protein"}


def _field_name(name):
    name = name.strip()
    wrapped = re.fullmatch(r"(?:characteristics|factor ?value|comment|parameter ?value)\s*\[(.+)\]", name, re.I)
    name = wrapped.group(1) if wrapped else name.rsplit(".", 1)[-1]
    name = re.sub(r"(?<=[a-z])(?=[A-Z])", " ", name).replace("_", " ")
    return " ".join(name.lower().split())


def _join_organisms(fields):
    return SEPARATOR.join(dict.fromkeys(value for name, value in fields if _field_name(name) == "organism"))


def _geo_assay(row):
    if row["type"] == "SRA":
        single_cell = "single cell" in row["library_source"]
        return f"{row['library_strategy']}{', single cell' if single_cell else ''}"
    return GEO_ARRAY_ASSAYS.get(row["type"], row["type"])


def _geo_datasets(samples_dir):
    with gzip.open(samples_dir / GEO_SAMPLES, "rt", newline="") as samples_file, \
            gzip.open(samples_dir / GEO_CHARACTERISTICS, "rt", newline="") as characteristics_file:
        characteristics = itertools.groupby(csv.DictReader(characteristics_file, delimiter="\t"),
                                            key=lambda row: (row["series"], row["sample"]))
        current = next(characteristics, None)
        for series, rows in itertools.groupby(csv.DictReader(samples_file, delimiter="\t"),
                                              key=lambda row: row["series"]):
            samples = []
            for row in rows:
                fields = [("organism", row["organism_ch1"])]
                if current and current[0] == (series, row["sample"]):
                    fields += [(item["key"], item["value"]) for item in current[1] if item["channel"] == "ch1"]
                    current = next(characteristics, None)
                samples.append({"sample": row["sample"], "assay": _geo_assay(row), "fields": fields})
            yield series, samples
        if current:
            raise SystemExit(f"characteristics of {current[0]} do not follow the order of the sample table")


def _sheet_datasets(sheets_path):
    datasets = defaultdict(dict)
    with gzip.open(sheets_path, "rt", newline="") as handle:
        for (dataset, sheet, row), cells in itertools.groupby(
                csv.DictReader(handle, delimiter="\t"), key=lambda cell: (cell["dataset"], cell["sheet"], cell["row"])):
            fields = [(cell["column"], cell["value"]) for cell in cells]
            named = defaultdict(list)
            for name, value in fields:
                named[_field_name(name)].append(value)
            sample = named["source name"][0] if named["source name"] else f"{sheet} row {row}"
            assay = SEPARATOR.join(dict.fromkeys(named["library strategy"] + named["technology type"]))
            entry = datasets[dataset].setdefault(sample, {"sample": sample, "assay": assay, "fields": []})
            entry["fields"] += [field for field in fields if field not in entry["fields"]]
    for dataset, samples in sorted(datasets.items()):
        yield dataset, list(samples.values())


def build_sample_table(out_dir):
    rows = 0
    with gzip.open(out_dir / SAMPLES_NAME, "wt", newline="") as out_file:
        writer = csv.writer(out_file, delimiter="\t", lineterminator="\n")
        writer.writerow(SAMPLE_COLUMNS)
        for repository, datasets in ((GEO_REPOSITORY, _geo_datasets(out_dir)),
                                     (ARRAYEXPRESS_REPOSITORY, _sheet_datasets(out_dir / SHEETS_NAME))):
            for dataset, samples in datasets:
                writer.writerows((dataset, repository, sample["sample"], _join_organisms(sample["fields"]),
                                  sample["assay"]) for sample in samples)
                rows += len(samples)
    print(f"{rows} samples written to {out_dir / SAMPLES_NAME}")


LABEL_SERIES_FIELDS = ("title", "summary", "overall_design")
LABEL_SAMPLE_FIELDS = (
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
IDF_PATTERN = "*.idf.txt"
IDF_SERIES_ROWS = (("title", "Investigation Title"), ("summary", "Experiment Description"),
                   ("overall_design", "Experimental Design"))
SHEET_FIELDS = (("organism", "Characteristics[organism]"), ("description", "Description"),
                ("library_source", "Comment[LIBRARY_SOURCE]"), ("library_selection", "Comment[LIBRARY_SELECTION]"))
PROTOCOL_FIELDS = (("treatment_protocol", "treatment protocol"), ("growth_protocol", "growth protocol"),
                   ("sample_collection_protocol", "sample collection protocol"))
ARRAYEXPRESS_SAMPLE_FIELDS = ("sample", "characteristics", "title", "source_name", "organism", "molecule",
                              "treatment_protocol", "growth_protocol", "sample_collection_protocol", "description",
                              "library_source", "library_selection")
_DESCRIBING_COLUMN = re.compile(r"(?:characteristics|factor ?value|unit) ?\[", re.I)
SOURCE_NAME = "sourcename"
EXTRACT_NAME = "extractname"
MATERIAL_TYPE = "materialtype"
PROTOCOL_REF = "protocolref"


def _read_rows(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def _column_key(column):
    return column.replace(" ", "").lower()


def _join(values):
    return SEPARATOR.join(dict.fromkeys(value for value in values if value))


def read_idf(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return {row[0].strip(): [value.strip() for value in row[1:]]
                for row in csv.reader(handle, delimiter="\t") if row and row[0].strip()}


def arrayexpress_series(idf):
    return {field: _join(idf.get(tag, [])) for field, tag in IDF_SERIES_ROWS}


def _protocols(idf):
    rows = (idf.get(tag, []) for tag in ("Protocol Name", "Protocol Type", "Protocol Description"))
    return {name: (kind, description) for name, kind, description in itertools.zip_longest(*rows, fillvalue="")
            if name}


def _extract_material_types(row):
    keys = [_column_key(column) for column, _ in row]
    if EXTRACT_NAME not in keys:
        return []
    start = keys.index(EXTRACT_NAME)
    return [value for key, (_, value) in zip(keys[start:], row[start:]) if key == MATERIAL_TYPE]


def arrayexpress_sample(sample, rows, idf):
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


def prepare_label_inputs(dataset_labels, out_dir, arrayexpress_dir):
    datasets = [row["dataset"] for row in _read_rows(dataset_labels)]
    wanted = set(datasets)
    samples = {dataset: {} for dataset in datasets}
    arrayexpress = set()
    for row in _read_rows(out_dir / SAMPLES_NAME):
        if row["dataset"] in wanted and row["assay"] in BULK_RNA_SEQ_ASSAYS:
            if row["sample"] in samples[row["dataset"]]:
                raise SystemExit(f"{row['dataset']}: sample {row['sample']} has two rows in the sample table")
            samples[row["dataset"]][row["sample"]] = {"sample": row["sample"], "characteristics": []}
            if row["repository"] == ARRAYEXPRESS_REPOSITORY:
                arrayexpress.add(row["dataset"])
    print(f"PROGRESS 1/5 sample table read: {sum(map(len, samples.values()))} samples of {len(datasets)} datasets",
          flush=True)

    series = {row["series"]: {field: row[field] for field in LABEL_SERIES_FIELDS}
              for row in _read_rows(out_dir / GEO_SERIES) if row["series"] in wanted}
    print(f"PROGRESS 2/5 series fields read: {len(series)} series", flush=True)

    seen = set()
    for row in _read_rows(out_dir / GEO_SAMPLES):
        sample = samples.get(row["series"], {}).get(row["sample"])
        if sample is not None:
            sample.update({name: row[column] for name, column in LABEL_SAMPLE_FIELDS})
            seen.add((row["series"], row["sample"]))
    print(f"PROGRESS 3/5 sample fields read: {len(seen)} samples", flush=True)

    for row in _read_rows(out_dir / GEO_CHARACTERISTICS):
        sample = samples.get(row["series"], {}).get(row["sample"])
        if sample is not None and row["channel"] == CHANNEL:
            sample["characteristics"].append([row["key"], row["value"]])
    print("PROGRESS 4/5 characteristics read", flush=True)

    sheet_rows = read_sheet_rows(out_dir / SHEETS_NAME, arrayexpress)
    for dataset in sorted(arrayexpress):
        idf = read_idf(find_idf(arrayexpress_dir / dataset))
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

    with gzip.open(out_dir / LABEL_INPUTS_NAME, "wt") as handle:
        for dataset in datasets:
            record = {"dataset": dataset, "series": series[dataset], "samples": list(samples[dataset].values())}
            handle.write(json.dumps(record, ensure_ascii=False) + "\n")
    print(f"{len(datasets)} datasets written to {out_dir / LABEL_INPUTS_NAME}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--downloads-dir", type=Path, required=True, help="the folder of download.sh")
    parser.add_argument("--dataset-labels", type=Path, required=True,
                        help="table with a column `dataset`: the datasets whose label inputs are written")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    harvest_geo_records(args.downloads_dir / GEO_FOLDER, args.out_dir, args.workers)
    harvest_sample_sheets(args.downloads_dir / ARRAYEXPRESS_FOLDER, args.out_dir)
    build_sample_table(args.out_dir)
    prepare_label_inputs(args.dataset_labels, args.out_dir, args.downloads_dir / ARRAYEXPRESS_FOLDER)


if __name__ == "__main__":
    main()
