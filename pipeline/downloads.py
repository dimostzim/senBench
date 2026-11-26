#!/usr/bin/env python3
"""The list of files that the benchmark is built from (data/downloads.tsv): which to download, what was downloaded,
and whether a data folder holds what a build needs.

One row per file: `kind` is the folder of the data folder that holds it, `accession` the folder inside that, `file` its
path inside the accession folder; `bytes`, `checksum_type` and `checksum` are what the source states (empty where it
states nothing); `our_bytes` and `our_md5` describe the file that our own build used.
"""

import argparse
import csv
import sys
from pathlib import Path

GEO = "geo"
GEO_UNIFORM_COUNTS = "geo_uniform_counts"
ARRAYEXPRESS = "biostudies_arrayexpress"
REFERENCE = "reference"
RAW_READS = "raw_reads"
KINDS = (GEO, GEO_UNIFORM_COUNTS, ARRAYEXPRESS, REFERENCE, RAW_READS)
# Only the counting of raw reads needs the reads and the reference they are aligned to.
COUNTING_KINDS = (REFERENCE, RAW_READS)
MANIFESTS_FOLDER = "manifests"
# The record that download_manifest.py keeps in each folder it downloads into.
RECORD_NAME = "downloads.tsv"
SUMMARY_NAME = "download_summary.tsv"
SUMMARY_COLUMNS = ("kind", "accession", "file", "status", "bytes", "md5", "same_as_ours")
GIB = 2**30

csv.field_size_limit(sys.maxsize)


def read_rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path, columns, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def file_path(data_dir, row):
    """Return where a listed file is in the data folder, as download_manifest.py places it."""
    return data_dir / row["kind"] / row["accession"] / row["file"]


def select_kinds(without_raw_reads):
    return [kind for kind in KINDS if not (without_raw_reads and kind in COUNTING_KINDS)]


def describe_size(rows):
    return f"{len(rows)} files, {sum(int(row['our_bytes'] or 0) for row in rows) / GIB:,.2f} GiB"


def write_manifests(rows, data_dir, kinds):
    """Write one manifest per kind for download_manifest.py and remove that of a kind left out; return the rows per
    kind."""
    for kind in set(KINDS) - set(kinds):
        (data_dir / MANIFESTS_FOLDER / f"{kind}.tsv").unlink(missing_ok=True)
    by_kind = {kind: [row for row in rows if row["kind"] == kind] for kind in kinds}
    for kind, kind_rows in by_kind.items():
        write_rows(data_dir / MANIFESTS_FOLDER / f"{kind}.tsv", list(rows[0]), kind_rows)
    return by_kind


def summarise_downloads(rows, data_dir):
    """Return one summary row per listed file of the kinds that have a manifest in the data folder."""
    kinds = [kind for kind in KINDS if (data_dir / MANIFESTS_FOLDER / f"{kind}.tsv").exists()]
    recorded = {}
    for kind in kinds:
        record = data_dir / kind / RECORD_NAME
        for got in read_rows(record) if record.exists() else []:
            # The last row of a file counts.
            recorded[(kind, got["accession"], got["file"])] = got
    summary = []
    for row in rows:
        if row["kind"] not in kinds:
            continue
        got = recorded.get((row["kind"], row["accession"], row["file"]))
        present = got is not None and file_path(data_dir, row).exists()
        md5 = got["checksum"] if present and got["checksum_type"] == "md5" else ""
        same = "unknown" if not (md5 and row["our_md5"]) else "yes" if md5 == row["our_md5"] else "no"
        summary.append({"kind": row["kind"], "accession": row["accession"], "file": row["file"],
                        "status": "ok" if present else "missing", "bytes": got["bytes"] if present else "",
                        "md5": md5, "same_as_ours": same})
    return summary


def list_missing(rows, data_dir, count_datasets):
    """Return the listed files that a build needs and the data folder lacks: every file of the kinds that are not
    for counting, and with datasets to count, the reference and the reads of those datasets."""
    needed = [row for row in rows if row["kind"] not in COUNTING_KINDS
              or (count_datasets and row["kind"] == REFERENCE)
              or (row["kind"] == RAW_READS and row["accession"] in count_datasets)]
    return [row for row in needed if not file_path(data_dir, row).exists()]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--downloads", type=Path, required=True, help="data/downloads.tsv")
    parser.add_argument("--data-dir", type=Path, required=True)
    commands = parser.add_subparsers(dest="command", required=True)
    listing = commands.add_parser("list", help="write the manifests of the files to download and print their size")
    listing.add_argument("--no-raw-reads", action="store_true")
    commands.add_parser("summarise", help=f"write {SUMMARY_NAME}; fails if a file of a listed kind is missing")
    check = commands.add_parser("check", help="fail unless the data folder holds every file a build needs")
    check.add_argument("--count-datasets", nargs="*", default=[], help="datasets whose raw reads will be counted")
    arguments = parser.parse_args()
    rows = read_rows(arguments.downloads)
    data_dir = arguments.data_dir

    if arguments.command == "list":
        kinds = select_kinds(arguments.no_raw_reads)
        by_kind = write_manifests(rows, data_dir, kinds)
        for kind, kind_rows in by_kind.items():
            print(f"{kind}: {describe_size(kind_rows)}")
        print(f"to download: {describe_size([row for kind_rows in by_kind.values() for row in kind_rows])}, "
              f"measured on our own downloads")
        left_out = [row for row in rows if row["kind"] not in kinds]
        if left_out:
            datasets = sorted({row["accession"] for row in left_out if row["kind"] == RAW_READS})
            print(f"--no-raw-reads: not downloading the raw reads of {len(datasets)} datasets and the reference they "
                  f"are aligned to ({describe_size(left_out)}); build.sh can then build their count tables only "
                  f"from tables given with --counted-reads")
    elif arguments.command == "summarise":
        summary = summarise_downloads(rows, data_dir)
        write_rows(data_dir / SUMMARY_NAME, SUMMARY_COLUMNS, summary)
        missing = [row for row in summary if row["status"] == "missing"]
        differ = [row for row in summary if row["same_as_ours"] == "no"]
        print(f"{len(summary) - len(missing)} of {len(summary)} files downloaded; {len(differ)} differ from the "
              f"files of our build (by md5); the table is {data_dir / SUMMARY_NAME}")
        for row in differ:
            print(f"differs from ours: {row['kind']}/{row['accession']}/{row['file']}")
        if missing:
            sys.exit(f"{len(missing)} files missing, the first: {missing[0]['kind']}/{missing[0]['accession']}/"
                     f"{missing[0]['file']}")
    else:
        missing = list_missing(rows, data_dir, set(arguments.count_datasets))
        if missing:
            sys.exit(f"{data_dir} lacks {len(missing)} files that the build needs, the first: "
                     f"{file_path(data_dir, missing[0])}; run download.sh")


if __name__ == "__main__":
    main()
