#!/usr/bin/env python3
"""Compare the outputs of a build with the sha256 of our own build (checksums.tsv) and write one row per output:
`equal`, `differs` or `missing`. A gzip file is compared by its uncompressed content, because gzip stores the time of
writing. Exits with 1 if an output differs or is missing."""

import argparse
import csv
import gzip
import hashlib
import sys
from collections import Counter
from pathlib import Path

UNCOMPRESSED = "uncompressed content"
EQUAL = "equal"
DIFFERS = "differs"
MISSING = "missing"
CHUNK_BYTES = 1 << 20


def compute_sha256(path, hashed):
    """Return the sha256 of a file, or of its uncompressed content."""
    digest = hashlib.sha256()
    with (gzip.open if hashed == UNCOMPRESSED else open)(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def compare_outputs(checksums, build_dir):
    """Return the rows of the checksum table with the sha256 of the build and the outcome."""
    rows = []
    for row in checksums:
        path = build_dir / row["path"]
        built = compute_sha256(path, row["hashed"]) if path.exists() else ""
        outcome = MISSING if not built else EQUAL if built == row["sha256"] else DIFFERS
        rows.append({**row, "sha256_of_this_build": built, "outcome": outcome})
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checksums", type=Path, required=True, help="checksums.tsv")
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    arguments = parser.parse_args()
    with arguments.checksums.open(newline="") as handle:
        checksums = list(csv.DictReader(handle, delimiter="\t"))
    rows = compare_outputs(checksums, arguments.build_dir)
    with arguments.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    by_folder = Counter((row["path"].split("/", 1)[0], row["outcome"]) for row in rows)
    for (folder, outcome), count in sorted(by_folder.items()):
        print(f"{folder}: {count} {outcome}")
    for row in rows:
        if row["outcome"] != EQUAL:
            print(f"{row['outcome']}: {row['path']}")
    outcomes = Counter(row["outcome"] for row in rows)
    print(f"{outcomes[EQUAL]} of {len(rows)} outputs equal those of our build; {outcomes[DIFFERS]} differ, "
          f"{outcomes[MISSING]} missing; the table is {arguments.out}")
    if outcomes[DIFFERS] or outcomes[MISSING]:
        sys.exit(1)


if __name__ == "__main__":
    main()
