#!/usr/bin/env python3
"""Download the files of file manifests, one folder per dataset, and check each file's size and checksum."""

import argparse
import csv
import hashlib
import re
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

DOWNLOADS_NAME = "downloads.tsv"
UNITS_NAME = "units.tsv"
DOWNLOADS_COLUMNS = ("accession", "file", "bytes", "checksum_type", "checksum", "url", "status")
DEFAULT_CHECKSUM = "md5"
CLIENT = "senBench-download"
CURL_CANNOT_RESUME = 33
# The HTTP status of a server that refuses a client for asking too often.
TOO_MANY_REQUESTS = "429"
MISMATCH_SUFFIX = ".checksum_mismatch"
# The EBI file servers write an error document into the body of a file where they cannot finish a transfer.
SERVER_ERROR_START = b"<Error><Code>"
SERVER_ERROR_MARKER = b"<ErrorMessage/><RequestId/></Error>"
SERVER_ERROR_MAX_BYTES = 1024
CHUNK_BYTES = 1 << 20

csv.field_size_limit(sys.maxsize)

_write_lock = threading.Lock()
_free_space_reached = threading.Event()


class _FileProblem(Exception):
    """A file that could not be downloaded or does not match its manifest row."""


def dataset_folder(accession):
    """Return the folder name of a dataset; a DOI holds characters that a folder name should not."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", accession)


def _digest(path, kind):
    """Return the checksum of a file and whether the file holds a server's error document."""
    digest = hashlib.new(kind)
    holds_error = False
    tail = b""
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_BYTES), b""):
            digest.update(chunk)
            holds_error = holds_error or SERVER_ERROR_MARKER in tail + chunk
            tail = chunk[-len(SERVER_ERROR_MARKER):]
    return digest.hexdigest(), holds_error


def _append_row(path, row):
    with _write_lock, path.open("a", newline="") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerow(row)


def _fetch(url, target, resume):
    """Return curl's exit code and the HTTP status of the last response."""
    # No retry inside curl: after a cut transfer its retry discards what the attempt had written and starts again.
    # A transfer slower than 1 kB/s for two minutes is cut off, so a stalled connection cannot hang the run.
    command = ["curl", "-fsSL", "-A", CLIENT, "--speed-limit", "1024", "--speed-time", "120", "-w", "%{http_code}"]
    result = subprocess.run(command + (["-C", "-"] if resume else []) + ["-o", str(target), url],
                            stdout=subprocess.PIPE, text=True)
    return result.returncode, result.stdout.strip()


def _remove_trailing_error_document(target):
    """Cut off an error document at the end of a file; return whether there was one."""
    with target.open("rb+") as handle:
        size = handle.seek(0, 2)
        start = max(size - SERVER_ERROR_MAX_BYTES, 0)
        handle.seek(start)
        tail = handle.read()
        if not tail.endswith(SERVER_ERROR_MARKER) or SERVER_ERROR_START not in tail:
            return False
        handle.truncate(start + tail.rindex(SERVER_ERROR_START))
        return True


def _fetch_until_complete(url, target, size, retries, retry_delay, rate_limit_delay):
    """Download a file, continuing after every cut transfer; give up after `retries` attempts in a row without new bytes.
    After an attempt that the server refused with HTTP 429, wait `rate_limit_delay` seconds instead of `retry_delay`."""
    without_progress = 0
    while True:
        before = target.stat().st_size if target.exists() else 0
        code, status = _fetch(url, target, True)
        if code == CURL_CANNOT_RESUME:
            target.unlink()
            code, status = _fetch(url, target, False)
        # Continuing after the error document would leave it inside the file.
        cut = target.exists() and _remove_trailing_error_document(target)
        have = target.stat().st_size if target.exists() else 0
        # A server can end a transfer early and still report success, so the stated size decides when one is known.
        if code == 0 and not cut and (size is None or have >= size):
            return
        without_progress = 0 if have > before else without_progress + 1
        if without_progress > retries:
            expected = "" if size is None else f", {size} expected"
            raise _FileProblem(f"curl exit {code}, HTTP {status}; {have} bytes on disk{expected}")
        time.sleep(rate_limit_delay if status == TOO_MANY_REQUESTS else retry_delay)


def _line_ends(path):
    count = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_BYTES), b""):
            count += chunk.count(b"\n")
    return count


def _download_file(row, target, retries, retry_delay, rate_limit_delay, larger_by_line_count):
    """Download one file and check it; return the checksum type and the checksum of the file on disk."""
    stated = int(row["bytes"]) if row["bytes"] else None
    # The size to download up to, which a stated size that may be larger by one byte per line does not give.
    size = None if larger_by_line_count else stated
    if size is None or not target.exists() or target.stat().st_size != size:
        if size is None and target.exists():
            # Without a known size, a part of a file cannot be told from the whole file.
            target.unlink()
        _fetch_until_complete(row["url"], target, size, retries, retry_delay, rate_limit_delay)
    on_disk = target.stat().st_size
    if on_disk == 0 and stated != 0:
        # A server can answer with success and an empty body for a file that it fails to read.
        target.unlink()
        expected = "the list states no size" if stated is None else f"{stated} bytes expected"
        raise _FileProblem(f"the server sent an empty file; {expected}")
    if stated is not None and on_disk != stated and not (larger_by_line_count
                                                         and stated - on_disk == _line_ends(target)):
        raise _FileProblem(f"{on_disk} bytes on disk, {stated} expected")
    kind = row["checksum_type"] or DEFAULT_CHECKSUM
    checksum, holds_error = _digest(target, kind)
    if holds_error or (row["checksum"] and checksum != row["checksum"].lower()):
        # Kept under another name, so that the next attempt downloads the file from the start.
        target.replace(target.with_name(target.name + MISMATCH_SUFFIX))
        if holds_error:
            raise _FileProblem("the file holds an error document of the server")
        raise _FileProblem(f"{kind} is {checksum}, {row['checksum']} expected")
    return kind, checksum


def _download_unit(rows, out_dir, min_free_bytes, retries, retry_delay, rate_limit_delay, delay, larger_by_line_count,
                   downloads, recorded):
    """Download the files of one unit; return (status, note) for the units table."""
    done = 0
    problems = []
    for row in rows:
        accession = row["accession"]
        target = out_dir / dataset_folder(accession) / row["file"]
        if (accession, row["file"]) in recorded and target.exists():
            done += 1
            continue
        if _free_space_reached.is_set() or shutil.disk_usage(out_dir).free - int(row["bytes"] or 0) < min_free_bytes:
            _free_space_reached.set()
            problems.append("stopped by the free-space limit")
            break
        target.parent.mkdir(parents=True, exist_ok=True)
        try:
            kind, checksum = _download_file(row, target, retries, retry_delay, rate_limit_delay,
                                            any(re.search(pattern, row["url"]) for pattern in larger_by_line_count))
        except _FileProblem as problem:
            problems.append(f"{row['file']}: {problem}")
            continue
        finally:
            time.sleep(delay)
        _append_row(downloads, (accession, row["file"], target.stat().st_size, kind, checksum, row["url"], "ok"))
        done += 1
    note = f"{done} of {len(rows)} files"
    if problems:
        return "failed", f"{note}; {len(problems)} problems, the first: {problems[0]}"
    return "ok", note


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, nargs="+", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--units-dir", type=Path, required=True)
    parser.add_argument("--min-free-gb", type=float, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--retries", type=int, default=8)
    # A server that refuses requests while it is busy needs minutes before the next attempt on the same file.
    parser.add_argument("--retry-delay", type=float, default=10)
    parser.add_argument("--rate-limit-delay", type=float, default=300)
    # Some servers refuse a client that asks for files in quick succession; the pause follows every request.
    parser.add_argument("--delay", type=float, default=0)
    # Some repositories state a size that is one byte per line larger than the file they serve; a pattern names
    # their addresses.
    parser.add_argument("--stated-size-larger-by-line-count", nargs="*", default=[])
    # A unit is downloaded by one worker, file after file; "file" lets the files of one dataset download in parallel.
    parser.add_argument("--unit", choices=("dataset", "file"), default="dataset")
    args = parser.parse_args()

    units_rows = {}
    for path in args.manifest:
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if not row["url"] or row["access"] != "open":
                    raise SystemExit(f"{path}: {row['accession']} {row['file']} has no open address")
                name = row["accession"] if args.unit == "dataset" else f"{row['accession']}/{row['file']}"
                units_rows.setdefault(name, []).append(row)

    args.out.mkdir(parents=True, exist_ok=True)
    downloads = args.out / DOWNLOADS_NAME
    if not downloads.exists():
        _append_row(downloads, DOWNLOADS_COLUMNS)
    with downloads.open(newline="") as handle:
        recorded = {(row["accession"], row["file"]) for row in csv.DictReader(handle, delimiter="\t")}
    units = args.units_dir / UNITS_NAME
    units.write_text("")
    min_free_bytes = int(args.min_free_gb * 1e9)

    def run_unit(name):
        status, note = _download_unit(units_rows[name], args.out, min_free_bytes, args.retries, args.retry_delay,
                                      args.rate_limit_delay, args.delay, args.stated_size_larger_by_line_count,
                                      downloads, recorded)
        _append_row(units, (name, status, note))
        return name, status, note

    failed = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for done, (name, status, note) in enumerate(pool.map(run_unit, units_rows), start=1):
            print(f"PROGRESS {done}/{len(units_rows)} {name} {status}: {note}", flush=True)
            if status != "ok":
                failed.append(name)
    if failed:
        print(f"failed units ({len(failed)}): {' '.join(failed)}", flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
