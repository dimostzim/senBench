#!/usr/bin/env python3
"""Download the files that the benchmark is built from (data/downloads.tsv) into a data folder, or check that a data
folder holds them all.

One row per file: `kind` is the folder of the data folder that holds it, `accession` the folder inside that, `file` its
path inside the accession folder; `bytes` is the size that the source states (empty where it states none); `our_bytes`
and `our_md5` describe the file that our own build used.

`download` fetches every listed file that the data folder does not hold yet, a few at a time, continues a cut transfer
from the bytes on disk where the source states the size, and stops when the disk has less than 10 GB free. It writes
download_summary.tsv (every file: whether it was downloaded, and whether its md5 equals that of the file of our build)
and exits with 1 if a file is missing. `check` exits with 1 unless the data folder holds every listed file.
"""

import argparse
import csv
import hashlib
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

DOWNLOADS = Path(__file__).resolve().parent / "data" / "downloads.tsv"
# The files downloaded so far, one row per file, appended as each file is complete; the last row of a file counts.
RECORD_NAME = "download_record.tsv"
RECORD_COLUMNS = ("kind", "accession", "file", "bytes", "md5")
SUMMARY_NAME = "download_summary.tsv"
SUMMARY_COLUMNS = ("kind", "accession", "file", "status", "bytes", "md5", "same_as_ours")
GIB = 2**30
CLIENT = "senBench-download"
CURL_CANNOT_RESUME = 33
# The HTTP status of a server that refuses a client for asking too often.
TOO_MANY_REQUESTS = "429"
# Attempts in a row without new bytes before a file counts as failed, and the waits after a refused or cut attempt.
RETRIES = 8
RETRY_DELAY_SECONDS = 60
RATE_LIMIT_DELAY_SECONDS = 300
# The pause after each file, so that no server sees more than a few requests a second.
DELAY_SECONDS = 1
MIN_FREE_BYTES = 10 * 10**9
ERROR_DOCUMENT_SUFFIX = ".error_document"
# The EBI file servers write an error document into the body of a file where they cannot finish a transfer.
SERVER_ERROR_START = b"<Error><Code>"
SERVER_ERROR_MARKER = b"<ErrorMessage/><RequestId/></Error>"
SERVER_ERROR_MAX_BYTES = 1024
CHUNK_BYTES = 1 << 20

csv.field_size_limit(sys.maxsize)

_write_lock = threading.Lock()
_free_space_reached = threading.Event()


class _FileProblem(Exception):
    """A file that could not be downloaded or does not match its row."""


def read_rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path, columns, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def file_path(data_dir, row):
    """Return where a listed file is in the data folder."""
    return data_dir / row["kind"] / row["accession"] / row["file"]


def describe_size(rows):
    return f"{len(rows)} files, {sum(int(row['our_bytes'] or 0) for row in rows) / GIB:,.2f} GiB"


def _digest(path):
    """Return the md5 of a file and whether the file holds a server's error document."""
    digest = hashlib.md5()
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


def _fetch_until_complete(url, target, size):
    """Download a file, continuing after every cut transfer; give up after RETRIES attempts in a row without new bytes.
    After an attempt that the server refused with HTTP 429, wait RATE_LIMIT_DELAY_SECONDS instead of
    RETRY_DELAY_SECONDS."""
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
        if without_progress > RETRIES:
            expected = "" if size is None else f", {size} expected"
            raise _FileProblem(f"curl exit {code}, HTTP {status}; {have} bytes on disk{expected}")
        time.sleep(RATE_LIMIT_DELAY_SECONDS if status == TOO_MANY_REQUESTS else RETRY_DELAY_SECONDS)


def _download_file(row, target):
    """Download one file and check it; return the md5 of the file on disk."""
    size = int(row["bytes"]) if row["bytes"] else None
    if size is None or not target.exists() or target.stat().st_size != size:
        if size is None and target.exists():
            # Without a known size, a part of a file cannot be told from the whole file.
            target.unlink()
        _fetch_until_complete(row["url"], target, size)
    on_disk = target.stat().st_size
    if on_disk == 0 and size != 0:
        # A server can answer with success and an empty body for a file that it fails to read.
        target.unlink()
        expected = "the list states no size" if size is None else f"{size} bytes expected"
        raise _FileProblem(f"the server sent an empty file; {expected}")
    if size is not None and on_disk != size:
        raise _FileProblem(f"{on_disk} bytes on disk, {size} expected")
    md5, holds_error = _digest(target)
    if holds_error:
        # Kept under another name, so that the next attempt downloads the file from the start.
        target.replace(target.with_name(target.name + ERROR_DOCUMENT_SUFFIX))
        raise _FileProblem("the file holds an error document of the server")
    return md5


def download_files(rows, data_dir, workers):
    """Download every listed file that the record of the data folder does not name as complete; return the problems
    of the files that failed."""
    record = data_dir / RECORD_NAME
    if not record.exists():
        _append_row(record, RECORD_COLUMNS)
    recorded = {(got["kind"], got["accession"], got["file"]) for got in read_rows(record)}

    def download_row(row):
        target = file_path(data_dir, row)
        if (row["kind"], row["accession"], row["file"]) in recorded and target.exists():
            return "skipped", "downloaded by an earlier run"
        if _free_space_reached.is_set() or shutil.disk_usage(data_dir).free - int(row["bytes"] or 0) < MIN_FREE_BYTES:
            _free_space_reached.set()
            return "failed", "stopped by the free-space limit"
        target.parent.mkdir(parents=True, exist_ok=True)
        try:
            md5 = _download_file(row, target)
        except _FileProblem as problem:
            return "failed", str(problem)
        finally:
            time.sleep(DELAY_SECONDS)
        _append_row(record, (row["kind"], row["accession"], row["file"], target.stat().st_size, md5))
        return "ok", ""

    failed = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for done, (row, (status, note)) in enumerate(zip(rows, pool.map(download_row, rows)), start=1):
            name = f"{row['kind']}/{row['accession']}/{row['file']}"
            print(f"PROGRESS {done}/{len(rows)} {name} {status}{': ' + note if note else ''}", flush=True)
            if status == "failed":
                failed.append(f"{name}: {note}")
    return failed


def summarise_downloads(rows, data_dir):
    """Return one summary row per listed file."""
    record = data_dir / RECORD_NAME
    recorded = {(got["kind"], got["accession"], got["file"]): got for got in read_rows(record)}
    summary = []
    for row in rows:
        got = recorded.get((row["kind"], row["accession"], row["file"]))
        present = got is not None and file_path(data_dir, row).exists()
        md5 = got["md5"] if present else ""
        same = "unknown" if not (md5 and row["our_md5"]) else "yes" if md5 == row["our_md5"] else "no"
        summary.append({"kind": row["kind"], "accession": row["accession"], "file": row["file"],
                        "status": "ok" if present else "missing", "bytes": got["bytes"] if present else "",
                        "md5": md5, "same_as_ours": same})
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=("download", "check"))
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=3, help="files downloaded at a time")
    arguments = parser.parse_args()
    rows = read_rows(DOWNLOADS)
    data_dir = arguments.data_dir

    if arguments.command == "check":
        missing = [row for row in rows if not file_path(data_dir, row).exists()]
        if missing:
            sys.exit(f"{data_dir} lacks {len(missing)} files that the build needs, the first: "
                     f"{file_path(data_dir, missing[0])}; run download.sh")
        return
    data_dir.mkdir(parents=True, exist_ok=True)
    print(f"to download: {describe_size(rows)}, measured on our own downloads", flush=True)
    failed = download_files(rows, data_dir, arguments.workers)
    summary = summarise_downloads(rows, data_dir)
    write_rows(data_dir / SUMMARY_NAME, SUMMARY_COLUMNS, summary)
    missing = [row for row in summary if row["status"] == "missing"]
    differ = [row for row in summary if row["same_as_ours"] == "no"]
    print(f"{len(summary) - len(missing)} of {len(summary)} files downloaded; {len(differ)} differ from the files of "
          f"our build (by md5); the table is {data_dir / SUMMARY_NAME}")
    for row in differ:
        print(f"differs from ours: {row['kind']}/{row['accession']}/{row['file']}")
    for problem in failed:
        print(f"failed: {problem}")
    if missing:
        sys.exit(f"{len(missing)} files missing, the first: {missing[0]['kind']}/{missing[0]['accession']}/"
                 f"{missing[0]['file']}")


if __name__ == "__main__":
    main()
