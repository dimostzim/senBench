#!/usr/bin/env bash
# Download every public file that the build of round all980_v10 needs (the list data/downloads.tsv): GEO's series and
# sample records, GEO's uniform count tables, the authors' supplementary count files, the ArrayExpress files, the
# reference that raw reads are aligned to, and the raw reads of the datasets counted from reads.
#
# Usage: download.sh [--no-raw-reads] [--workers N] [DATA_DIR]
#   DATA_DIR        where the files go (default ./downloads). A file recorded there as complete is not downloaded
#                   again, and a cut transfer continues from the bytes on disk where the source states the size.
#   --no-raw-reads  leave out the raw reads and the reference they are aligned to
#   --workers       files downloaded at a time (default 3)
# Needs python3 and curl; the shell's proxy variables apply. Writes DATA_DIR/download_summary.tsv (every file: whether
# it was downloaded, and whether its md5 equals that of the file of our build) and exits with 1 if a file is missing.
set -euo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
workers=3
listing=()
positional=()
while [ $# -gt 0 ]; do
    case $1 in
        --no-raw-reads) listing+=(--no-raw-reads); shift ;;
        --workers) workers=$2; shift 2 ;;
        -h|--help) sed -n '2,12p' "$0"; exit 0 ;;
        -*) echo "unknown option $1" >&2; exit 2 ;;
        *) positional+=("$1"); shift ;;
    esac
done
[ ${#positional[@]} -le 1 ] || { echo "one folder at most: DATA_DIR" >&2; exit 2; }
data_dir=$(realpath -m "${positional[0]:-downloads}")
mkdir -p "$data_dir"
downloads=(python3 "$repo/pipeline/downloads.py" --downloads "$repo/data/downloads.tsv" --data-dir "$data_dir")

"${downloads[@]}" list "${listing[@]}"
status=0
for kind in geo geo_uniform_counts biostudies_arrayexpress reference raw_reads; do
    manifest=$data_dir/manifests/$kind.tsv
    [ -e "$manifest" ] || continue
    echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) $kind"
    mkdir -p "$data_dir/manifests/$kind"
    # One second after each file, so that no server sees more than a few requests a second.
    python3 "$repo/pipeline/download_manifest.py" --manifest "$manifest" --out "$data_dir/$kind" \
        --units-dir "$data_dir/manifests/$kind" --min-free-gb 10 --workers "$workers" --unit file --retries 8 \
        --retry-delay 60 --delay 1 || status=1
done
"${downloads[@]}" summarise || status=1
exit $status
