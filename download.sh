#!/usr/bin/env bash
# Download every public file that build.sh needs (the list pipeline/data/downloads.tsv, 18.76 GiB): GEO's series records
# and the authors' supplementary count files, GEO's uniform count tables and its gene annotation table, and the files of
# the ArrayExpress experiments.
#
# Usage: bash download.sh [--workers N] [DATA_DIR]
#   DATA_DIR   where the files go (default: downloads/ in this folder). A file recorded there as complete is not
#              downloaded again, and a cut transfer continues from the bytes on disk where the source states the size.
#   --workers  files downloaded at a time (default 3)
# Needs python3 and curl; the shell's proxy variables apply. Writes DATA_DIR/download_summary.tsv (every file: whether
# it was downloaded, and whether its md5 equals that of the file of our build) and exits with 1 if a file is missing.
set -euo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
workers=3
positional=()
while [ $# -gt 0 ]; do
    case $1 in
        --workers) workers=$2; shift 2 ;;
        -h|--help) sed -n '2,11p' "$0"; exit 0 ;;
        -*) echo "unknown option $1" >&2; exit 2 ;;
        *) positional+=("$1"); shift ;;
    esac
done
[ ${#positional[@]} -le 1 ] || { echo "one folder at most: DATA_DIR" >&2; exit 2; }
data_dir=$(realpath -m "${positional[0]:-$repo/downloads}")
exec python3 "$repo/pipeline/downloads.py" download --workers "$workers" --data-dir "$data_dir"
