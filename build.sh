#!/usr/bin/env bash
# Build round all980_v10 of the senescence benchmark from the files that download.sh fetched and the tables of
# pipeline/data/, and compare every output with the sha256 of our own build (pipeline/checksums.tsv).
#
# Usage: bash build.sh [--threads N] [--conda PROGRAM] [--envs-dir DIR] [DATA_DIR [BUILD_DIR]]
#   DATA_DIR    the folder of download.sh (default: downloads/ in this folder)
#   BUILD_DIR   where the outputs go (default: build/ in this folder); a step whose folder exists there is not run
#               again, so a build that stopped continues where it stopped
#   --threads   CPUs to use (default 8)
#   --conda     conda, mamba or micromamba, to create the environment of pipeline/environment.yml (default conda)
#   --envs-dir  where the environment is created (default BUILD_DIR/environments)
set -euo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
pipeline=$repo/pipeline
data=$pipeline/data
threads=8
conda=conda
envs=
positional=()
while [ $# -gt 0 ]; do
    case $1 in
        --threads) threads=$2; shift 2 ;;
        --conda) conda=$2; shift 2 ;;
        --envs-dir) envs=$2; shift 2 ;;
        -h|--help) sed -n '2,11p' "$0"; exit 0 ;;
        -*) echo "unknown option $1" >&2; exit 2 ;;
        *) positional+=("$1"); shift ;;
    esac
done
[ ${#positional[@]} -le 2 ] || { echo "at most two folders: DATA_DIR and BUILD_DIR" >&2; exit 2; }
downloads=$(realpath -m "${positional[0]:-$repo/downloads}")
out=$(realpath -m "${positional[1]:-$repo/build}")
envs=$(realpath -m "${envs:-$out/environments}")
[ -d "$downloads" ] || { echo "no data folder $downloads; run download.sh first" >&2; exit 1; }
mkdir -p "$out"

annotation=$downloads/geo_uniform_counts/annotation/Human.GRCh38.p13.annot.tsv.gz
round=all980_v10

log() { echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) $*"; }

# Run a step into BUILD_DIR/NAME.partial and rename it to BUILD_DIR/NAME when it succeeds; skip a step done before.
run_step() {
    local name=$1
    shift
    if [ -e "$out/$name" ]; then
        log "$name: done by an earlier run, kept"
        return
    fi
    rm -rf "$out/$name.partial"
    log "$name"
    "$@" "$out/$name.partial"
    mv "$out/$name.partial" "$out/$name"
}

# Create the environment once; a folder left by a creation that failed is created again.
prefix=$envs/count-tables
if [ ! -e "$prefix/.senbench_complete" ]; then
    rm -rf "$prefix"
    log "environment with $conda"
    CONDA_ALWAYS_YES=true MAMBA_ALWAYS_YES=true "$conda" env create --prefix "$prefix" \
        --file "$pipeline/environment.yml"
    touch "$prefix/.senbench_complete"
fi
py=$prefix/bin/python
"$py" "$pipeline/downloads.py" check --data-dir "$downloads"
log "threads $threads"

build_sample_table() {
    "$py" "$pipeline/build_sample_table.py" --downloads-dir "$downloads" --inventory "$data/inventory.tsv" \
        --dataset-labels "$data/dataset_labels.tsv" --out-dir "$1" --workers "$threads"
}

# The count tables of the datasets whose table the round reads; the records of the other datasets are added, so that
# the table of count tables lists every dataset.
build_count_tables() {
    "$py" "$pipeline/build_count_tables.py" --datasets "$data/rebuilt_datasets.tsv" \
        --uniform-state "$data/geo_uniform_counts_state.tsv" --label-inputs "$out/sample_table/label_inputs.jsonl.gz" \
        --downloads-dir "$downloads" --annotation "$annotation" --layouts "$data/count_layouts.json" \
        --records-not-rebuilt "$data/count_records_not_rebuilt.json" \
        --column-reassignments "$data/count_column_reassignments.tsv" --counted-dir "$data/raw_read_counts" \
        --out-dir "$1" --workers "$threads"
}

build_round() {
    local dir=$1
    "$py" "$pipeline/build_benchmark_bulk.py" --sample-labels "$data/sample_labels.tsv.gz" \
        --dataset-labels "$data/dataset_labels.tsv" --samples "$out/sample_table/samples.tsv.gz" \
        --count-tables "$out/counts/count_tables.tsv" --count-columns "$out/counts/count_table_columns.tsv.gz" \
        --counts-dir "$out/counts" --annotation "$annotation" --min-genes 10000 --max-median-genes-with-reads 35000 \
        --out-dir "$dir"
    "$py" "$pipeline/build_study_split.py" --samples "$dir/benchmark_bulk_samples.tsv" \
        --contacts "$data/geo_series_contacts.tsv" --lists "$data/method_source_datasets.tsv" \
        --main-scores "$data/main_scores.tsv" --accession-pubmed-ids "$data/accession_pubmed_ids.tsv" \
        --laboratory-group-links "$data/laboratory_group_links.tsv" \
        --split-amendments "$data/split_amendments.tsv" --out-dir "$dir/split" --classes-dir "$dir"
    "$py" "$pipeline/write_final_sets.py" --split-dir "$dir/split" --samples "$dir/benchmark_bulk_samples.tsv" \
        --label-corrections "$data/label_corrections.tsv"
}

build_agentomics_dataset() {
    local dir=$1
    mkdir -p "$dir"
    "$py" "$pipeline/build_agentomics_dataset.py" --counts "$out/$round/benchmark_bulk_counts.tsv.gz" \
        --training "$out/$round/split/training_data_samples.tsv" --test "$out/$round/split/test_samples.tsv" \
        --annotation "$annotation" --validation-share 0.5 --dataset-dir "$dir/senescence_bulk_$round" \
        --key-dir "$dir/senescence_bulk_${round}_key"
}

run_step sample_table build_sample_table
run_step counts build_count_tables
run_step "$round" build_round
run_step agentomics build_agentomics_dataset

log "comparison with our build"
"$py" "$pipeline/compare_checksums.py" --checksums "$pipeline/checksums.tsv" --build-dir "$out" \
    --out "$out/checksum_comparison.tsv"
