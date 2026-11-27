#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/pipeline/env.sh"
pipeline=$repo/pipeline
data=$pipeline/data
downloads=$repo/downloads
out=$repo/build
round=$out/all980_v10
annotation=$downloads/geo_uniform_counts/annotation/Human.GRCh38.p13.annot.tsv.gz
python=$envs/pipeline/bin/python
workers=$(nproc)

run_step() {
    local name=$1
    shift
    [ ! -e "$out/$name" ] || return 0
    rm -rf "$out/$name.partial"
    echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) $name"
    "$@" "$out/$name.partial"
    mv "$out/$name.partial" "$out/$name"
}

build_sample_table() {
    "$python" "$pipeline/build_sample_table.py" --downloads-dir "$downloads" \
        --dataset-labels "$data/dataset_labels.tsv" --out-dir "$1" --workers "$workers"
}

build_count_tables() {
    "$python" "$pipeline/build_count_tables.py" --datasets "$data/rebuilt_datasets.tsv" \
        --uniform-state "$data/geo_uniform_counts_state.tsv" --label-inputs "$out/sample_table/label_inputs.jsonl.gz" \
        --downloads-dir "$downloads" --annotation "$annotation" --layouts "$data/count_layouts.json" \
        --records-not-rebuilt "$data/count_records_not_rebuilt.json" \
        --column-reassignments "$data/count_column_reassignments.tsv" --counted-dir "$data/raw_read_counts" \
        --out-dir "$1" --workers "$workers"
}

build_round() {
    "$python" "$pipeline/build_benchmark_bulk.py" --sample-labels "$data/sample_labels.tsv.gz" \
        --dataset-labels "$data/dataset_labels.tsv" --samples "$out/sample_table/samples.tsv.gz" \
        --count-tables "$out/counts/count_tables.tsv" --count-columns "$out/counts/count_table_columns.tsv.gz" \
        --counts-dir "$out/counts" --annotation "$annotation" --min-genes 10000 --max-median-genes-with-reads 35000 \
        --out-dir "$1"
    "$python" "$pipeline/build_study_split.py" --samples "$1/benchmark_bulk_samples.tsv" \
        --contacts "$data/geo_series_contacts.tsv" --lists "$data/method_source_datasets.tsv" \
        --list-methods "$data/list_methods.tsv" --accession-pubmed-ids "$data/accession_pubmed_ids.tsv" \
        --laboratory-group-links "$data/laboratory_group_links.tsv" --out-dir "$1/split"
    "$python" "$pipeline/write_final_sets.py" --split-dir "$1/split" --samples "$1/benchmark_bulk_samples.tsv" \
        --label-corrections "$data/label_corrections.tsv"
}

build_agentomics_dataset() {
    "$python" "$pipeline/build_agentomics_dataset.py" --counts "$round/benchmark_bulk_counts.tsv.gz" \
        --training "$round/split/training_data_samples.tsv" --test "$round/split/test_samples.tsv" \
        --annotation "$annotation" --validation-share 0.5 --dataset-dir "$1/senescence_bulk_all980_v10" \
        --key-dir "$1/senescence_bulk_all980_v10_key"
}

python3 "$pipeline/downloads.py" check
install pipeline micromamba env create --yes --prefix "$envs/pipeline" --file "$pipeline/environment.yml"
run_step sample_table build_sample_table
run_step counts build_count_tables
run_step all980_v10 build_round
run_step agentomics build_agentomics_dataset
"$python" "$pipeline/compare_checksums.py" --checksums "$pipeline/checksums.tsv" --build-dir "$out" \
    --out "$out/checksum_comparison.tsv"
