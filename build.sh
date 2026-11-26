#!/usr/bin/env bash
# Build round all980_v10 of the senescence benchmark from the files that download.sh fetched, and compare every output
# with the sha256 of our own build (checksums/all980_v10.tsv).
#
# Usage: build.sh [--threads N] [--conda PROGRAM] [--envs-dir DIR] [--counted-reads DIR] [DATA_DIR [BUILD_DIR]]
#   DATA_DIR       the folder of download.sh (default ./downloads)
#   BUILD_DIR      where the outputs go (default ./build); a step whose folder exists there is not run again, so a
#                  build that stopped continues where it stopped
#   --threads      CPUs to use (default 8)
#   --conda        conda, mamba or micromamba, to create the environments of pipeline/environments (default conda)
#   --envs-dir     where the environments are created (default BUILD_DIR/environments)
#   --counted-reads  a folder of count tables of raw reads, <dataset>.tsv.gz as pipeline/count_raw_reads.py writes
#                  them; a dataset with a table there is not counted from its reads
set -euo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
pipeline=$repo/pipeline
data=$repo/data
threads=8
conda=conda
envs=
counted_reads=
positional=()
while [ $# -gt 0 ]; do
    case $1 in
        --threads) threads=$2; shift 2 ;;
        --conda) conda=$2; shift 2 ;;
        --envs-dir) envs=$2; shift 2 ;;
        --counted-reads) counted_reads=$(realpath "$2"); shift 2 ;;
        -h|--help) sed -n '2,13p' "$0"; exit 0 ;;
        -*) echo "unknown option $1" >&2; exit 2 ;;
        *) positional+=("$1"); shift ;;
    esac
done
[ ${#positional[@]} -le 2 ] || { echo "at most two folders: DATA_DIR and BUILD_DIR" >&2; exit 2; }
downloads=$(realpath "${positional[0]:-downloads}")
out=$(realpath -m "${positional[1]:-build}")
envs=$(realpath -m "${envs:-$out/environments}")
[ -d "$downloads" ] || { echo "no data folder $downloads; run download.sh first" >&2; exit 1; }
mkdir -p "$out"

annotation=$downloads/geo_uniform_counts/annotation/Human.GRCh38.p13.annot.tsv.gz
round=all980_v10
# featureCounts' option set of our counts: every alignment of a read counted, a share of a read on each of its genes.
option_set=multimapping_overlapping_fraction

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

# Create an environment of pipeline/environments once; a folder left by a creation that failed is created again.
make_environment() {
    local name=$1 prefix=$envs/$1
    [ -e "$prefix/.senbench_complete" ] && return
    rm -rf "$prefix"
    log "environment $name with $conda"
    CONDA_ALWAYS_YES=true MAMBA_ALWAYS_YES=true "$conda" env create --prefix "$prefix" \
        --file "$pipeline/environments/$name.yml"
    touch "$prefix/.senbench_complete"
}

# The datasets counted from raw reads, in the order of the table of runs, and those whose reads this build counts.
mapfile -t raw_read_datasets < <(tail -n +2 "$data/counts/raw_read_runs.tsv" | cut -f1 | uniq)
to_count=()
for dataset in "${raw_read_datasets[@]}"; do
    [ -n "$counted_reads" ] && [ -e "$counted_reads/$dataset.tsv.gz" ] && continue
    to_count+=("$dataset")
done

make_environment count-tables
py=$envs/count-tables/bin/python
"$py" "$pipeline/downloads.py" --downloads "$data/downloads.tsv" --data-dir "$downloads" check \
    --count-datasets "${to_count[@]}"
if [ ${#to_count[@]} -gt 0 ]; then
    make_environment raw-read-counts
    make_environment hisat2-2.1.0
fi
log "threads $threads; ${#raw_read_datasets[@]} datasets counted from raw reads, ${#to_count[@]} of them by this build"

prepare_reference() {
    local dir=$1 reference=$downloads/reference
    mkdir -p "$dir/hisat2_index"
    tar -xzf "$reference/analysis_set/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.hisat2_index.tar.gz" \
        -C "$dir/hisat2_index"
    "$py" "$pipeline/make_geo_gene_features.py" \
        --gtf "$reference/annotation_release/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz" \
        --assembly-report "$reference/annotation_release/GCF_000001405.39_GRCh38.p13_assembly_report.txt" \
        --analysis-set-index "$reference/analysis_set/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai" \
        --out "$dir/geo_gene_features_109.20190905_no_alt.saf"
}

# Each sample is aligned with the number of threads of our count (the column `threads` of the table of runs), because
# HISAT2's choices depend on it; as many samples as --threads allows are counted at a time.
count_dataset() {
    local dataset=$1 dir=$2 indexes sample_threads
    mapfile -t indexes < <(find "$out/reference/hisat2_index" -name '*.1.ht2')
    [ ${#indexes[@]} -eq 1 ] || { echo "${#indexes[@]} HISAT2 indexes in $out/reference/hisat2_index" >&2; exit 1; }
    sample_threads=$(awk -F '\t' -v dataset="$dataset" 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "threads") column = i }
        $1 == dataset { print $column }' "$data/counts/raw_read_runs.tsv" | sort -n | tail -1)
    mkdir -p "$dir"
    PATH="$envs/raw-read-counts/bin:$PATH" "$envs/raw-read-counts/bin/python" "$pipeline/count_raw_reads.py" \
        --runs "$data/counts/raw_read_runs.tsv" --datasets "$dataset" --reads-dir "$downloads/raw_reads" \
        --hisat2 "$envs/hisat2-2.1.0/bin/hisat2" --index "${indexes[0]%.1.ht2}" \
        --features "$out/reference/geo_gene_features_109.20190905_no_alt.saf" --annotation "$annotation" \
        --option-sets "$pipeline/raw_read_counts_options.json" --use "$option_set" --out-dir "$dir/counts" \
        --units-dir "$dir" --work-dir "$dir/alignments" --workers $(( threads > sample_threads ? threads / sample_threads : 1 ))
}

build_sample_table() {
    local dir=$1
    mkdir -p "$dir/harvest_geo" "$dir/harvest_sheets"
    "$py" "$pipeline/harvest_geo_samples.py" --geo-root "$downloads/geo" --out-dir "$dir" --units-dir "$dir/harvest_geo" \
        --workers "$threads"
    "$py" "$pipeline/harvest_sample_sheets.py" --data-root "$downloads" --out "$dir/sample_sheets_downloads.tsv.gz" \
        --units-dir "$dir/harvest_sheets"
    "$py" "$pipeline/build_sample_table.py" --samples-dir "$dir" --inventory "$data/sample_table/inventory.tsv" \
        --out "$dir/samples.tsv.gz" --coverage "$dir/sample_table_coverage.tsv" \
        --states "$dir/sample_table_senescence_states.tsv"
    "$py" "$pipeline/prepare_label_inputs.py" --datasets "$data/labels/dataset_labels.tsv" \
        --samples "$dir/samples.tsv.gz" --geo-series "$dir/geo_series.tsv" --geo-samples "$dir/geo_samples.tsv.gz" \
        --geo-characteristics "$dir/geo_sample_characteristics.tsv.gz" \
        --sample-sheets "$dir/sample_sheets_downloads.tsv.gz" --arrayexpress-dir "$downloads/biostudies_arrayexpress" \
        --out "$dir/label_inputs.jsonl.gz"
}

# The count tables of GEO series, then of ArrayExpress experiments, whose processed files are in another folder; the
# records of the tables that the round never reads are copied, so that the table of count tables lists every dataset.
build_count_tables() {
    local dir=$1 repository
    mkdir -p "$dir/records"
    cp -p "$data/counts/records_not_rebuilt/"*.json "$dir/records/"
    for repository in geo biostudies_arrayexpress; do
        mkdir -p "$dir/units_$repository"
        "$py" "$pipeline/build_count_tables.py" --datasets "$data/counts/rebuilt_datasets_$repository.tsv" \
            --uniform-state "$data/counts/geo_uniform_counts_state.tsv" \
            --label-inputs "$out/sample_table/label_inputs.jsonl.gz" --geo-dir "$downloads/$repository" \
            --uniform-dir "$downloads/geo_uniform_counts" --annotation "$annotation" \
            --layouts-dir "$data/counts/layouts" --column-reassignments "$data/counts/count_column_reassignments.tsv" \
            --counted-dir "$out/raw_read_counts/counted" --out-dir "$dir" --units-dir "$dir/units_$repository" \
            --workers "$threads"
    done
}

build_round() {
    local dir=$1
    "$py" "$pipeline/build_benchmark_bulk.py" --sample-labels "$data/labels/sample_labels.tsv.gz" \
        --dataset-labels "$data/labels/dataset_labels.tsv" --samples "$out/sample_table/samples.tsv.gz" \
        --count-tables "$out/counts/count_tables.tsv" --count-columns "$out/counts/count_table_columns.tsv.gz" \
        --counts-dir "$out/counts" --annotation "$annotation" --min-genes 10000 --max-median-genes-with-reads 35000 \
        --out-dir "$dir"
    "$py" "$pipeline/benchmark_bulk_classes.py" --samples "$dir/benchmark_bulk_samples.tsv" --out-dir "$dir"
    "$py" "$pipeline/build_study_split.py" --samples "$dir/benchmark_bulk_samples.tsv" \
        --contacts "$data/split/geo_series_contacts.tsv" --lists "$data/split/method_source_datasets.tsv" \
        --main-scores "$data/split/main_scores.tsv" --accession-pubmed-ids "$data/split/accession_pubmed_ids.tsv" \
        --laboratory-group-links "$data/split/laboratory_group_links.tsv" \
        --split-amendments "$data/split/split_amendments.tsv" --out-dir "$dir/split"
    "$py" "$pipeline/write_final_sets.py" --split-dir "$dir/split" --samples "$dir/benchmark_bulk_samples.tsv" \
        --label-corrections "$data/labels/label_corrections.tsv"
}

build_agentomics_dataset() {
    local dir=$1
    mkdir -p "$dir"
    "$py" "$pipeline/build_agentomics_dataset.py" --counts "$out/$round/benchmark_bulk_counts.tsv.gz" \
        --training "$out/$round/split/training_data_samples.tsv" --test "$out/$round/split/test_samples.tsv" \
        --annotation "$annotation" --validation-share 0.5 --dataset-dir "$dir/senescence_bulk_$round" \
        --key-dir "$dir/senescence_bulk_${round}_key"
}

if [ ${#to_count[@]} -gt 0 ]; then
    run_step reference prepare_reference
fi
mkdir -p "$out/raw_read_counts/counted"
for dataset in "${raw_read_datasets[@]}"; do
    counted=$out/raw_read_counts/counted/$dataset.tsv.gz
    if [ -n "$counted_reads" ] && [ -e "$counted_reads/$dataset.tsv.gz" ]; then
        ln -sfn "$counted_reads/$dataset.tsv.gz" "$counted"
        continue
    fi
    run_step "raw_read_counts/$dataset" count_dataset "$dataset"
    ln -sfn "../$dataset/counts/tables/$option_set/$dataset.tsv.gz" "$counted"
done
run_step sample_table build_sample_table
run_step counts build_count_tables
run_step "$round" build_round
run_step agentomics build_agentomics_dataset

log "comparison with our build"
not_built=()
[ -e "$out/reference" ] || not_built+=(reference)
"$py" "$pipeline/compare_checksums.py" --checksums "$repo/checksums/$round.tsv" --build-dir "$out" \
    --out "$out/checksum_comparison.tsv" --not-built "${not_built[@]}"
