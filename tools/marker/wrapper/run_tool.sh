#!/usr/bin/env bash
# Score the samples of one study with the senescence gene sets that the markeR package ships.
# Usage: bash run_tool.sh --counts <counts.tsv.gz> --annotation <annotation.tsv.gz> --out-dir <new folder>
#            [--env-prefix <environment made by setup_env.sh>]
# --out-dir must not exist.
set -euo pipefail

env_prefix=/SCRATCH/dtzim01/senescence/tools/envs/marker
counts=
annotation=
out_dir=
while [ $# -gt 0 ]; do
    case $1 in
        --counts) counts=$2 ;;
        --annotation) annotation=$2 ;;
        --out-dir) out_dir=$2 ;;
        --env-prefix) env_prefix=$2 ;;
        *) echo "error: unknown argument $1" >&2; exit 2 ;;
    esac
    shift 2
done
[ -n "$counts" ] && [ -n "$annotation" ] && [ -n "$out_dir" ] || {
    echo "error: --counts, --annotation and --out-dir are required" >&2; exit 2; }
[ ! -e "$out_dir" ] || { echo "error: $out_dir exists" >&2; exit 1; }

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# R_LIBS_USER points to a folder that does not exist, so that only the packages of the environment are loaded.
export R_LIBS_USER=$env_prefix/no_user_library
# One thread: the scores do not depend on the linear algebra library, and the node is shared.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
"$env_prefix/bin/Rscript" --vanilla "$script_dir/score_gene_sets.R" \
    "$counts" "$annotation" "$out_dir" "$env_prefix/markeR_source.tsv" "$script_dir/prepare_marker_input.R"
