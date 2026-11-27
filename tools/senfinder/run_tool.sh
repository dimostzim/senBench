#!/usr/bin/env bash
# Classify the samples of one study with SenFinder, the released model of the Mahmud 2024 classifier.
# Usage: bash run_tool.sh --counts <counts.tsv.gz> --annotation <annotation.tsv.gz> --out-dir <new folder>
#            [--set-absent-genes-to-zero] [--env-prefix <environment made by setup_env.sh>]
#            [--tool-repository <clone of XiaoDongLab/SenFinder>]
# --out-dir must not exist. The folder that holds the environment must also hold env.sh.
# Without --set-absent-genes-to-zero the run stops when the count table has no gene for a gene of the model.
set -euo pipefail

SENFINDER_COMMIT=5f56506300baab2ed4c20af5a52c6252c4578b50

env_prefix=/SCRATCH/dtzim01/senescence/tools/envs/senfinder
tool_repository=/SCRATCH/dtzim01/senescence/tools/SenFinder
set_absent_genes_to_zero=FALSE
counts=
annotation=
out_dir=
while [ $# -gt 0 ]; do
    case $1 in
        --counts) counts=$2; shift ;;
        --annotation) annotation=$2; shift ;;
        --out-dir) out_dir=$2; shift ;;
        --env-prefix) env_prefix=$2; shift ;;
        --tool-repository) tool_repository=$2; shift ;;
        --set-absent-genes-to-zero) set_absent_genes_to_zero=TRUE ;;
        *) echo "error: unknown argument $1" >&2; exit 2 ;;
    esac
    shift
done
[ -n "$counts" ] && [ -n "$annotation" ] && [ -n "$out_dir" ] || {
    echo "error: --counts, --annotation and --out-dir are required" >&2; exit 2; }
[ ! -e "$out_dir" ] || { echo "error: $out_dir exists" >&2; exit 1; }

commit=$(git -C "$tool_repository" rev-parse HEAD)
[ "$commit" = "$SENFINDER_COMMIT" ] || {
    echo "error: $tool_repository is at commit $commit, not $SENFINDER_COMMIT" >&2; exit 1; }
[ -z "$(git -C "$tool_repository" status --porcelain)" ] || {
    echo "error: $tool_repository differs from commit $SENFINDER_COMMIT" >&2; exit 1; }

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# env.sh keeps temporary files off the system disk.
source "$(dirname "$env_prefix")/env.sh"
# The released script is run as it is, so the library outside the environment is switched off here.
unset R_LIBS
export R_LIBS_USER=NULL R_LIBS_SITE=NULL
# Without this the linear algebra library starts one thread per CPU.
export OPENBLAS_NUM_THREADS=1
"$env_prefix/bin/Rscript" --vanilla "$script_dir/score_senfinder.R" \
    "$counts" "$annotation" "$tool_repository" "$commit" "$script_dir/model_genes.tsv" \
    "$set_absent_genes_to_zero" "$out_dir"
