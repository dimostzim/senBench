#!/usr/bin/env bash
# Score the samples of one study with hUSI (human universal senescence index).
# Usage: bash run_tool.sh --counts <counts.tsv.gz> --annotation <annotation.tsv.gz> --out-dir <new folder>
#            [--env-prefix <environment made by setup_env.sh>] [--tool-repository <clone of WJPina/HUSI>]
# --out-dir must not exist. The folder that holds the environment must also hold env.sh.
set -euo pipefail

HUSI_COMMIT=300da32f111f8c8603e822eca8e84c9b5a76db48

env_prefix=/SCRATCH/dtzim01/senescence/tools/envs/husi
tool_repository=/SCRATCH/dtzim01/senescence/tools/HUSI
counts=
annotation=
out_dir=
while [ $# -gt 0 ]; do
    case $1 in
        --counts) counts=$2 ;;
        --annotation) annotation=$2 ;;
        --out-dir) out_dir=$2 ;;
        --env-prefix) env_prefix=$2 ;;
        --tool-repository) tool_repository=$2 ;;
        *) echo "error: unknown argument $1" >&2; exit 2 ;;
    esac
    shift 2
done
[ -n "$counts" ] && [ -n "$annotation" ] && [ -n "$out_dir" ] || {
    echo "error: --counts, --annotation and --out-dir are required" >&2; exit 2; }
[ ! -e "$out_dir" ] || { echo "error: $out_dir exists" >&2; exit 1; }

commit=$(git -C "$tool_repository" rev-parse HEAD)
[ "$commit" = "$HUSI_COMMIT" ] || {
    echo "error: $tool_repository is at commit $commit, not $HUSI_COMMIT" >&2; exit 1; }
[ -z "$(git -C "$tool_repository" status --porcelain)" ] || {
    echo "error: $tool_repository differs from commit $HUSI_COMMIT" >&2; exit 1; }

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# env.sh keeps temporary files off the system disk.
source "$(dirname "$env_prefix")/env.sh"
# Scoring does no linear algebra; without this the linear algebra library starts one thread per CPU.
export OPENBLAS_NUM_THREADS=1
"$env_prefix/bin/Rscript" --vanilla "$script_dir/score_husi.R" \
    "$counts" "$annotation" "$tool_repository" "$commit" "$out_dir"
