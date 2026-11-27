#!/usr/bin/env bash
# Run SenCID on the raw counts of one study.
# Usage: bash run_tool.sh --counts <counts.tsv.gz> --annotation <annotation.tsv.gz> --out-dir <new folder>
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/config.sh"
source "$ENVS_SETUP"
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=$MAX_THREADS OPENBLAS_NUM_THREADS=$MAX_THREADS MKL_NUM_THREADS=$MAX_THREADS
export NUMEXPR_NUM_THREADS=$MAX_THREADS NUMBA_NUM_THREADS=$MAX_THREADS

exec micromamba run -p "$ENV_DIR" python "$script_dir/run_sencid.py" "$@"
