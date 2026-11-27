#!/usr/bin/env bash
# Applies the released SenCat RNA marker weights to the raw counts of one study.
# Usage: bash run_tool.sh --counts <counts.tsv.gz> --annotation <annotation.tsv.gz> --out-dir <new folder>
set -euo pipefail

# The folder that holds the unchanged tool clones and envs/ (with env.sh and the sencat environment).
TOOLS_DIR="${SENESCENCE_TOOLS_DIR:-/SCRATCH/dtzim01/senescence/tools}"

source "$TOOLS_DIR/envs/env.sh"
# No package from the home directory, and no .pyc written into the read-only clone.
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
export NUMBA_CACHE_DIR="$XDG_CACHE_HOME/numba-sencat"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMBA_NUM_THREADS=8

exec "$TOOLS_DIR/envs/sencat/bin/python" "$(dirname "$(readlink -f "$0")")/run_sencat.py" \
    --repository "$TOOLS_DIR/wf-ml-markers-senescence" "$@"
