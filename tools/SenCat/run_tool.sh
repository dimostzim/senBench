#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$here/../common.sh"
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
export NUMBA_CACHE_DIR=$XDG_CACHE_HOME/numba-sencat
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMBA_NUM_THREADS=8
exec "$envs/SenCat/bin/python" "$here/run_sencat.py" "$1" "$2" "$clones/wf-ml-markers-senescence" "$3"
