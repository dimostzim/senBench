#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$here/../common.sh"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$envs/SENCAN/bin/Rscript" --vanilla "$here/score_sencan.R" "$1" "$2" "$clones/HUSI" "$3"
