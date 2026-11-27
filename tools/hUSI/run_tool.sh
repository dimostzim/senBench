#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$here/../common.sh"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$envs/hUSI/bin/Rscript" --vanilla "$here/score_husi.R" "$1" "$2" "$clones/HUSI" "$3"
