#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$here/../common.sh"
export R_LIBS_USER=$envs/markeR/no_user_library
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$envs/markeR/bin/Rscript" --vanilla "$here/score_markeR.R" "$1" "$2" "$3"
