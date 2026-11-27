#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$here/../common.sh"
unset R_LIBS
export R_LIBS_USER=NULL R_LIBS_SITE=NULL
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$envs/SenFinder/bin/Rscript" --vanilla "$here/score_senfinder.R" "$1" "$2" "$clones/SenFinder" \
    "$here/model_genes.tsv" "$3"
