#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$here/../common.sh"
prefix=$envs/hUSI

clone_tool https://github.com/WJPina/HUSI 300da32f111f8c8603e822eca8e84c9b5a76db48 "$clones/HUSI"
micromamba create --yes --prefix "$prefix" --file "$here/conda_explicit_lock.txt"
micromamba run --prefix "$prefix" Rscript --vanilla "$here/install_r_packages.R" "$prefix/r_packages_lock.tsv"
diff "$here/r_packages_lock.tsv" "$prefix/r_packages_lock.tsv"
