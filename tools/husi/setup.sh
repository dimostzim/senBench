#!/usr/bin/env bash
# Clones hUSI at the commit the benchmark used and builds the environment of its wrapper, in an install folder.
# Usage: bash setup.sh <install folder>
# The environment is built from the wrapper's conda lock; the R packages are built from the dated CRAN snapshot that
# the wrapper names, and their versions must equal the wrapper's r_packages_lock.tsv.
set -euo pipefail

install_dir=$(realpath --canonicalize-missing "$1")
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$here/../common.sh"
clone=$install_dir/HUSI
commit=300da32f111f8c8603e822eca8e84c9b5a76db48
env_prefix=$install_dir/envs/husi

prepare_install_dir "$install_dir"
source "$install_dir/envs/env.sh"
clone_tool https://github.com/WJPina/HUSI "$commit" "$clone"
create_env_from_lock "$env_prefix" "$here/wrapper/conda_explicit_lock.txt"
# micromamba run activates the environment, which puts its compilers on the path for the source installs.
micromamba run --prefix "$env_prefix" Rscript --vanilla "$here/wrapper/install_r_packages.R" \
    "$env_prefix/r_packages_lock.tsv"
diff "$here/wrapper/r_packages_lock.tsv" "$env_prefix/r_packages_lock.tsv"
check_clone "$clone" "$commit"
echo "hUSI is ready: clone $clone, environment $env_prefix"
