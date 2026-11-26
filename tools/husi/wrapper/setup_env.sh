#!/usr/bin/env bash
# Creates the environment of the hUSI wrapper and writes its two lock files next to this script.
# Usage: source <envs folder>/env.sh; bash setup_env.sh <environment folder that does not exist>
set -euo pipefail

env_prefix=$1
[ ! -e "$env_prefix" ] || { echo "error: $env_prefix exists" >&2; exit 1; }
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

micromamba create --yes --prefix "$env_prefix" --file "$script_dir/environment.yml"
# micromamba run activates the environment, which puts its compilers on the path for the source installs.
micromamba run --prefix "$env_prefix" Rscript --vanilla "$script_dir/install_r_packages.R" \
    "$script_dir/r_packages_lock.tsv"
micromamba env export --prefix "$env_prefix" --explicit --md5 > "$script_dir/conda_explicit_lock.txt"
