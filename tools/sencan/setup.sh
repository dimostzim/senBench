#!/usr/bin/env bash
# Clones the hUSI repository, which holds the SENCAN classifier, at the commit the benchmark used and builds the
# environment of the SENCAN wrapper, in an install folder.
# Usage: bash setup.sh [install folder]   (default: build/tools/ in the repository's folder)
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
install_dir=$(realpath --canonicalize-missing "${1:-$here/../../build/tools}")
source "$here/../common.sh"
clone=$install_dir/HUSI
commit=300da32f111f8c8603e822eca8e84c9b5a76db48
env_prefix=$install_dir/envs/sencan

prepare_install_dir "$install_dir"
source "$install_dir/envs/env.sh"
clone_tool https://github.com/WJPina/HUSI "$commit" "$clone"
create_env_from_lock "$env_prefix" "$here/conda_explicit_lock.txt"
check_clone "$clone" "$commit"
echo "SENCAN is ready: clone $clone, environment $env_prefix"
