#!/usr/bin/env bash
# Clones SenFinder at the commit the benchmark used and builds the environment of its wrapper, in an install folder.
# Usage: bash setup.sh <install folder>
set -euo pipefail

install_dir=$(realpath --canonicalize-missing "$1")
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$here/../common.sh"
clone=$install_dir/SenFinder
commit=5f56506300baab2ed4c20af5a52c6252c4578b50
env_prefix=$install_dir/envs/senfinder

prepare_install_dir "$install_dir"
source "$install_dir/envs/env.sh"
clone_tool https://github.com/XiaoDongLab/SenFinder "$commit" "$clone"
create_env_from_lock "$env_prefix" "$here/wrapper/conda_explicit_lock.txt"
check_clone "$clone" "$commit"
echo "SenFinder is ready: clone $clone, environment $env_prefix"
