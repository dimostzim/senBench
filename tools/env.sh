# Source this before installing or running a tool environment. It keeps micromamba's packages, the caches and the
# temporary files in the folder that holds this file, and puts its micromamba on the path. Proxy variables come from
# the calling shell. setup.sh copies this file to <install folder>/envs/env.sh.
# shellcheck shell=bash

senbench_envs_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export MAMBA_ROOT_PREFIX=$senbench_envs_dir/micromamba
export CONDA_PKGS_DIRS=$senbench_envs_dir/micromamba/pkgs
export PIP_CACHE_DIR=$senbench_envs_dir/cache/pip
export XDG_CACHE_HOME=$senbench_envs_dir/cache/xdg
export TMPDIR=$senbench_envs_dir/tmp
export PATH=$senbench_envs_dir/bin:$PATH
unset senbench_envs_dir
