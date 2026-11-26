#!/usr/bin/env bash
# Clones wf-ml-markers-senescence (the SenCat RNA score) at the commit the benchmark used and builds the environment
# of its wrapper, in an install folder.
# Usage: bash setup.sh <install folder>
# The conda packages come from the wrapper's lock and the four pip packages from its README; the installed Python
# packages must equal the wrapper's pip_freeze.txt.
set -euo pipefail

install_dir=$(realpath --canonicalize-missing "$1")
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$here/../common.sh"
clone=$install_dir/wf-ml-markers-senescence
commit=0e23fcc926f31fc2206631ff9b64a852fb81fe89
env_prefix=$install_dir/envs/sencat

prepare_install_dir "$install_dir"
source "$install_dir/envs/env.sh"
export PYTHONNOUSERSITE=1
clone_tool https://github.com/maragkakislab/wf-ml-markers-senescence "$commit" "$clone"
create_env_from_lock "$env_prefix" "$here/wrapper/conda_explicit_lock.txt"
"$env_prefix/bin/python" -m pip install --no-user --no-deps \
    pydeseq2==0.5.4 formulaic-contrasts==1.0.0 session-info==1.0.1 stdlib-list==0.12.0
"$env_prefix/bin/python" -m pip list --format=freeze > "$env_prefix/pip_freeze.txt"
diff "$here/wrapper/pip_freeze.txt" "$env_prefix/pip_freeze.txt"
check_clone "$clone" "$commit"
echo "SenCat is ready: clone $clone, environment $env_prefix"
