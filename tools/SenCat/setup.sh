#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$here/../common.sh"
prefix=$envs/SenCat
export PYTHONNOUSERSITE=1

clone_tool https://github.com/maragkakislab/wf-ml-markers-senescence 0e23fcc926f31fc2206631ff9b64a852fb81fe89 \
    "$clones/wf-ml-markers-senescence"
micromamba create --yes --prefix "$prefix" --file "$here/conda_explicit_lock.txt"
"$prefix/bin/python" -m pip install --no-user --no-deps \
    pydeseq2==0.5.4 formulaic-contrasts==1.0.0 session-info==1.0.1 stdlib-list==0.12.0
"$prefix/bin/python" -m pip list --format=freeze > "$prefix/pip_freeze.txt"
diff "$here/pip_freeze.txt" "$prefix/pip_freeze.txt"
