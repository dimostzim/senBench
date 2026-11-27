#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$here/../common.sh"

clone_tool https://github.com/XiaoDongLab/SenFinder 5f56506300baab2ed4c20af5a52c6252c4578b50 "$clones/SenFinder"
micromamba create --yes --prefix "$envs/SenFinder" --file "$here/conda_explicit_lock.txt"
