#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$here/../common.sh"

clone_tool https://github.com/WJPina/HUSI 300da32f111f8c8603e822eca8e84c9b5a76db48 "$clones/HUSI"
micromamba create --yes --prefix "$envs/SENCAN" --file "$here/conda_explicit_lock.txt"
