#!/usr/bin/env bash
# Creates the environment of the SenFinder wrapper and writes its lock file next to this script.
# Usage: source <envs folder>/env.sh; bash setup_env.sh <environment folder that does not exist>
set -euo pipefail

env_prefix=$1
[ ! -e "$env_prefix" ] || { echo "error: $env_prefix exists" >&2; exit 1; }
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

micromamba create --yes --prefix "$env_prefix" --override-channels --strict-channel-priority \
    --file "$script_dir/environment.yml"
micromamba env export --prefix "$env_prefix" --explicit --md5 > "$script_dir/conda_explicit_lock.txt"
