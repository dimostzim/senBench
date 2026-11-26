#!/usr/bin/env bash
# Create the environment of the markeR wrapper and install markeR into it from the unchanged clone.
# The clone is not written to: its checked-out commit is exported with `git archive` to a build folder
# inside the environment, installed from there, and the build folder is removed.
# Writes the two lock files next to this script.
# Usage: bash setup_env.sh <env_sh> <markeR_clone> <new_env_prefix>
set -euo pipefail

env_sh=$1
clone=$2
env_prefix=$3
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
build_dir=$env_prefix/markeR_build

[ ! -e "$env_prefix" ] || { echo "error: $env_prefix exists" >&2; exit 1; }

source "$env_sh"
micromamba create --yes --no-rc --strict-channel-priority --prefix "$env_prefix" --file "$script_dir/environment.yml"

commit=$(git -C "$clone" rev-parse HEAD)
mkdir "$build_dir"
git -C "$clone" archive "$commit" | tar -x -C "$build_dir"

# R_LIBS_USER points to a folder that does not exist, so that no package of a home directory is read or written.
export R_LIBS_USER=$env_prefix/no_user_library
"$env_prefix/bin/R" --vanilla CMD INSTALL --library="$env_prefix/lib/R/library" "$build_dir"
rm -r "$build_dir"

printf 'url\tcommit\n%s\t%s\n' "$(git -C "$clone" remote get-url origin)" "$commit" > "$env_prefix/markeR_source.tsv"

micromamba env export --explicit --prefix "$env_prefix" > "$script_dir/conda_explicit_lock.txt"
"$env_prefix/bin/Rscript" --vanilla -e \
    'p <- installed.packages()[, c("Package", "Version")]; write.table(p[order(p[, "Package"]), ], sep = "\t", quote = FALSE, row.names = FALSE)' \
    > "$script_dir/r_package_versions.tsv"
