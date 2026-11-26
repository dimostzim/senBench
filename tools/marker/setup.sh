#!/usr/bin/env bash
# Clones markeR at the commit the benchmark used and builds the environment of its wrapper, in an install folder.
# Usage: bash setup.sh <install folder>
# The conda packages come from the wrapper's lock. markeR is installed from its commit as the wrapper's setup_env.sh
# installs it (through git archive, so the clone is not written to), and the versions of the R packages must equal
# the wrapper's r_package_versions.tsv.
set -euo pipefail

install_dir=$(realpath --canonicalize-missing "$1")
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$here/../common.sh"
clone=$install_dir/markeR
commit=c6528f0fa1570602b077e4688865dcd68af3d058
env_prefix=$install_dir/envs/marker
build_dir=$env_prefix/markeR_build

prepare_install_dir "$install_dir"
source "$install_dir/envs/env.sh"
clone_tool https://github.com/DiseaseTranscriptomicsLab/markeR "$commit" "$clone"
create_env_from_lock "$env_prefix" "$here/wrapper/conda_explicit_lock.txt"

mkdir "$build_dir"
git -C "$clone" archive "$commit" | tar -x -C "$build_dir"
# R_LIBS_USER points to a folder that does not exist, so that no package of a home directory is read or written.
export R_LIBS_USER=$env_prefix/no_user_library
"$env_prefix/bin/R" --vanilla CMD INSTALL --library="$env_prefix/lib/R/library" "$build_dir"
rm -r "$build_dir"
# run_tool.sh writes this file's url and commit into run_info.json.
printf 'url\tcommit\n%s\t%s\n' "$(git -C "$clone" remote get-url origin)" "$commit" > "$env_prefix/markeR_source.tsv"

"$env_prefix/bin/Rscript" --vanilla -e \
    'p <- installed.packages()[, c("Package", "Version")]; write.table(p[order(p[, "Package"]), ], sep = "\t", quote = FALSE, row.names = FALSE)' \
    > "$env_prefix/r_package_versions.tsv"
diff "$here/wrapper/r_package_versions.tsv" "$env_prefix/r_package_versions.tsv"
check_clone "$clone" "$commit"
echo "markeR is ready: clone $clone, environment $env_prefix"
