#!/usr/bin/env bash
# Clones SenCID at the commit the benchmark used, builds the environment of its wrapper and installs SenCID into it
# from the clone, in an install folder.
# Usage: bash setup.sh [install folder]   (default: build/tools/ in the repository's folder)
# The wrapper reads its paths from config.sh, which names the folders of the benchmark's host. So the wrapper runs
# from a copy of its three scripts in <install folder>/wrappers/sencid/, whose config.sh names this install folder;
# the other two equal the wrapper's. The conda and pip packages come from the wrapper's two locks; SenCID is installed
# from the clone as the benchmark installed it.
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
install_dir=$(realpath --canonicalize-missing "${1:-$here/../../build/tools}")
source "$here/../common.sh"
clone=$install_dir/SenCID
commit=f93cba19b9d44e84f89504ecc6aaeb108ce8e508
env_prefix=$install_dir/envs/sencid
wrapper_copy=$install_dir/wrappers/sencid

[ ! -e "$wrapper_copy" ] || { echo "error: $wrapper_copy exists; delete it to set up again" >&2; exit 1; }
prepare_install_dir "$install_dir"
source "$install_dir/envs/env.sh"
export PYTHONNOUSERSITE=1
clone_tool https://github.com/JackieHanLab/SenCID "$commit" "$clone"
create_env_from_lock "$env_prefix" "$here/environment.explicit.lock"
# pip writes a temporary file into the current folder.
cd "$TMPDIR"
micromamba run --prefix "$env_prefix" python -m pip install --no-deps -r "$here/environment.pip.lock"

# pip builds inside the source folder, so SenCID is installed from a copy of the commit.
build_dir=$(mktemp -d "$TMPDIR/sencid-build.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT
git -C "$clone" archive --prefix=SenCID/ "$commit" | tar -x -C "$build_dir"
micromamba run --prefix "$env_prefix" python -m pip install --no-deps --no-build-isolation "$build_dir/SenCID"
micromamba run --prefix "$env_prefix" python -m pip check
# The installed code and model files must be the files of the clone.
installed_dir=$("$env_prefix/bin/python" -I -c 'import os, SenCID; print(os.path.dirname(SenCID.__file__))')
diff -r --exclude=__pycache__ --exclude=demo "$clone/SenCID" "$installed_dir"
printf '{"repository": "https://github.com/JackieHanLab/SenCID", "commit": "%s"}\n' "$commit" \
    > "$env_prefix/sencid_source.json"

mkdir -p "$wrapper_copy"
cp "$here/config.sh" "$here/run_tool.sh" "$here/run_sencid.py" "$wrapper_copy/"
sed -i -e "s|^ENVS_SETUP=.*|ENVS_SETUP=$install_dir/envs/env.sh|" -e "s|^ENV_DIR=.*|ENV_DIR=$env_prefix|" \
    -e "s|^TOOL_CLONE=.*|TOOL_CLONE=$clone|" "$wrapper_copy/config.sh"
grep --quiet "^TOOL_COMMIT=$commit$" "$wrapper_copy/config.sh"
check_clone "$clone" "$commit"
echo "SenCID is ready: clone $clone, environment $env_prefix, wrapper $wrapper_copy/run_tool.sh"
