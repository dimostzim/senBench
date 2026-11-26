#!/usr/bin/env bash
# Build the environment of the SenCID wrapper and install SenCID into it from the unchanged clone.
# Usage: bash install_env.sh
# Refuses to run when the environment folder exists. Writes the locks next to this script.
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/config.sh"
source "$ENVS_SETUP"
export PYTHONNOUSERSITE=1

if [[ -e "$ENV_DIR" ]]; then
    echo "error: $ENV_DIR exists" >&2
    exit 1
fi
clone_commit=$(git -C "$TOOL_CLONE" rev-parse HEAD)
if [[ "$clone_commit" != "$TOOL_COMMIT" ]]; then
    echo "error: $TOOL_CLONE is at $clone_commit, expected $TOOL_COMMIT" >&2
    exit 1
fi
if [[ -n "$(git -C "$TOOL_CLONE" status --porcelain)" ]]; then
    echo "error: $TOOL_CLONE has local changes" >&2
    exit 1
fi

# micromamba writes a temporary requirements file into the current folder.
cd "$TMPDIR"
micromamba create -y -p "$ENV_DIR" -f "$script_dir/environment.yml"

# pip builds inside the source folder, so SenCID is installed from a copy of the clone.
build_dir=$(mktemp -d "$TMPDIR/sencid-build.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT
git -C "$TOOL_CLONE" archive --prefix=SenCID/ "$TOOL_COMMIT" | tar -x -C "$build_dir"
micromamba run -p "$ENV_DIR" python -m pip install --no-deps --no-build-isolation "$build_dir/SenCID"
micromamba run -p "$ENV_DIR" python -m pip check

# The installed code and model files must be the files of the clone.
installed_dir=$("$ENV_DIR/bin/python" -I -c 'import os, SenCID; print(os.path.dirname(SenCID.__file__))')
diff -r --exclude=__pycache__ --exclude=demo "$TOOL_CLONE/SenCID" "$installed_dir"

printf '{"repository": "https://github.com/JackieHanLab/SenCID", "commit": "%s"}\n' "$TOOL_COMMIT" \
    > "$ENV_DIR/sencid_source.json"

micromamba env export --explicit -p "$ENV_DIR" > "$script_dir/environment.explicit.lock"
# SenCID itself is left out: pip would list the temporary folder; its commit is TOOL_COMMIT.
micromamba run -p "$ENV_DIR" python -m pip freeze --all --exclude SenCID > "$script_dir/environment.pip.lock"
echo "installed SenCID $TOOL_COMMIT into $ENV_DIR"
