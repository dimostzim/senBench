#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$here/../common.sh"
prefix=$envs/SenCID
clone=$clones/SenCID
commit=f93cba19b9d44e84f89504ecc6aaeb108ce8e508
export PYTHONNOUSERSITE=1

clone_tool https://github.com/JackieHanLab/SenCID "$commit" "$clone"
micromamba create --yes --prefix "$prefix" --file "$here/environment.explicit.lock"
cd "$TMPDIR"
micromamba run --prefix "$prefix" python -m pip install --no-deps -r "$here/environment.pip.lock"
build_dir=$(mktemp -d "$TMPDIR/sencid-build.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT
git -C "$clone" archive --prefix=SenCID/ "$commit" | tar -x -C "$build_dir"
micromamba run --prefix "$prefix" python -m pip install --no-deps --no-build-isolation "$build_dir/SenCID"
micromamba run --prefix "$prefix" python -m pip check
installed_dir=$("$prefix/bin/python" -I -c 'import os, SenCID; print(os.path.dirname(SenCID.__file__))')
diff -r --exclude=__pycache__ --exclude=demo "$clone/SenCID" "$installed_dir"
