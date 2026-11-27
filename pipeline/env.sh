repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
envs=$repo/build/envs
export MAMBA_ROOT_PREFIX=$envs/micromamba
export CONDA_PKGS_DIRS=$envs/micromamba/pkgs
export PIP_CACHE_DIR=$envs/cache/pip
export XDG_CACHE_HOME=$envs/cache/xdg
export TMPDIR=$envs/tmp
export PATH=$envs/bin:$PATH

require() {
    [ -e "$envs/$1/.complete" ] || { echo "error: $1 is not installed; run: bash install.sh $1" >&2; exit 1; }
}
