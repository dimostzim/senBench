repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
envs=$repo/build/envs
export MAMBA_ROOT_PREFIX=$envs/micromamba
export CONDA_PKGS_DIRS=$envs/micromamba/pkgs
export PIP_CACHE_DIR=$envs/cache/pip
export XDG_CACHE_HOME=$envs/cache/xdg
export TMPDIR=$envs/tmp
export PATH=$envs/bin:$PATH

install_micromamba() {
    mkdir -p "$envs/bin" "$envs/cache" "$envs/tmp"
    if [ ! -e "$envs/bin/micromamba" ]; then
        curl --fail --silent --show-error --location https://micro.mamba.pm/api/micromamba/linux-64/2.3.0 |
            tar --extract --bzip2 --to-stdout bin/micromamba > "$envs/tmp/micromamba"
        chmod +x "$envs/tmp/micromamba"
        mv "$envs/tmp/micromamba" "$envs/bin/micromamba"
    fi
    echo "83f38b81295dade6b242b40dbab81bcd9fb03bb196928ab15f552cbafe82386b  $envs/bin/micromamba" |
        sha256sum --check --quiet
}

install() {
    local name=$1
    shift
    [ ! -e "$envs/$name/.complete" ] || return 0
    install_micromamba
    rm -rf "${envs:?}/$name"
    "$@"
    touch "$envs/$name/.complete"
}
