# Functions shared by the setup.sh of every tool: the install folder with micromamba, the clone of a tool at its
# commit, and an environment built from the exact package builds of a wrapper's lock file.
# shellcheck shell=bash

# micromamba 2.3.0, the version that built the benchmark's environments; the sha256 is that of its binary.
MICROMAMBA_URL=https://micro.mamba.pm/api/micromamba/linux-64/2.3.0
MICROMAMBA_SHA256=83f38b81295dade6b242b40dbab81bcd9fb03bb196928ab15f552cbafe82386b
TOOLS_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# Creates <install folder>/envs with env.sh and micromamba, and checks the micromamba binary.
prepare_install_dir() {
    local envs=$1/envs
    mkdir -p "$envs/bin" "$envs/cache" "$envs/tmp"
    cp "$TOOLS_DIR/env.sh" "$envs/env.sh"
    if [ ! -e "$envs/bin/micromamba" ]; then
        local archive
        archive=$(mktemp "$envs/tmp/micromamba.XXXXXX")
        curl --fail --silent --show-error --location --output "$archive" "$MICROMAMBA_URL"
        tar --extract --bzip2 --to-stdout --file "$archive" bin/micromamba > "$archive.bin"
        rm "$archive"
        chmod +x "$archive.bin"
        mv "$archive.bin" "$envs/bin/micromamba"
    fi
    echo "$MICROMAMBA_SHA256  $envs/bin/micromamba" | sha256sum --check --quiet
}

# Stops unless <clone> is at <commit> and holds no changed or new file.
check_clone() {
    local clone=$1 commit=$2
    [ "$(git -C "$clone" rev-parse HEAD)" = "$commit" ] || {
        echo "error: $clone is not at commit $commit" >&2; exit 1; }
    [ -z "$(git -C "$clone" status --porcelain)" ] || {
        echo "error: $clone differs from commit $commit" >&2; exit 1; }
}

# Clones <url> at <commit> into <clone> (only that commit, without the history), or checks an existing clone.
clone_tool() {
    local url=$1 commit=$2 clone=$3
    if [ ! -e "$clone" ]; then
        git init --quiet "$clone"
        git -C "$clone" remote add origin "$url"
        git -C "$clone" fetch --quiet --depth 1 origin "$commit"
        git -C "$clone" checkout --quiet --detach FETCH_HEAD
    fi
    check_clone "$clone" "$commit"
}

# Creates the environment <prefix> from an explicit lock file, which names the exact build of every conda package.
create_env_from_lock() {
    local prefix=$1 lock=$2
    [ ! -e "$prefix" ] || { echo "error: $prefix exists; delete it to build it again" >&2; exit 1; }
    micromamba create --yes --prefix "$prefix" --file "$lock"
}
