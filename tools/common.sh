source "$(dirname "${BASH_SOURCE[0]}")/../pipeline/env.sh"
clones=$repo/build/tools

clone_tool() {
    local url=$1 commit=$2 clone=$3
    if [ ! -e "$clone" ]; then
        git init --quiet "$clone"
        git -C "$clone" remote add origin "$url"
        git -C "$clone" fetch --quiet --depth 1 origin "$commit"
        git -C "$clone" checkout --quiet --detach FETCH_HEAD
    fi
    [ "$(git -C "$clone" rev-parse HEAD)" = "$commit" ] || { echo "error: $clone is not at commit $commit" >&2; exit 1; }
    [ -z "$(git -C "$clone" status --porcelain)" ] || { echo "error: $clone differs from commit $commit" >&2; exit 1; }
}
