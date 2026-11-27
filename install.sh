#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/pipeline/env.sh"
names=(pipeline Agentomics hUSI SenCID SenCat SENCAN markeR SenFinder)
[ $# -eq 0 ] || names=("$@")

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

install_micromamba
for name in "${names[@]}"; do
    case $name in
        pipeline) installer=(micromamba env create --yes --prefix "$envs/pipeline" --file "$repo/pipeline/environment.yml") ;;
        Agentomics) installer=(micromamba env create --yes --prefix "$envs/Agentomics" --file "$repo/agentomics_model/environment.yml") ;;
        *)
            [ -f "$repo/tools/$name/setup.sh" ] || { echo "error: no method $name" >&2; exit 1; }
            installer=(bash "$repo/tools/$name/setup.sh")
            ;;
    esac
    [ ! -e "$envs/$name/.complete" ] || continue
    echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) $name"
    rm -rf "${envs:?}/$name"
    "${installer[@]}"
    touch "$envs/$name/.complete"
done
