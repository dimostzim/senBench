#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/pipeline/env.sh"
python=$envs/pipeline/bin/python
methods=(Agentomics hUSI SenCID SenCat SENCAN markeR SenFinder)
[ $# -eq 0 ] || methods=("$@")

for method in "${methods[@]}"; do
    if [ "$method" = Agentomics ]; then
        install Agentomics micromamba env create --yes --prefix "$envs/Agentomics" \
            --file "$repo/agentomics_model/environment.yml"
    else
        install "$method" bash "$repo/tools/$method/setup.sh"
    fi
    "$python" "$repo/pipeline/predict.py" "$method"
done
"$python" "$repo/pipeline/score.py"
