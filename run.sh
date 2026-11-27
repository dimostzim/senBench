#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/pipeline/env.sh"
python=$envs/pipeline/bin/python
methods=(Agentomics hUSI SenCID SenCat SENCAN markeR SenFinder)
[ $# -eq 0 ] || methods=("$@")

require pipeline
for method in "${methods[@]}"; do
    require "$method"
    "$python" "$repo/pipeline/predict.py" "$method"
done
"$python" "$repo/pipeline/score.py"
