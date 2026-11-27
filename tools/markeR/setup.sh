#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$here/../common.sh"
prefix=$envs/markeR
commit=c6528f0fa1570602b077e4688865dcd68af3d058

clone_tool https://github.com/DiseaseTranscriptomicsLab/markeR "$commit" "$clones/markeR"
micromamba create --yes --prefix "$prefix" --file "$here/conda_explicit_lock.txt"
mkdir "$prefix/markeR_build"
git -C "$clones/markeR" archive "$commit" | tar -x -C "$prefix/markeR_build"
export R_LIBS_USER=$prefix/no_user_library
"$prefix/bin/R" --vanilla CMD INSTALL --library="$prefix/lib/R/library" "$prefix/markeR_build"
rm -r "$prefix/markeR_build"
"$prefix/bin/Rscript" --vanilla -e \
    'p <- installed.packages()[, c("Package", "Version")]; write.table(p[order(p[, "Package"]), ], sep = "\t", quote = FALSE, row.names = FALSE)' \
    > "$prefix/r_package_versions.tsv"
diff "$here/r_package_versions.tsv" "$prefix/r_package_versions.tsv"
