# Paths and machine limits of the SenCID wrapper. The only place where they are set.
# shellcheck shell=bash

# Sets micromamba, the package and cache folders, TMPDIR and the proxy.
ENVS_SETUP=/SCRATCH/dtzim01/senescence/tools/envs/env.sh
# The environment of the wrapper.
ENV_DIR=/SCRATCH/dtzim01/senescence/tools/envs/sencid
# The unchanged clone of https://github.com/JackieHanLab/SenCID. Read only.
TOOL_CLONE=/SCRATCH/dtzim01/senescence/tools/SenCID
# The commit that install_env.sh expects in TOOL_CLONE.
TOOL_COMMIT=f93cba19b9d44e84f89504ecc6aaeb108ce8e508
# Upper limit of threads for the numeric libraries.
MAX_THREADS=8
