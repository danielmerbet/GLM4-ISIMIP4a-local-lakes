#!/bin/bash
# Fetch the pinned GLM-AED source into upstream/ (needed by build_glm4.sh).
#   bash get_upstream.sh
# glm-aed commit 83e20ed (2026-09-08) and the submodule commits it pins:
#   GLM ade47de, libaed-api d62ed0d, libaed-benthic c6d4b0e, libaed-demo d35b05e,
#   libaed-water 86b892e, libplot 2a540c7, libutil c8fe70e
# build_glm4.sh then applies patches/*.patch to the submodules and builds bin/glm-built.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
COMMIT=83e20ed4312fe530f73a733fac5246e4e174bb7c
git clone --filter=blob:none --no-checkout https://github.com/AquaticEcoDynamics/glm-aed.git "$HERE/upstream"
cd "$HERE/upstream"
git sparse-checkout set glm-source
git checkout "$COMMIT"
git submodule update --init --recursive glm-source
git submodule status --recursive
