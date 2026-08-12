#!/bin/bash
# Builds pyaces via its own existing build.sh (f2py + a dozen-plus hard-won
# objcopy symbol-collision fixes) -- deliberately NOT ported into CMake
# line-by-line (see NOTES.md for why). This only drives it, exactly like
# phase-1's build.sh orchestrator did. pyaces/build.sh loads its own
# modules and derives its aces2 dependency path relative to its own
# location ($SCRIPT_DIR/../aces2), so it needs no arguments here beyond its
# own directory -- it just needs aces2 already built (see the CMake
# `DEPENDS aces2` on this step in the top-level CMakeLists.txt).
set -euo pipefail

PYACES_DIR="$1"
bash "$PYACES_DIR/build.sh"
