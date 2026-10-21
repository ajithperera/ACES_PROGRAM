#!/bin/bash
# Runs the aces2 fast regression suite, driven by ctest. Does not
# reimplement the suite's own pass/fail logic -- only loads the shared
# toolchain, puts the freshly-built aces2/bin on PATH, and invokes the
# suite's own runner, same as phase-1's run_tests.sh did. Always exits 0 --
# whether the suite's numbers match the expected baseline is judged by
# cmake/rollup.sh from the log this writes, not by this script's exit code
# (the baseline intentionally includes known non-regression FAIL/ERROR
# cases, so "suite ran" and "suite matches baseline" are different
# questions).
set -uo pipefail

MODULES="$1"
ACES2_DIR="$2"
NP="$3"
LOG="$4"

module purge
# shellcheck disable=SC2086
module load $MODULES

mkdir -p "$(dirname "$LOG")"
PATH="$ACES2_DIR/bin:$PATH" bash -c "cd '$ACES2_DIR' && python3 tests/fast/run_tests.py --jobs $NP" > "$LOG" 2>&1
echo "aces2 fast suite exit code: $?" >> "$LOG"
exit 0
