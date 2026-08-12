#!/bin/bash
# Runs the aces3 fast regression suite, driven by ctest. See
# run_suite_aces2.sh for why this always exits 0 -- cmake/rollup.sh judges
# baseline parity from the log, not this script's exit code.
set -uo pipefail

MODULES="$1"
ACES3_DIR="$2"
NP="$3"
LOG="$4"

module purge
# shellcheck disable=SC2086
module load $MODULES

mkdir -p "$(dirname "$LOG")"
( cd "$ACES3_DIR" && python3 tests/fast/run_tests.py --jobs "$NP" --np "$NP" ) > "$LOG" 2>&1
echo "aces3 fast suite exit code: $?" >> "$LOG"
exit 0
