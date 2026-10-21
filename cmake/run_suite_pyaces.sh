#!/bin/bash
# Runs the pyaces regression suite, driven by ctest. See run_suite_aces2.sh
# for why this always exits 0 -- cmake/rollup.sh judges baseline parity
# from the log, not this script's exit code.
set -uo pipefail

MODULES="$1"
PYACES_DIR="$2"
LOG="$3"

module purge
# shellcheck disable=SC2086
module load $MODULES

mkdir -p "$(dirname "$LOG")"
( cd "$PYACES_DIR" && python3 pyaces_suite_runner.py ) > "$LOG" 2>&1
echo "pyaces suite exit code: $?" >> "$LOG"
exit 0
