#!/bin/bash
# Thin test orchestrator for ACES_PROGRAM. Runs each subsystem's OWN
# regression suite with the shared toolchain from ./configure, then prints
# one rollup table -- so "are we still at the zero-point baseline
# (zero-point-20260812: aces2 121/1/3, aces3 37/0, pyaces 37/16/4/1)" is a
# single command instead of three separately-remembered invocations. Does
# not reimplement any suite's own pass/fail logic, only parses each one's
# existing summary line.
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ ! -f "$ROOT/build_config.sh" ]; then
  echo "no build_config.sh -- run ./configure first" >&2
  exit 1
fi
# shellcheck source=/dev/null
source "$ROOT/build_config.sh"

# How the toolchain gets onto PATH (modules / a setvars script / nothing at
# all) is whatever ./configure probed on this host; build_config.sh carries the
# loader itself. Older build_config.sh files predate it and were modules-only.
if type aces_program_load_env >/dev/null 2>&1; then
  aces_program_load_env || exit 1
else
  module purge
  # shellcheck disable=SC2086
  if ! module load $ACES_PROGRAM_MODULES; then
    echo "module load failed: $ACES_PROGRAM_MODULES -- re-run ./configure" >&2
    exit 1
  fi
fi

mkdir -p "$ROOT/build_logs"

# The aces2 binaries the suites run (xaces2, xjoda, ...): wherever the gmake
# build installed them -- this checkout's aces2/bin, or the DIR_INSTBIN an
# untracked aces2/Makefiles/local.mk sets (e.g. $(HOME)/Develop/bin). ACES2_BIN
# overrides both.
if [ -z "${ACES2_BIN:-}" ]; then
  if [ -f "$ROOT/aces2/Makefiles/local.mk" ]; then
    ACES2_BIN=$(printf 'include local.mk\nall:;@echo $(DIR_INSTBIN)\n' |
                make -s --no-print-directory -C "$ROOT/aces2/Makefiles" -f - 2>/dev/null)
  fi
  ACES2_BIN=${ACES2_BIN:-$ROOT/aces2/bin}
fi

TARGETS="$*"
[ -z "$TARGETS" ] && TARGETS="aces2 ring aces3 pyaces"

for t in $TARGETS; do
  echo "=== testing $t ==="
  LOG="$ROOT/build_logs/unified_test_$t.log"
  case "$t" in
    aces2)
      PATH="$ACES2_BIN:$PATH" bash -c \
        "cd '$ROOT/aces2' && python3 tests/fast/run_tests.py --jobs $ACES_PROGRAM_NP" \
        > "$LOG" 2>&1
      ;;
    aces3)
      # xaces3 runs the ACES II xjoda, found on PATH.
      ( cd "$ROOT/aces3" && PATH="$ACES2_BIN:$PATH" \
          python3 tests/fast/run_tests.py --jobs "$ACES_PROGRAM_NP" --np "$ACES_PROGRAM_NP" ) \
        > "$LOG" 2>&1
      ;;
    pyaces)
      # The suite gives every case its own working directory. The runner keeps
      # them under pyaces/suite_run by default; on HiPerGator they go to the
      # /blue scratch filesystem as before -- unless the caller already named
      # somewhere better.
      if [ -z "${PYACES_SUITE_RUN_DIR:-}" ] && [ -d /blue/ufhpc ]; then
        export PYACES_SUITE_RUN_DIR=/blue/ufhpc/perera/aces_program_pyaces_suite_test
      fi
      ( cd "$ROOT/pyaces" && PATH="$ACES2_BIN:$PATH" python3 pyaces_suite_runner.py ) \
        > "$LOG" 2>&1
      ;;
    ring)
      # The ACES II ring-CC/RPA suites, each with its own runner: ringcc,
      # ringcl (RHF vs UHF DRCCD Lambda), rpa, jcp2020_eom_rccd (published
      # EOM ring-CC excitation energies). Seconds each.
      : > "$LOG"
      for s in ringcc ringcl rpa jcp2020_eom_rccd; do
        echo "=== $s" >> "$LOG"
        PATH="$ACES2_BIN:$PATH" bash -c \
          "cd '$ROOT/aces2' && python3 tests/$s/run_tests.py --jobs $ACES_PROGRAM_NP" \
          >> "$LOG" 2>&1
      done
      # The same ground-state and EOM cases through xaces2's ring-CC route.
      for s in ringcc jcp2020_eom_rccd; do
        echo "=== $s --xaces2" >> "$LOG"
        PATH="$ACES2_BIN:$PATH" bash -c \
          "cd '$ROOT/aces2' && python3 tests/$s/run_tests.py --xaces2 --jobs $ACES_PROGRAM_NP" \
          >> "$LOG" 2>&1
      done
      ;;
    *)
      echo "unknown target: $t (expected aces2, ring, aces3, pyaces)" >&2
      exit 1
      ;;
  esac
  echo "  done (log: build_logs/unified_test_$t.log)"
done

echo
echo "=== unified suite rollup ==="
printf "%-8s %8s %8s %8s %8s\n" "SUITE" "TOTAL" "PASS" "FAIL" "OTHER"
printf "%-8s %8s %8s %8s %8s\n" "-----" "-----" "----" "----" "-----"

for t in $TARGETS; do
  LOG="$ROOT/build_logs/unified_test_$t.log"
  TOTAL="?"; PASS="?"; FAIL="?"; OTHER="?"
  case "$t" in
    aces2)
      # -a: a crashed case can leave non-UTF-8 bytes in the log, and grep would
      # then answer only "binary file matches" and the rollup would show "?".
      LINE=$(grep -aoE '[0-9]+ case\(s\), [0-9]+ passed, [0-9]+ failed, [0-9]+ errored' "$LOG" | tail -1)
      if [ -n "$LINE" ]; then
        TOTAL=$(echo "$LINE" | grep -oE '^[0-9]+')
        PASS=$(echo "$LINE"  | grep -oE '[0-9]+ passed' | grep -oE '[0-9]+')
        FAIL=$(echo "$LINE"  | grep -oE '[0-9]+ failed' | grep -oE '[0-9]+')
        OTHER=$(echo "$LINE" | grep -oE '[0-9]+ errored' | grep -oE '[0-9]+')
      fi
      ;;
    ring)
      # One summary line per suite ("N case(s)" or "N pair(s), P passed, F
      # failed, E errored[, X known-crash ...]"); OTHER is errors + known crashes.
      LINES=$(grep -aoE '[0-9]+ (case|pair)\(s\), [0-9]+ passed, [0-9]+ failed' "$LOG")
      if [ -n "$LINES" ]; then
        TOTAL=0; PASS=0; FAIL=0
        while read -r N _ P _ F _; do
          TOTAL=$((TOTAL + N)); PASS=$((PASS + P)); FAIL=$((FAIL + F))
        done <<< "$LINES"
        OTHER=$((TOTAL - PASS - FAIL))
      fi
      ;;
    aces3)
      LINE=$(grep -aoE '[0-9]+ case\(s\), [0-9]+ failed,' "$LOG" | tail -1)
      if [ -n "$LINE" ]; then
        TOTAL=$(echo "$LINE" | grep -oE '^[0-9]+')
        FAIL=$(echo "$LINE"  | grep -oE '[0-9]+ failed' | grep -oE '[0-9]+')
        PASS=$((TOTAL - FAIL))
        OTHER=0
      fi
      ;;
    pyaces)
      # Read every "  STATUS: n" line of the runner's Summary block rather than
      # a fixed list: MISMATCH and TIMEOUT used to fall out of the totals, and
      # a run with no PASS at all printed "?" everywhere.
      SUMMARY=$(sed -n '/^=== Summary ===$/,/^$/p' "$LOG" | grep -aE '^  [A-Z_]+: [0-9]+$')
      if [ -n "$SUMMARY" ]; then
        PASS=0; FAIL=0; OTHER=0
        while read -r STATUS N; do
          case "${STATUS%:}" in
            PASS)                PASS=$((PASS + N)) ;;
            MISMATCH|FAIL_CLEAN) FAIL=$((FAIL + N)) ;;
            *)                   OTHER=$((OTHER + N)) ;;
          esac
        done <<< "$SUMMARY"
        TOTAL=$(( PASS + FAIL + OTHER ))
      fi
      ;;
  esac
  printf "%-8s %8s %8s %8s %8s\n" "$t" "$TOTAL" "$PASS" "$FAIL" "$OTHER"
done
