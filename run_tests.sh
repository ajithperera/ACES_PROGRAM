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

module purge
# shellcheck disable=SC2086
module load $ACES_PROGRAM_MODULES

mkdir -p "$ROOT/build_logs"

TARGETS="$*"
[ -z "$TARGETS" ] && TARGETS="aces2 aces3 pyaces"

for t in $TARGETS; do
  echo "=== testing $t ==="
  LOG="$ROOT/build_logs/unified_test_$t.log"
  case "$t" in
    aces2)
      PATH="$ROOT/aces2/bin:$PATH" bash -c \
        "cd '$ROOT/aces2' && python3 tests/fast/run_tests.py --jobs $ACES_PROGRAM_NP" \
        > "$LOG" 2>&1
      ;;
    aces3)
      ( cd "$ROOT/aces3" && python3 tests/fast/run_tests.py --jobs "$ACES_PROGRAM_NP" --np "$ACES_PROGRAM_NP" ) \
        > "$LOG" 2>&1
      ;;
    pyaces)
      ( cd "$ROOT/pyaces" && python3 pyaces_suite_runner.py ) \
        > "$LOG" 2>&1
      ;;
    *)
      echo "unknown target: $t (expected aces2, aces3, pyaces)" >&2
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
      LINE=$(grep -oE '[0-9]+ case\(s\), [0-9]+ passed, [0-9]+ failed, [0-9]+ errored' "$LOG" | tail -1)
      if [ -n "$LINE" ]; then
        TOTAL=$(echo "$LINE" | grep -oE '^[0-9]+')
        PASS=$(echo "$LINE"  | grep -oE '[0-9]+ passed' | grep -oE '[0-9]+')
        FAIL=$(echo "$LINE"  | grep -oE '[0-9]+ failed' | grep -oE '[0-9]+')
        OTHER=$(echo "$LINE" | grep -oE '[0-9]+ errored' | grep -oE '[0-9]+')
      fi
      ;;
    aces3)
      LINE=$(grep -oE '[0-9]+ case\(s\), [0-9]+ failed,' "$LOG" | tail -1)
      if [ -n "$LINE" ]; then
        TOTAL=$(echo "$LINE" | grep -oE '^[0-9]+')
        FAIL=$(echo "$LINE"  | grep -oE '[0-9]+ failed' | grep -oE '[0-9]+')
        PASS=$((TOTAL - FAIL))
        OTHER=0
      fi
      ;;
    pyaces)
      P=$(grep -oE 'PASS: [0-9]+' "$LOG" | grep -oE '[0-9]+')
      FC=$(grep -oE 'FAIL_CLEAN: [0-9]+' "$LOG" | grep -oE '[0-9]+')
      CR=$(grep -oE 'CRASH: [0-9]+' "$LOG" | grep -oE '[0-9]+')
      NR=$(grep -oE 'NO_RESULT: [0-9]+' "$LOG" | grep -oE '[0-9]+')
      if [ -n "$P" ]; then
        PASS=$P; FAIL=${FC:-0}
        OTHER=$(( ${CR:-0} + ${NR:-0} ))
        TOTAL=$(( PASS + FAIL + OTHER ))
      fi
      ;;
  esac
  printf "%-8s %8s %8s %8s %8s\n" "$t" "$TOTAL" "$PASS" "$FAIL" "$OTHER"
done
