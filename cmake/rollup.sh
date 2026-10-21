#!/bin/bash
# Parses each subsystem's fast-suite log (written by cmake/run_suite_*.sh)
# into TOTAL/PASS/FAIL/OTHER, prints the unified rollup table (same shape
# phase-1's run_tests.sh printed), and gates on the zero-point-20260812
# baseline: aces2 121 pass/1 fail/3 error out of 125, aces3 37/0 out of 37,
# pyaces 37 pass/16 fail_clean/4 crash/1 no_result out of 58. Does not
# reimplement any suite's own pass/fail logic -- only parses each one's
# existing summary line, exactly like run_tests.sh did.
set -uo pipefail

ACES2_LOG="$1"
ACES3_LOG="$2"
PYACES_LOG="$3"

# zero-point-20260812 baseline, as TOTAL:PASS:FAIL:OTHER. Overridable via
# environment for a deliberate, reviewed re-baseline -- not meant to be
# edited casually to make a regression "pass".
BASELINE_ACES2="${ACES_BASELINE_ACES2:-125:121:1:3}"
BASELINE_ACES3="${ACES_BASELINE_ACES3:-37:37:0:0}"
BASELINE_PYACES="${ACES_BASELINE_PYACES:-58:37:16:5}"

parse_aces2() {
  local log="$1" line total pass fail other
  line=$(grep -oE '[0-9]+ case\(s\), [0-9]+ passed, [0-9]+ failed, [0-9]+ errored' "$log" 2>/dev/null | tail -1)
  [ -z "$line" ] && { echo "?:?:?:?"; return; }
  total=$(echo "$line" | grep -oE '^[0-9]+')
  pass=$(echo "$line"  | grep -oE '[0-9]+ passed' | grep -oE '[0-9]+')
  fail=$(echo "$line"  | grep -oE '[0-9]+ failed' | grep -oE '[0-9]+')
  other=$(echo "$line" | grep -oE '[0-9]+ errored' | grep -oE '[0-9]+')
  echo "$total:$pass:$fail:$other"
}

parse_aces3() {
  local log="$1" line total fail pass
  line=$(grep -oE '[0-9]+ case\(s\), [0-9]+ failed,' "$log" 2>/dev/null | tail -1)
  [ -z "$line" ] && { echo "?:?:?:?"; return; }
  total=$(echo "$line" | grep -oE '^[0-9]+')
  fail=$(echo "$line"  | grep -oE '[0-9]+ failed' | grep -oE '[0-9]+')
  pass=$((total - fail))
  echo "$total:$pass:$fail:0"
}

parse_pyaces() {
  local log="$1" p fc cr nr other total
  p=$(grep -oE 'PASS: [0-9]+' "$log" 2>/dev/null | grep -oE '[0-9]+')
  fc=$(grep -oE 'FAIL_CLEAN: [0-9]+' "$log" 2>/dev/null | grep -oE '[0-9]+')
  cr=$(grep -oE 'CRASH: [0-9]+' "$log" 2>/dev/null | grep -oE '[0-9]+')
  nr=$(grep -oE 'NO_RESULT: [0-9]+' "$log" 2>/dev/null | grep -oE '[0-9]+')
  [ -z "$p" ] && { echo "?:?:?:?"; return; }
  other=$(( ${cr:-0} + ${nr:-0} ))
  total=$(( p + ${fc:-0} + other ))
  echo "$total:$p:${fc:-0}:$other"
}

print_row() {
  local name="$1" result="$2" total pass fail other
  IFS=: read -r total pass fail other <<< "$result"
  printf "%-8s %8s %8s %8s %8s\n" "$name" "$total" "$pass" "$fail" "$other"
}

RESULT_ACES2=$(parse_aces2 "$ACES2_LOG")
RESULT_ACES3=$(parse_aces3 "$ACES3_LOG")
RESULT_PYACES=$(parse_pyaces "$PYACES_LOG")

echo "=== unified suite rollup ==="
printf "%-8s %8s %8s %8s %8s\n" "SUITE" "TOTAL" "PASS" "FAIL" "OTHER"
printf "%-8s %8s %8s %8s %8s\n" "-----" "-----" "----" "----" "-----"
print_row "aces2"  "$RESULT_ACES2"
print_row "aces3"  "$RESULT_ACES3"
print_row "pyaces" "$RESULT_PYACES"

STATUS=0
check_baseline() {
  local name="$1" actual="$2" expected="$3"
  if [ "$actual" != "$expected" ]; then
    echo "MISMATCH: $name actual=$actual expected(baseline)=$expected" >&2
    STATUS=1
  else
    echo "OK: $name matches baseline ($expected)"
  fi
}

echo
check_baseline "aces2"  "$RESULT_ACES2"  "$BASELINE_ACES2"
check_baseline "aces3"  "$RESULT_ACES3"  "$BASELINE_ACES3"
check_baseline "pyaces" "$RESULT_PYACES" "$BASELINE_PYACES"

exit $STATUS
