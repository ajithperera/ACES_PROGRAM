#!/bin/bash
# Thin build orchestrator for ACES_PROGRAM. Drives each subsystem's OWN,
# already-working build recipe with the shared toolchain from ./configure --
# does not reimplement or replace any of them:
#   aces2:  Makefiles/GNUmakefile via gmake
#   aces3:  autotools, configure_line.ufhpc + gmake
#   pyaces: pyaces/build.sh (f2py)
#
# aces2 must build before pyaces (pyaces links against aces2/lib/*.a).
# aces3 is fully independent of both.
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

declare -A RESULT

for t in $TARGETS; do
  echo "=== building $t ==="
  LOG="$ROOT/build_logs/unified_build_$t.log"
  case "$t" in
    aces2)
      ( cd "$ROOT/aces2" && \
        gmake FFLAGS_EXTRA=-fPIC F9XFLAGS_EXTRA=-fPIC CFLAGS_EXTRA=-fPIC CXXFLAGS_EXTRA=-fPIC all install \
      ) > "$LOG" 2>&1
      RESULT[$t]=$?
      ;;
    aces3)
      ( cd "$ROOT/aces3" && \
        bash configure_line.ufhpc && gmake clean && gmake all \
      ) > "$LOG" 2>&1
      RESULT[$t]=$?
      ;;
    pyaces)
      bash "$ROOT/pyaces/build.sh" > "$LOG" 2>&1
      RESULT[$t]=$?
      ;;
    *)
      echo "unknown target: $t (expected aces2, aces3, pyaces)" >&2
      exit 1
      ;;
  esac
  echo "  $t exit code: ${RESULT[$t]}  (log: build_logs/unified_build_$t.log)"
done

echo
echo "=== build summary ==="
FAIL=0
for t in $TARGETS; do
  if [ "${RESULT[$t]}" = "0" ]; then
    echo "  $t: OK"
  else
    echo "  $t: FAILED (exit ${RESULT[$t]})"
    FAIL=1
  fi
done
exit $FAIL
