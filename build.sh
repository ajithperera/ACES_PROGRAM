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
if ! module load $ACES_PROGRAM_MODULES; then
  echo "module load failed: $ACES_PROGRAM_MODULES" >&2
  echo "This tree's toolchain is Intel compilers + Intel MKL + MPI, provided by those" >&2
  echo "modules on HiPerGator. If this host names them differently, re-run" >&2
  echo "  ./configure --modules=\"...\"" >&2
  exit 1
fi

# `module load` succeeding is not by itself proof the toolchain is usable: a
# differently-named module set can load cleanly and still leave no Intel
# compiler on PATH. None of the three subsystems reports that as the cause --
# they die deep inside gmake with a per-source-file "ifx: command not found"
# or an empty -L$MKLROOT link line instead -- so check up front, per target,
# against exactly what that target's own build recipe hardcodes.
require_tools() {
  local target="$1"; shift
  local missing=""
  local tool
  for tool in "$@"; do
    command -v "$tool" >/dev/null 2>&1 || missing="$missing $tool"
  done
  if [ -n "$missing" ]; then
    echo "  cannot build $target -- not on PATH after 'module load $ACES_PROGRAM_MODULES':$missing" >&2
    return 1
  fi
  return 0
}

mkdir -p "$ROOT/build_logs"

TARGETS="$*"
[ -z "$TARGETS" ] && TARGETS="aces2 aces3 pyaces"

declare -A RESULT

for t in $TARGETS; do
  echo "=== building $t ==="
  LOG="$ROOT/build_logs/unified_build_$t.log"
  case "$t" in
    aces2)
      # Makefiles/GNUmakefile selects its own compilers from its MACHNAME map
      # (and errors out by itself, clearly, if this host has no entry there),
      # so gmake is the only tool to check for here.
      if require_tools aces2 gmake; then
        ( cd "$ROOT/aces2" && \
          gmake FFLAGS_EXTRA=-fPIC F9XFLAGS_EXTRA=-fPIC CFLAGS_EXTRA=-fPIC CXXFLAGS_EXTRA=-fPIC all install \
        ) > "$LOG" 2>&1
        RESULT[$t]=$?
      else
        RESULT[$t]=1
      fi
      ;;
    aces3)
      # configure_line.ufhpc hardcodes these MPI wrappers.
      if require_tools aces3 gmake mpif77 mpicc mpicxx; then
        ( cd "$ROOT/aces3" && \
          bash configure_line.ufhpc && gmake clean && gmake all \
        ) > "$LOG" 2>&1
        RESULT[$t]=$?
      else
        RESULT[$t]=1
      fi
      ;;
    pyaces)
      # pyaces/build.sh hardcodes ifort/ifx/mpicc and links against
      # $MKLROOT/lib/intel64 -- an unset MKLROOT silently becomes "-L/lib/intel64".
      if require_tools pyaces ifort ifx mpicc python3 && [ -n "${MKLROOT:-}" ]; then
        bash "$ROOT/pyaces/build.sh" > "$LOG" 2>&1
        RESULT[$t]=$?
      else
        [ -n "${MKLROOT:-}" ] || echo "  cannot build pyaces -- MKLROOT is unset (no Intel MKL in the loaded modules)" >&2
        RESULT[$t]=1
      fi
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
