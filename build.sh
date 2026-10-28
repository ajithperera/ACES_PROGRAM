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
# Defaults for configs written by an older ./configure (which had neither knob).
ACES_PROGRAM_MACHNAME="${ACES_PROGRAM_MACHNAME:-}"
ACES_PROGRAM_ACES3_CONFIG="${ACES_PROGRAM_ACES3_CONFIG:-configure_line.ufhpc}"
# pyaces/build.sh does its own `module purge` + load, so it has to be told
# which module set this build is using rather than inheriting the loaded one.
export ACES_PROGRAM_MODULES

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
    echo "  cannot build $target -- not on PATH after loading the configured toolchain:$missing" >&2
    return 1
  fi
  return 0
}

# aces2 has no "did the build succeed" signal of its own: with -k (needed, see
# the aces2 branch) gmake always ends nonzero because of the modules that never
# build, while without it a half-finished tree can still look fine. What the
# tree does record is the enumerated list of built binaries in aces2/.gitignore
# -- explicit paths, no globs, by that file's own stated convention. Diff it
# against what is actually installed in aces2/bin.
#
# xrspdens1 is excluded: ACES_LINK_LISTS names an rspdens1 module, but no
# rspdens1/ directory exists in this repo (it is not in git history either), so
# that .gitignore line is aspirational and nothing can ever satisfy it.
aces2_missing_binaries() {
  comm -23 \
    <(grep -oE '^[a-z0-9_]+/x[a-z0-9_]+$' "$ROOT/aces2/.gitignore" \
        | cut -d/ -f2 | grep -vx xrspdens1 | sort -u) \
    <(ls "$ROOT/aces2/bin" 2>/dev/null | grep '^x' | sort -u)
}

# aces2's makefile makes each module's binary depend on that module's OWN
# library and nothing else -- not on the shared libraries it links against
# (acescore, libr, tools, ...). Editing one of those therefore rebuilds the
# archive and leaves all ~76 executables in place, silently linked against the
# old code: the build reports OK, the binary set is complete, and the suite
# quietly tests yesterday's build. (Hit for real while fixing the BW-CC
# file-name bug in acescore/aces_io_init.F -- libacescore.a had the fix, xvcc
# was five hours old, and the case still failed exactly as before.)
#
# A from-scratch build trips this too, and legitimately: modules link in
# alphabetical order, so an early module's binary can be linked before a later
# module's library that it needs even exists (pdens needs libvmol2ja.a).
#
# List the binaries older than the newest installed library.
aces2_stale_binaries() {
  local newest_lib lib_mtime b
  newest_lib=$(ls -t "$ROOT/aces2/lib"/*.a 2>/dev/null | head -1)
  [ -n "$newest_lib" ] || return 0
  lib_mtime=$(stat -c %Y "$newest_lib" 2>/dev/null) || return 0
  for b in "$ROOT/aces2/bin"/x*; do
    [ -f "$b" ] || continue
    if [ "$(stat -c %Y "$b" 2>/dev/null || echo 0)" -lt "$lib_mtime" ]; then
      basename "$b"
    fi
  done
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
      # MACHNAME is passed on the command line (not exported) so it overrides
      # the GNUmakefile's own `uname -n` assignment and reaches every
      # recursive sub-make; empty means "don't pass it", i.e. HiPerGator's
      # own hostname lookup, unchanged.
      MACHNAME_ARG=""
      [ -n "$ACES_PROGRAM_MACHNAME" ] && MACHNAME_ARG="MACHNAME=$ACES_PROGRAM_MACHNAME"
      if require_tools aces2 gmake; then
        : > "$LOG"
        # A fresh checkout has no aces2/bin or aces2/lib (both gitignored) --
        # they are created by aces2's own xprep, and every module's install
        # step writes into them. Without it the first module that links dies
        # with "ld: cannot find -llb". xprep is idempotent but also chmods the
        # whole tree, so only run it when the directories are actually absent.
        ( cd "$ROOT/aces2" && { { [ -d lib ] && [ -d bin ]; } || sh ./xprep -f; } ) >> "$LOG" 2>&1
        # Why -k, and why more than one pass:
        #  * 6 module directories (ccsdtqf, copt, oombpt2, pdens, psi4dbg,
        #    tests) do not build in this tree and never have -- none of them
        #    appears in aces2/.gitignore's enumerated list of built binaries.
        #    Without -k the first one to fail aborts the other ~110 modules,
        #    which is only survivable on a tree whose lib/ is already full
        #    from some earlier build.
        #  * The makefile carries no .mod dependency graph, so a module whose
        #    Fortran modules are defined in a later-compiled file fails on one
        #    pass and succeeds on the next (alice_nwchem and get_acesmo each
        #    need two).
        # So iterate until the set of missing binaries stops shrinking.
        MISSING=""
        PREV_N=-1
        for pass in 1 2 3; do
          echo "=== aces2 build pass $pass ===" >> "$LOG"
          # shellcheck disable=SC2086
          ( cd "$ROOT/aces2" && \
            gmake -k $MACHNAME_ARG FFLAGS_EXTRA=-fPIC F9XFLAGS_EXTRA=-fPIC CFLAGS_EXTRA=-fPIC CXXFLAGS_EXTRA=-fPIC all install \
          ) >> "$LOG" 2>&1
          MISSING=$(aces2_missing_binaries)
          N=$(echo "$MISSING" | wc -w)
          echo "=== aces2 pass $pass: $N expected binaries still missing ===" >> "$LOG"
          [ "$N" = "0" ] && break
          [ "$N" = "$PREV_N" ] && break
          PREV_N=$N
        done
        # Relink whatever the shared-library dependency gap left behind (see
        # aces2_stale_binaries). One pass is enough: relink is binclean + all +
        # install, so afterwards every binary post-dates every library, and an
        # otherwise-unchanged rebuild finds nothing to do here.
        STALE=$(aces2_stale_binaries)
        if [ -n "$STALE" ]; then
          echo "=== aces2: $(echo "$STALE" | wc -w) binaries older than the libraries -- relinking ===" >> "$LOG"
          # shellcheck disable=SC2086
          ( cd "$ROOT/aces2" && \
            gmake -k $MACHNAME_ARG FFLAGS_EXTRA=-fPIC F9XFLAGS_EXTRA=-fPIC CFLAGS_EXTRA=-fPIC CXXFLAGS_EXTRA=-fPIC relink \
          ) >> "$LOG" 2>&1
          MISSING=$(aces2_missing_binaries)
          STALE=$(aces2_stale_binaries)
          [ -n "$STALE" ] && echo "  aces2: still older than the libraries after relink: $STALE" >&2
        fi
        if [ -z "$MISSING" ]; then
          RESULT[$t]=0
        else
          echo "  aces2: expected binaries still missing after $pass pass(es): $MISSING" >&2
          RESULT[$t]=1
        fi
      else
        RESULT[$t]=1
      fi
      ;;
    aces3)
      # Which MPI wrappers to check for depends on which configure wrapper is
      # selected (they differ by site: openmpi's mpif77/mpicc/mpicxx on
      # HiPerGator, Intel MPI's mpiifort/mpiicx/mpiicpx locally), so read the
      # names out of the file that is actually about to run.
      A3CONF="$ROOT/aces3/$ACES_PROGRAM_ACES3_CONFIG"
      if [ ! -f "$A3CONF" ]; then
        echo "  cannot build aces3 -- no such configure wrapper: $A3CONF" >&2
        RESULT[$t]=1
        echo "  $t exit code: ${RESULT[$t]}  (log: build_logs/unified_build_$t.log)"
        continue
      fi
      A3TOOLS=$(grep -oE '\b(F77|FC|CC|CXX)=[^ ]+' "$A3CONF" | cut -d= -f2 | sort -u | tr '\n' ' ')
      # shellcheck disable=SC2086
      if require_tools aces3 gmake $A3TOOLS; then
        # aces3/bin (and bin/sio, where the compiled .sio jobflows land) is
        # gitignored, so a fresh checkout does not have it -- the sial
        # compiler's own "cp sial ../../../bin/sial" then fails. lib/ gets
        # created by the build itself; bin/ does not. aces2 has xprep for
        # exactly this; aces3 has no equivalent, so create them here.
        ( cd "$ROOT/aces3" && \
          mkdir -p bin/sio lib && \
          bash "$ACES_PROGRAM_ACES3_CONFIG" && gmake clean && gmake all \
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
