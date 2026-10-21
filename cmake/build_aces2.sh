#!/bin/bash
# Builds aces2 (ACES II) via its own GNUmakefile-based build system, driven
# by the CMake superbuild. Does not reimplement aces2's build -- only loads
# the shared toolchain and invokes its existing `gmake ... all install`
# recipe (same flags phase-1's build.sh used).
#
# DIR_INSTBIN/DIR_INSTLIB are passed explicitly on the gmake command line so
# this always installs into THIS checkout's own aces2/bin,lib -- never the
# hardcoded main-tree default baked into Makefiles/GNUmakefile.src (see the
# 2026-08-12 comment there). GNU make passes command-line variable
# assignments to every recursive submake automatically via MAKEFLAGS, so one
# override here covers the whole (deeply recursive) build.
#
# 2026-08-12: aces2's own install logic (Makefiles/GNUmakefile's
# dir_writeable check) silently SKIPS installing a library/binary if
# DIR_INSTBIN/DIR_INSTLIB don't already exist on disk -- it never mkdirs
# them itself, and treats "directory missing" the same as "don't install
# here, fall through" rather than an error. On a fresh worktree these
# directories have never been created, so every single install silently
# no-ops and later links fail with "ld: cannot find -llb" etc. -- with no
# error at the point something actually went wrong. mkdir them first so the
# existing (unmodified) install logic actually engages.
set -euo pipefail

MODULES="$1"
ACES2_DIR="$2"
INSTALL_BIN="$3"
INSTALL_LIB="$4"

module purge
# shellcheck disable=SC2086
module load $MODULES

mkdir -p "$INSTALL_BIN" "$INSTALL_LIB"

cd "$ACES2_DIR"
gmake FFLAGS_EXTRA=-fPIC F9XFLAGS_EXTRA=-fPIC CFLAGS_EXTRA=-fPIC CXXFLAGS_EXTRA=-fPIC \
      DIR_INSTBIN="$INSTALL_BIN" DIR_INSTLIB="$INSTALL_LIB" \
      all install
