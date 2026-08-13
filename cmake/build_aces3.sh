#!/bin/bash
# Builds aces3 (ACES III) via its own existing autotools-generated
# per-directory Makefiles, driven by the CMake superbuild -- does not
# reimplement any of them. configure_line.ufhpc pins the compilers/flags
# (including the AMD-EPYC-safe -march=core-avx2, NOT -xcore-avx2 -- see
# commit 07938684) and is left untouched.
#
# 2026-08-12: aces3's own TOP-LEVEL `Makefile.in` ("FILES:= src
# src/aces/aces_sial test_compare tests", filtered to dirs that literally
# have their own Makefile) only ever reaches src/aces/aces_sial and
# test_compare -- "src" itself has no Makefile, so plain `gmake clean &&
# gmake all` from aces3/ NEVER builds any of the ~20 other library
# directories (framelib, sip1/2/3, the special_directory/* integral libs,
# a2driver's aces2/geopt/symcor, the sial compiler itself) OR src/main
# (xaces3 itself!). This was never caught before because every prior build
# in this project's history ran against a tree that already had those
# directories' outputs sitting in aces3/lib,bin from some earlier, more
# complete bootstrap (git-ignored, so invisible in the repo) -- a truly
# fresh checkout (this worktree) hits it immediately: `gmake all` fails at
# aces_sial's very first .sio compile with "../../../bin/sial: No such
# file or directory" because nothing ever built bin/sial.
#
# Rather than rewrite/restructure aces3's autotools files (real risk of
# subtly changing compilation flags none of this task is meant to touch),
# this supplies the missing top-level orchestration explicitly, in
# dependency order, confirmed via static analysis of each directory's own
# Makefile.in (grep for what library/binary each depends on) and
# cross-checked against aces3/nohup.out (a git-tracked historical build
# log from the original acesiii repo showing this same directory
# recursion). Every directory listed here still builds via its own
# unmodified Makefile -- nothing here compiles anything itself.
set -euo pipefail

MODULES="$1"
ACES3_DIR="$2"

module purge
# shellcheck disable=SC2086
module load $MODULES

cd "$ACES3_DIR"
# aces3/lib is tracked (via .gitkeep) but aces3/bin is not -- the SIAL
# compiler step's `cp sial ../../../bin/sial` assumes it already exists
# (same "must already exist, no mkdir of its own" pattern as aces2's
# install dirs) and fails with a plain "No such file or directory"
# otherwise. aces_sial/Makefile.in's SIO_DIR=../../../bin/sio has the
# identical assumption one level deeper: its rules are bare `cp $@
# $(SIO_DIR)` with no directory creation of their own, so when bin/sio
# doesn't already exist, cp silently creates a *plain file* named "sio"
# instead (each subsequent .sio compile just clobbers that same file) --
# no error, but every fast-suite case then fails at runtime with "Cannot
# open object file .../bin/sio/scf_rhf.sio" because the path component is
# a file, not a directory. Same latent bug as bin/sial, same reason it was
# never caught before this fresh worktree (see 2026-08-12 note above).
mkdir -p bin bin/sio lib

bash configure_line.ufhpc

# Level 0: independent leaf libraries -- each depends only on its own
# local sources (confirmed: none of their `all:` targets reference another
# aces3 directory's .a as a prerequisite).
LEVEL0_DIRS="
src/sia/framelib
src/sia/sip_shared
src/sia/sip/aces_instructions
src/sia/sip/aces_instructions_extnd
src/sia/sip/sip_instructions
src/sia/worker
src/sia/manager
src/aces/aces_library/special_directory/oed
src/aces/aces_library/special_directory/dup
src/aces/aces_library/special_directory/ecp
src/aces/aces_library/special_directory/nlo
src/aces/aces_library/special_directory/erd
src/aces/aces_library/special_directory/erd_F12
src/aces/aces_library/special_directory/oed_F12
src/aces/aces_library/special_directory/F12_grid
src/aces/aces_library/special_directory/dkh
src/aces/aces_library/special_directory/direct_integral
src/aces/aces_library/special_directory/prop_ints
src/aces/aces_library/a2driver/geopt
src/aces/aces_library/a2driver/aces2
src/aces/aces_library/a2driver/symcor
src/sial_compiler/sial
"

# Level 1: the SIAL compiler itself -- needs Level 0's libsial.a and
# libsip_shared.a. Produces bin/sial.
LEVEL1_DIRS="src/sial_compiler/compiler"

# Level 2: needs Level 1's bin/sial (aces_sial compiles .sial->.sio) or
# Level 0's libaces2.a (test_compare) -- these are the only two directories
# aces3's own top-level Makefile.in actually reaches, so `gmake -C .. all`
# here is the exact same recipe phase-1 used, just no longer run before its
# own prerequisites exist.
LEVEL2_DIRS="src/aces/aces_sial test_compare"

# Level 3: src/main (xaces3 itself) -- needs virtually everything above
# (including worker/manager's .o files, gathered via a plain `ls` at
# Makefile-parse time, so those two must already be built). Not reached by
# aces3's own top-level Makefile.in at all.
LEVEL3_DIRS="src/main"

for dir in $LEVEL0_DIRS $LEVEL1_DIRS $LEVEL2_DIRS $LEVEL3_DIRS; do
  gmake -C "$dir" clean
done

for dir in $LEVEL0_DIRS; do
  gmake -C "$dir" all
done
for dir in $LEVEL1_DIRS; do
  gmake -C "$dir" all
done
for dir in $LEVEL2_DIRS; do
  gmake -C "$dir" all
done
for dir in $LEVEL3_DIRS; do
  gmake -C "$dir" all
done
