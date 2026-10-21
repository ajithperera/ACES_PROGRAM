# CMake superbuild -- status notes

Branch: `worktree-cmake-superbuild` (pushed to `origin/worktree-cmake-superbuild`)
Worktree: `.claude/worktrees/cmake-superbuild`

## Status: DONE

aces3's CMake integration is verified end-to-end and committed. All three
subsystems (aces2, aces3, pyaces) build cleanly from a truly fresh
checkout via the CMake superbuild, and all three fast test suites match
the zero-point-20260812 baseline exactly:

| suite  | total | pass | fail | other |
|--------|-------|------|------|-------|
| aces2  | 125   | 121  | 1    | 3     |
| aces3  | 37    | 37   | 0    | 0     |
| pyaces | 58    | 37   | 16   | 4 crash + 1 no_result |

`ctest --test-dir build` (the `rollup` test) reports 100% pass, 4/4.

## History this session

Commit `f0d9c0d3` (prior session) added the CMake superbuild itself and
fixed 6 pre-existing aces2 link/config bugs, discovered because this
worktree is a genuinely from-scratch checkout with none of the stale
build residue (`bin/`, `lib/` contents) that every prior build in this
project's history had accidentally relied on.

This session found and fixed two more bugs in the same class, both in
aces3, both invisible for the same reason (latent since forever, only
surfaced by a from-scratch build):

1. **`aces3/configure.ac` was missing `AC_CONFIG_FILES` for
   `prop_ints/Makefile`.** `prop_ints` is real integral-routine code (not
   a stub), and every other subsystem Makefile.in already links
   `-lprop_ints` -- but `./configure` never generated its Makefile, so
   `libprop_ints.a` never got built, and the final link failed with an
   undefined-symbol error downstream. Fixed by adding the
   `AC_CONFIG_FILES` line and regenerating `configure` via autoreconf.
   `cmake/build_aces3.sh` also needed `prop_ints` added to its `LEVEL0_DIRS`
   leaf-library build list (this part landed in the prior session, just
   uncommitted until now).

2. **`aces_sial/Makefile.in`'s `SIO_DIR=../../../bin/sio` rules assume
   the directory already exists** -- they're bare `cp $@ $(SIO_DIR)` with
   no `mkdir` of their own (same pattern as the already-fixed `bin/sial`
   case). When `bin/sio` doesn't pre-exist as a directory, `cp` silently
   creates a **plain file** named `sio` instead of failing loudly -- so
   the *build* reports success, but every fast-suite case then fails at
   *runtime* with `Cannot open object file .../bin/sio/scf_rhf.sio`
   because that path component is a file, not a directory. This is why
   the first full-suite run (job 39346718) showed aces2/pyaces matching
   baseline but aces3 at 0/37 (all FAIL). Fixed by adding `bin/sio` to
   the `mkdir -p` line in `cmake/build_aces3.sh`.

Both fixes are documented inline (dated comments) in
`cmake/build_aces3.sh` and `aces3/configure.ac`.

## Verification jobs (sbatch, this session)

- `39294168` -- full build (aces2+aces3+pyaces), confirmed aces3 built
  clean with the prop_ints fix (`ACES3_BUILD_EXIT: 0`), but no test suite
  had been run against it yet.
- `39346718` -- first full `ctest` run. aces2 and pyaces matched
  baseline immediately. aces3 came back 0/37 (all FAIL) -- this is what
  surfaced the `bin/sio` bug above.
- `39347164` -- after the `bin/sio` fix: rebuilt aces3, reran full
  `ctest`. All three suites + rollup pass, exact baseline match (table
  above). aces3 suite alone took ~36 min wall time (real physics runs,
  not the ~8s of instant failures from the broken run) -- much longer than
  aces2/pyaces, worth knowing for future turnaround expectations.

Logs: `build_logs/cmake_build.39294168.out`,
`build_logs/cmake_tests.39346718.out`,
`build_logs/cmake_aces3_retest.39347164.out`.

## Uncommitted / left as-is (intentional)

- `aces3/config.log`, `aces3/config.status` -- pure `./configure`
  regeneration byproducts (hostnames, PATH ordering, flag whitespace).
  Diffed by hand to confirm no real content changed. Matches this
  project's established pattern of leaving these two files out of
  version control.
- `aces3/src/aces/aces_library/a2driver/aces2/make.err`,
  `.../make.out` -- stray compiler-warning capture files from the
  aces_sial build's `1>>make.out 2>>make.err` redirects. Harmless, not
  gitignored, not committed (same treatment as other build byproducts in
  this repo).

## What's next (not blocking, just for whoever picks this up)

- No further aces3-specific work is needed -- the CMake superbuild's
  aces3 leg is done and verified.
- The two latent bugs found this session (and the earlier 6 in aces2)
  all share one root cause: **every build in this project's history
  before the CMake superbuild ran against a tree with stale bootstrap
  residue** (built libs/bins left over from some earlier, more manual
  build) that masked missing directory-creation and missing
  `AC_CONFIG_FILES` entries. It's plausible more of these exist elsewhere
  in aces3's ~20 build directories, only reachable by scenarios this
  fast suite doesn't exercise (e.g. slow/full test suite, or a build
  option combination not covered by `ACES_BUILD_ACES3=ON` default). If a
  *different* from-scratch build ever hits a new "no such file or
  directory" or silent-clobber failure, check first whether it's the
  same "directory assumed to pre-exist" pattern before assuming it's
  something else.
