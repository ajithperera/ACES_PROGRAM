# pyaces TODO

Work still to be done. When an item is finished, remove it from this file and record it in [DONE.md](DONE.md).
For background and rationale, see the [improvement plan](pyaces_improvements_en.pdf) and the [developer manual](pyaces_developer_manual.pdf).

## Conventions

- **Number:** `T-<n>`. A number is never reused (a finished item keeps its number in DONE.md).
- **Priority:** high / medium / low. **Size:** small (a few hours) / medium (about a day) / large (several days, needs a trial implementation).
- **Done when:** a change counts as finished only after both the local workstation and ccsd show aces2 143/143, pyaces 51/51 and a matching link audit. Add any item-specific criteria to the item.
- To add an item, copy the template below into the right stage.

```markdown
### T-00. Title
- **Priority / size:** medium / small
- **Problem:** what is wrong and why
- **To do:** the concrete work
- **Done when:** what has to be checked
- **Files:** `path`
```

## Summary

| No. | Item | Area | Priority | Size | Stage |
|---|---|---|---|---|---|
| T-01 | Per-module partial links that resolve names as the classic build does | correctness | high | large | 4 |
| T-08 | Incremental build | build | low | large | 4 |
| T-05 | Remaining known defects (Getlstpy, Get1ehpy, Runpccd) | correctness | medium | small–medium | any time |
| T-09 | Compiler and flag cleanup (mpicc, -qopenmp) | build | low | small | any time |
| T-10 | Fewer assumptions about the build environment (module, ifort check) | build | medium | small | any time |
| T-17 | API unit tests | testing | medium | medium | any time |

## Stage 1 — low risk, immediate payoff

The stage 1 items (T-03, T-04, T-07, T-16) were finished on 2026-10-23. See [DONE.md](DONE.md).

## Stage 2 — structural cleanup (no change in behaviour)

The stage 2 items (T-11, T-12, T-15) were finished on 2026-10-23. See [DONE.md](DONE.md).

## Stage 3 — usability and maintainability

The stage 3 items (T-06, T-13, T-14) were finished on 2026-10-23. The existing `aces2py` function names were kept, and the lower-case API is provided separately as the `python/pyaces` package. See [DONE.md](DONE.md).

## Stage 4 — large structural changes (decide after a trial implementation)

T-02 was finished together with T-14's `pyaces.run()`. See [DONE.md](DONE.md).

### T-01. Per-module partial links that resolve names as the classic build does
- **Priority / size:** high / large
- **Problem:** a single `.so` holding about 107 archives is held together by about 60 hand-written rename pairs. A new module or a toolchain change can introduce new collisions, and the audit baseline carries 296 known differences.
- **To do:** partially link each entry module with its own `ACES_LINK_LISTS` (`ld -r`) and hide every symbol other than the entry point with `objcopy`. Shared COMMON blocks stay global. Try it on one module (e.g. vcc) first.
- **Done when:** 51/51 with no rename list; the audit differences drop substantially; `.so` size and build time are compared and recorded.
- **Files:** `pyaces/build.sh`, `aces2/Makefiles/GNUmakefile` (`ACES_LINK_LISTS`)

### T-08. Incremental build
- **Priority / size:** low / large
- **Problem:** every build wipes `gen/` and `build/` and unpacks about 37 archives again.
- **To do:** split the steps and their dependencies with make or meson.
- **Done when:** rebuilding after a change to one glue file is much faster, and a clean build gives the same result.
- **Files:** `pyaces/build.sh`

## Any time — fold into other work

### T-05. Remaining known defects
- **Priority / size:** medium / small–medium
- **Problem:**
  - `Getlstpy` ignores `Irrep` and `Ndim` and uses `Left_index` as the irrep.
  - The symmetry-block offset in `Get1ehpy` (`get1ehpy.F:101`) looks wrong.
  - The runner calls `Init()` before `Runpccd`, but Runpccd calls Joda itself (never tested).
- **To do:** fix each one or document its meaning, and verify Get1ehpy on a molecule with symmetry.
- **Done when:** each function's result matches an independent calculation (e.g. the output of the classic executables).
- **Files:** `pyaces/a2_pyf/getlstpy.F`, `pyaces/a2_pyf/get1ehpy.F`, `pyaces/pyaces_run_case.py`

### T-09. Compiler and flag cleanup
- **Priority / size:** low / small
- **Problem:** the C wrappers are compiled with `mpicc` although they need no MPI. It is unclear whether `-qopenmp` on the glue code is needed with sequential MKL.
- **To do:** switch to `icx`/`gcc`, and check whether `-qopenmp` is needed and tidy it up.
- **Done when:** build and test results are unchanged.
- **Files:** `pyaces/build.sh` (lines 79–85), `build.sh` (top level)

### T-10. Fewer assumptions about the build environment
- **Priority / size:** medium / small
- **Problem:** `pyaces/build.sh` always runs `module purge/load`, so it fails in a setvars/path environment. The top-level `build.sh` checks for `ifort` before the `ifort` shim is on PATH, so it rejects an ifx-only toolchain.
- **To do:** follow configure's `ACES_PROGRAM_ENV_MODE`, and move the `ifort` check after the shim or accept ifx.
- **Done when:** modules mode (local and ccsd) gives the same results, and a setvars-mode build succeeds.
- **Files:** `pyaces/build.sh` (lines 1–17), `build.sh` (top level, around lines 215–225), `configure`

### T-17. API unit tests
- **Priority / size:** medium / medium
- **Problem:** only energies are checked today.
- **To do:** add small tests for reading and writing records, `Getreclenpy`/`Cgetrecpy`, each branch of `Runhbar`, and list access, and include them in `run_tests.sh pyaces`.
- **Done when:** the tests pass on the current code, and the relevant tests fail on the code before f97a9759.
- **Files:** new `pyaces/tests/` (tentative), `run_tests.sh`
