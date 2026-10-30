# pyaces DONE

Record of finished work. New entries go at the top (newest first).
Items that came from [TODO.md](TODO.md) keep their number (`T-<n>`); work that was not in TODO is recorded as `D-<date>`.

## Conventions

```markdown
### T-00. Title (YYYY-MM-DD)
- **Commits:** `hash` (all of them, if several)
- **What changed:** what was changed and how
- **Verification:** local / ccsd results (aces2, pyaces, link audit), plus any extra checks
- **Open / notes:** follow-up work or caveats (omit if none)
```

---

### T-17. API unit tests (2026-10-30)
- **Commits:** (not committed yet)
- **What changed:**
  - New `pyaces/tests`. `api_cases.py` holds 8 tests and their ZMATs; `run_api_tests.py` runs each test in a fresh process and directory and prints a Summary block. `run_tests.sh pyaces` runs them after the suite (log `build_logs/unified_test_pyaces_api.log`) as a separate `py-api` roll-up row.
  - Tests:
    - `records`: `Getreclenpy`, missing records, float64/int64/character round trips, reopen.
    - `scf_rhf`, `scf_uhf`: energies against xaces2, and reading vscf's character record `EVCSYMAF`.
    - `hbar_from_scratch`, `hbar_after_scf`: CCSD energy and the `HBAR` record.
    - `lists`: list 46 by irrep and column vs `get_all_irreps`, error paths, `put_list`.
    - `one_electron_c2v`, `one_electron_c1`: `CᵀSC = 1`, block structure with an empty irrep, and E(SCF) rebuilt from `oneh` against xaces2.
  - The reference energies are from xaces2, not pyaces.
  - Runhbar's "CC done" branch is not tested: it needs a `CC_DONE` record without `SCFEVCA0`, and no ACES II program writes `CC_DONE`, so the API cannot reach it.
  - The tests found two character-record defects, now fixed:
    - `Cputrecpy` declared `Source` as `Character*8`, so only 8 characters were stored and the rest of the record was whatever followed in memory. It is now assumed-length, and `put_string_record` passes the bytes directly.
    - pyaces was compiled without aces2's `-DF_64BIT`, so `a2_getcrec`/`a2_putcrec` counted 4-byte words. `Cgetrecpy` on a record written by aces2 (e.g. `EVCSYMAF`) asked for twice its length and ended the process ("record length mismatch"). `FFLAGS` now has `-DF_64BIT`; in pyaces it only affects those two files.
- **Verification:**
  - Local: py-api 8/8; audit 244; pyaces 60/60 with every `got=` value unchanged.
  - With `f97a9759` reverted in a scratch copy, `hbar_after_scf` crashes (SIGSEGV), and the six tests that call `record_length` fail on `Getreclenpy`. The build's `--calls` check also reports the two `hbar_` calls.
  - `scf_uhf` still passes on that copy: the 4-byte `Iuhf` happened to arrive as 1, so the test cannot show that fix.
  - The new character-record checks fail on the code before this change.

### T-10. Fewer assumptions about the build environment (2026-10-30)
- **Commits:** `442efa80`
- **What changed:**
  - The top-level half was already done by Ajith's `f444e634`. `build.sh` checks for configure's `ACES_PROGRAM_FC`, and the ifort→ifx shim is gone; in setvars/path mode it passes `PYACES_MODULES=""` after loading the toolchain itself.
  - Left: `pyaces/build.sh` run on its own ignored `build_config.sh` and loaded the HiPerGator module set. On this workstation (setvars mode, but `module` exists) that failed with "Unable to locate a modulefile for 'intel/2025.1.0'".
  - It now takes, in order: `PYACES_MODULES` when set (CMake superbuild, `../build.sh`); a non-empty `ACES_PROGRAM_MODULES`; `build_config.sh`'s `aces_program_load_env` with its `ACES_PROGRAM_FC`/`ACES_PROGRAM_CC`; the HiPerGator set only without a `build_config.sh`. A config from an older `./configure` (no loader) still gets its module set.
  - The build log has a `compilers:` line.
- **Verification:**
  - Local, setvars mode: `bash pyaces/build.sh` on its own, `./build.sh pyaces`, and `PYACES_MODULES=""` with setvars loaded all build. The standalone run used ifort and icx 2024.0. Audit 244 each time.
  - pyaces 60/60 for the standalone and top-level builds, every `got=` value identical to the run before T-09.
  - ccsd, modules mode (ifx/icx 2023.2, PBS 8128): `bash pyaces/build.sh` on its own after `module purge` loads configure's modules itself, uses ifx and icx, audit 268, pyaces 60/60. Then `./build.sh pyaces`: audit 268, aces2 158/158, pyaces 60/60. Both runs give every `got=` value identical to the T-09 ccsd run.

### T-09. Compiler and flag cleanup (2026-10-30)
- **Commits:** `0213f6a7`
- **What changed:**
  - `pyaces/build.sh` compiles the f2py C module and `fortranobject.c` with `$PYACES_CC` instead of `mpicc`; neither uses MPI. The top-level `build.sh` passes configure's `ACES_PROGRAM_CC` (icx locally and on ccsd) and checks for it instead of `mpicc`. Standalone, the default is icx, else gcc. Locally `mpicc` was Intel MPI's wrapper around gcc.
  - `-qopenmp` is dropped from the Fortran flags. The pyaces sources have no OpenMP directives, aces2 is compiled without it (`-Zp8 -zero -assume dummy_aliases -i8 -O2`), and MKL is linked sequential. Its only effect on the glue code was to make every local variable automatic.
- **Verification:**
  - Local (ifort 2024.0, icx 2024.0): audit 244 with no NEW entries; pyaces 60/60.
  - All 60 `got=` values are identical to the run before the change.
  - The .so's `NEEDED` libraries are the same before and after, and it has no OpenMP symbols either way.
  - ccsd (ifx 2023.2, icx, PBS 8124): audit 268 with no NEW entries, aces2 158/158, pyaces 60/60, all 60 `got=` values identical to the previous ccsd run; the T-05 checks give the same numbers as before.
  - icx reports 7 `-Wformat` warnings (long vs long long) in the f2py-generated C; the two types have the same size on x86-64 Linux.

### T-05. Remaining known defects: Getlstpy and Get1ehpy (2026-10-30)
- **Commits:** `dc9b4cd7`
- **What changed:**
  - `Get1ehpy` returned zeros for every molecule. `Aces_malloc` can return a negative `J0` (the heap lies below the anchor `Jcore`), which the routine took for a failed allocation, so it printed "Unable to allocate" and returned before reading IIII. It now tests `Jcore(1) = 0`, as `aces2/acescore/aces_malloc.c` says to. With symmetry the irrep blocks also went to wrong offsets (`Index = Index + Nelms*(K-1)`, and the triangle offsets were one too far); it now returns the full block-diagonal matrix in the SO basis, functions in `NUMBASIR` order. Unknown names are rejected instead of reading an unset label.
  - `Getlstpy`'s arguments now have documented meanings: `Left_index` is the irrep of the columns and `Right_index` the list number (`MOIO(irrep, list)`, as in ACES II's `getlst`); `Irrep` is unused and kept for existing callers. `Getlstpy` and `Putlstpy` check that the list exists, the column range and `Ndim`, and print a message and do nothing instead of overrunning the array. The new f2py function `Getlstdimpy(left, right)` returns a list's rows and columns (27 exposed functions).
  - Python: `list_shape(lst, irrep)`; `get_list` sizes its array by default and, with `put_list`, raises `KeyError`/`ValueError` on a bad request; `get_1e_ints(name)` wraps `Get1ehpy` as an `(nbas, nbas)` array.
- **Verification:**
  - Local (ifort 2024.0): H2O 6-31G* CCSD in C2v (NUMBASIR 9/3/5/1) and with `SYMMETRY=OFF`. In both, `CᵀSC = 1` with `SCFEVCA0` to 5×10⁻¹⁵. The (H, S) generalized eigenvalues agree between the two runs, and so does the one-electron energy 2 tr(DH) = −123.0176492438 (to 10⁻¹¹; 10.000 electrons). List 46 read irrep by irrep equals `Getallpy`, and the T2 sum of squares agrees between the runs. The error paths and a `put_list` round trip were also checked.
  - Audit 244 with no NEW entries; pyaces 60/60.
  - ccsd (ifx 2023.2, PBS 8123): pyaces rebuilt, audit 268 with no NEW entries, aces2 158/158, pyaces 60/60; the same checks give `CᵀSC = 1` to 1×10⁻¹⁴ and the same one-electron energy in C2v and C1 (to 10⁻¹¹).
- **Open:** `get1epy.F` (linked, not exposed, never called) is an old copy of Get1ehpy with the same defects. The checks above are a scratch script; T-17 should turn them into tests.

### D-2026-10-30b. Runpccd fixed and a pCCD test case added (2026-10-30; the Runpccd part of T-05)
- **Commits:** `11762bef`
- **What changed:**
  - `Runpccd` read ZMAT again with its own `Joda_`/`aces_init` although the worker calls `Init()` first ("job archive subsystem is already initialized"), and called `Ints_()` with no arguments. It now runs after `Init()` like `Runcc` and stops at `Vcc`, as xrunpccd's HF/KS branch and xaces2 do for single-point pCCD. `Vlambda` and `Pccd` were dropped: `Vlambda` stopped on "List (3,144) does not exist".
  - New aces2 fast case `pccd_h2o_sto3g` (H2O STO-3G, `CALC=PCCD`, HF reference; reference −74.988498185196 from xaces2). It is the first case the TSV routes to `Runpccd`.
  - Dropping the two calls changed the link: hcmult's `check_t2_` started winning for vcc, then libc4's (whose `check_t2.o` also defines `checksum_`, clashing with lcct's). Both are renamed in `tools/localize.txt` (21st instance), leaving trp's as in the classic links. The new `symoff_` common-block shadow is latent (its only callers, dens's `intov5`/`xint7`, are not in the .so) and is recorded in `common_shadow_baseline.txt`.
- **Verification:**
  - Local (ifort 2024.0): `pccd_h2o_sto3g` passes in aces2 and pyaces (difference 4×10⁻¹²), pyaces 60/60, audit 244 with no NEW entries.
  - ccsd (ifx 2023.2, PBS 8116): pyaces rebuilt, audit 268 with no NEW entries, aces2 158/158, pyaces 60/60.
- **Open:** orbital-optimized pCCD (`PROGRAM=PCCD`) remains out of reach for pyaces; it is xrunpccd's own loop.

### D-2026-10-30a. Suite classifier repaired and TSV regenerated (2026-10-30)
- **Commits:** `847b1aeb`
- **What changed:** `classify_suite.py --check` stopped with a ValueError on the bare-name `test_results` entries upstream added on 2026-10-19 (`pccd_oo_rfo_*`); the parser now reads the first field only, as aces2's runner does. `PROGRAM=PCCD` cases are excluded (xaces2 hands them to xrunpccd). The regenerated TSV covers all 157 aces2 cases (32 were missing) and adds eight H2/H2O CCSD cases to the pyaces suite. The manual (chapter 2, 9.3, 10.2) and TODO no longer say to edit the TSV by hand.
- **Verification:** local pyaces 59/59; ccsd (PBS 8114) pyaces 59/59.

### T-02. An ACES error terminates the Python interpreter (2026-10-23)
- **Commits:** `ba6c8646` (together with T-14)
- **What changed:** `pyaces.run()` runs the calculation in a subprocess (`pyaces_run_case.py`) and returns a `Result` on success, or a `CalculationError` (status, returncode, workdir, log_tail) on failure (`aces_exit`, crash, timeout). `pyaces_run_case.py` now honours `PYACES_EXT_DIR`.
- **Verification:** a ZMAT changed to `CALC=NOSUCHMETHOD`: both local and ccsd get a `CalculationError` (CRASH, exit 1) and the calling process keeps running. The normal case 026 gives −76.27254545743 (within 7×10⁻¹¹ of the reference).
- **Open:** functions used in-process (`pyaces.init()` and so on) still end the process on an ACES II error. That follows from how `aces_exit` works and will not change; the documentation says so.

### T-14. Python package and API layer (2026-10-23)
- **Commits:** `ba6c8646`
- **What changed:**
  - Added `pyproject.toml` and the `pyaces` package under `pyaces/python/`. The existing `aces2py` names are unchanged (still available through `pyaces.ext()`), with lower-case functions layered on top.
  - `drivers.py`: `init()`, `run_scf()`, `run_cc()` and so on. `records.py`: `record_length()`, `get_record()` (length found automatically, `KeyError` for a missing record), `put_record()`, character records. `lists.py`: list access with argument names that match what they actually mean. `runner.py`: `run()` (T-02). `moio.py`: the old `Python/moio_len.py`, moved.
  - `aces2py.so` is found in `pyaces/Main/build` or `$PYACES_EXT_DIR`.
- **Verification:**
  - In a venv, `pip install --no-build-isolation -e pyaces/python`, then import succeeds.
  - Case 026 in-process: `init → run_scf → run_hbar`, SCF −76.0462837086, CCSD −76.27254545743, `record_length("NOSUCHRC")` → None, `get_record` on a missing record → KeyError, real and character records read back what was written.
  - ccsd (PBS 7717): `run("cc", 026)` −76.27254545743; the failing ZMAT gives `CalculationError`.
- **Open:** the list functions (`get_list` and so on) were not tested (to be covered by T-17).

### T-13. Rename list moved into a data file (2026-10-23)
- **Commits:** `ba6c8646`
- **What changed:** `EXTRA_LOCALIZE_SYMBOLS` (60 pairs) and `EXTRA_LOCALIZE_MODULES` (21 modules) became `rename <module> <symbol>  # reason` and `extract <module>` lines in `tools/localize.txt`, and the roughly 290 lines of discovery-history comments in build.sh moved to that file too. build.sh (about 510 → 247 lines) and the audit tool read the same file. The pairless `extract libr2` was kept (re-packing the archive would change member order and could change the link result).
- **Verification:** the two lists read from the file equal the old strings. `nm` output of all 37 newly built patched archives is identical to before. Audit 296 (matches).

### T-06. f2py module no longer generated twice (2026-10-23)
- **Commits:** `ba6c8646`
- **What changed:** f2py is no longer given `Main/aces2.pyf` and builds the module from the sources alone; `aces2.pyf` was deleted. Signatures come only from the `Cf2py` directives in the sources.
- **Verification:** `aces2pymodule.c` built without the `.pyf` is byte-identical to before, and the f2py log shows the module generated once. The 26 exposed functions are unchanged. Local and ccsd: aces2 125/125, pyaces 51/51, audit 296, `--calls` 0.
- **Notes:** the unused COMMON initialisation wrappers (`f2pyinit*` in `aces2py-f2pywrappers.f`) are gone. Because the Python names were kept, the `sed` lower-casing (first argument of `F_FUNC`) stays. f2py `--lower` would remove it but would change the Python names.

### T-15. `Python/moio_len.py` fixed (2026-10-23)
- **Commits:** `fabd383a`
- **What changed:** the broken `"key" = value` dictionaries became lists of `(label, number)` pairs, `MO_DISTRIBUTIONS` (21) and `MO_LISTS` (93), with lookup functions `distribution_numbers()`, `list_numbers()` and `list_labels()`. The values were left as in the original.
- **Verification:** imports locally and on ccsd. Example: `list_numbers("T2(Ij,Ab):Ab,Ij")` → `[46]`.
- **Later change:** moved to `python/pyaces/moio.py` (`pyaces.moio`) in T-14.
- **Open:** duplicates in the original were not fixed but are listed in the module docstring: distribution "ia" appears as both 14 and 16, "HB(MI):MI" as 911 and 912 (912 is probably "HB(mi):mi"), and "HB(Mb,Ej):MbEj" overlaps 56–59. The values were not checked one by one against the ACES II sources. Lists 7–10 and 27–30 being shared by W and HB is expected, since HBAR is built on top of the integral lists.

### T-12. Duplicate header copies removed (2026-10-23)
- **Commits:** `fabd383a`
- **What changed:** deleted the nine `Main/*.com` files, so `init.F` and `aces2_main.F` read aces2's `include` and `acescore/include` headers like every other pyaces source.
- **Verification:**
  - The COMMON block sizes in the compiled `init.o` match a build against the aces2 headers (`nm -S`). `aces2_main.o` is identical before and after.
  - Local and ccsd: aces2 125/125, pyaces 51/51, audit 296, `--calls` 0.
- **Notes:** the old copies differed in layout, not just size. The `lists.com` copy declared `/lists/` as 0x61a80 bytes with `maxref=1`, against 0x64b540 bytes in the aces2 header, and `io_ptrs` was 0xa0 against 0xa50. `init.F` only includes the headers and uses none of the variables, so there was no actual damage. The manual's statement that "the sizes are the same, so it is harmless" was wrong and has been corrected.

### T-11. Legacy files cleaned up (2026-10-23)
- **Commits:** `fabd383a` (state before the deletion: tag `pyaces-legacy-20261023` = `4e53d532`)
- **What changed:** deleted 6,415 files: `Fxx/`, `include/`, `a2_pyf/include_no/`, `Main/Archive/`, `a2_pyf/Archive/`, `Main/Pilot_tests/`, `Exp/`, `rename.py`, the old `.f` copies (35 `a2_pyf/*.f`, `Main/init.f`, `Main/aces2_main.f`) and `suite_classification.tsv.before_*`. Also deleted, though not in TODO, the parts of the same old build: `Main/compile_acespy`, `Main/GNUmakefile*` (symbolic links into `Fxx`) and `a2_pyf/GNUmakefile*`. The pyaces directory went from about 47 MB to 2.2 MB (excluding build output).
- **Verification:** no script refers to the deleted paths (`git grep`). Build and test results are as for T-12.
- **Notes:** `README`, `frequently_asked_questions` and `LICENSE` were kept. Chapter 2 and appendix A of the manual were updated. The tag was pushed to both ccsd and GitHub.

### T-16. Missing test summary and classifier problems (2026-10-23)
- **Commits:** `95cc01a2`
- **What changed:**
  - The pyaces summary in `run_tests.sh` reads every `STATUS: n` line of the runner's Summary block: PASS → PASS, MISMATCH and FAIL_CLEAN → FAIL, everything else (CRASH, NO_RESULT, TIMEOUT …) → OTHER.
  - The rules that had been hand-edited into the TSV moved into `classify()`: `SCF_TYPE=KS` → KS/Runintgrt, `FOCK=AO` MBPT(2) → Rundirmp2, and CCSD(TQ), `BRUECKNER` and `HFSTABILITY=FOLLOW` excluded. Added `classify_suite.py --check` (exit 1 if the TSV differs from the rules).
  - Regenerated the TSV. The 51 runnable rows are unchanged; only three description strings of excluded rows and case 020's iuhf (0 → 1, because of REF=UHF) changed.
  - The runner's default work directory became `build_logs/pyaces_suite_run` in the repository, and on HiPerGator `run_tests.sh` points it at `/blue`.
- **Verification:**
  - Compared on fabricated logs: MISMATCH 2, TIMEOUT 1, CRASH 1 → new summary `TOTAL=4 FAIL=2 OTHER=2` (the old code gave `?` throughout). CRASH 3 with no PASS → `TOTAL=3 OTHER=3` (previously `?`).
  - `classify_suite.py --check`: matches locally and on ccsd (125 rows).
  - For the local and ccsd suite results, see T-03.
- **Notes:** a mistyped command during testing overwrote the local `build_logs/unified_test_pyaces.log`; the pyaces tests were rerun to restore it (51/51).

### T-07. Python version no longer pinned (2026-10-23)
- **Commits:** `95cc01a2`
- **What changed:** `build.sh` names the extension module from `sysconfig.get_config_var('EXT_SUFFIX')` (the hard-coded `cpython-312` is gone). The audit tool already searched for `aces2py*.so`, so it needed no change.
- **Verification:** on 3.12 the name (`aces2py.cpython-312-x86_64-linux-gnu.so`) and the results are unchanged (local and ccsd).
- **Open:** not tested with any other Python version. Only 3.12 is available locally, and ccsd's python/3.10, 3.11 and 3.6 have no numpy.

### T-04. Build stops when a rename fails (2026-10-23)
- **Commits:** `95cc01a2`
- **What changed:** before patching, every pair in `EXTRA_LOCALIZE_SYMBOLS` (including the automatically added ones) is checked. If the module is not in `MAIN_COLLISION_MODULES`/`EXTRA_LOCALIZE_MODULES`, or the archive has no `T <symbol>`, the build prints the list and stops. The `2>/dev/null || true` after `objcopy --redefine-sym` is gone too.
- **Verification:**
  - All 60 current pairs are valid (pre-check).
  - A build with deliberately added `vcc:nosuch_` and `notamod:finish_`: both pairs reported, stopped with rc=1.

### T-03. Check for calls bound to COMMON blocks (2026-10-23)
- **Commits:** `95cc01a2`
- **What changed:** added `collision_audit.py --calls`. If a `U` symbol of a pyaces Fortran object (`Main/build/*.o`, excluding the two C objects) is data (B/b/D/d) in the `.so`, it is reported and the tool exits 1. It runs in build step 5 (right after the import test), so a hit fails the build.
- **Verification:**
  - Objects built from `runhbar.F`/`buildhbar.F` as they were before f97a9759: two `hbar_` hits reported, rc=1. Current objects: 0, rc=0.
  - Local: pyaces build OK, aces2 125/125, pyaces 51/51, audit 296 (matches the baseline), `--calls` 0.
  - ccsd (PBS 7713/7714): the same.

### D-2026-10-23b. pyaces documentation written (2026-10-23)
- **Commits:** `86a51cda`, `3a123417`, `daecf899`, `259a4fda`
- **What changed:**
  - Wrote the developer manual (`docs/developer_manual.html`, `docs/pyaces_developer_manual.pdf`, 25 pages, English).
  - Wrote the improvement plan in Korean and English (`docs/pyaces_improvements_ko.pdf`, `docs/pyaces_improvements_en.pdf`).
  - Turned the improvement plan into [TODO.md](TODO.md) and started this record.
- **Verification:** PDFs generated and rendering checked; confirmed no Korean text was left in the English version.

### D-2026-10-23a. Three defects found while writing the manual, fixed (2026-10-23)
- **Commits:** `f97a9759`
- **What changed:**
  - `Call Hbar` → `Call Hbar_` in `Runhbar`/`Buildhbar`. `hbar_` is not a subroutine but lambda's `COMMON /HBAR/` (`B hbar_` in the .so), so the branch taken after SCF and CC had finished could jump into BSS.
  - `Integer Iuhf` → `Integer*8 Iuhf` in `run*.F`, removing the mismatch where f2py passed 4 bytes and the `-i8` Fortran read 8.
  - Added `Cf2py intent(out)` to `Getreclenpy` and `Cgetrecpy` so that the length and the string are returned to Python. Removed the debug output from `Cgetrecpy`.
- **Verification:**
  - Local: pyaces 51/51, audit 296 (matches the baseline). Case 026 run as `Init → Runscf → Runhbar`; CCSD energy within 7×10⁻¹¹ of the reference.
  - ccsd (PBS 7700/7701): pyaces 51/51, audit 296.
- **API change:** `Length = Getreclenpy(Name)`, `Dest = Cgetrecpy(Name, Length)`.
- **Notes:** calling `Runhbar` after `Runcc` takes the SCF branch and calls Vtran twice. The manual says so.

### D-2026-10-14b. Link result changing from build to build, fixed (2026-10-14)
- **Commits:** `eab3c4e6`
- **What changed:**
  - Sorted the source lists in the aces2 makefile with `$(sort …)` (GNU Make 4.2 does not sort wildcard results; 4.3 sorts them in locale order).
  - Added `export LC_ALL=C` to `pyaces/build.sh`.
  - Added `getdt1` (vcc), `fndmaxs` (vee) and `t2l2rng`/`r1l2y2c` (lcct), whose winner depended on link order, to the rename list. The baseline entry where lcct used vee's `r1l2y2c` went away, taking the baseline from 297 to 296.
- **Verification:** local and ccsd agree: aces2 125/125, pyaces 51/51, audit 296.

### D-2026-10-14a. ccsd server installation and oneAPI 2023.2 support (2026-10-14)
- **Commits:** `73944a6d`, `dbed9377`
- **What changed:**
  - Installed on ccsd (`/home/ycpark/program/ACES`). Updates go through `git push ccsd master`; builds and tests run only through PBS.
  - Moved the `#include` lines in `pyaces/Main/init.F` to column 1 (ifort 2023.2 rejects indented includes).
  - Added `-assume dummy_aliases` to the aces2 Intel build (ifx 2023.2 miscompiled vmol's `VGAMMA2` at `-O2`, making the SCF energy 0.2 Eh wrong).
- **Verification:** on ccsd, aces2 125/125, aces3 37/37, pyaces 51/51.
