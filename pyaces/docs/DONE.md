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
