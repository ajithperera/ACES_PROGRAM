# ACES_PROGRAM

Unified repository for the ACES quantum chemistry program system:

- `aces2/` — ACES II, the serial code (also the correctness reference for ACES III)
- `aces3/` — ACES III, the parallel SIAL/SIP code
- `pyaces/` — pyaces, the Python (F2PY) interface to ACES II
- `docs/` — the manual (user guide, keyword reference, programmer manual)

Each subdirectory's history was preserved via `git subtree` from its original standalone repository
(`ajithperera/ACES-II`, `ajithperera/acesiii`, `ajithperera/acespy`), which remain on GitHub as
historical archives.

## Building and testing

```
cmake -S . -B build              # -DACES_MODULES="" to use the current environment
cmake --build build              # builds aces2, aces3, pyaces
ctest --test-dir build           # runs the three fast suites and a rollup table
```

Requirements: Intel oneAPI compilers and MKL, an MPI library, Python 3 with NumPy, CMake 3.18+.
On HiPerGator, build and test through Slurm, and run the ACES III suite on a single node
(`--nodes=1`). See `docs/getting_started/index.md` for options and details.

**Baseline**: current numbers, measured 2026-10-28 after merging upstream/master, on the local workstation
(`build_logs/unified_test_*.log`) -- treat these as the reference point, not something to re-derive
from scratch each session:

| suite  | result                | notes                                                          |
|--------|-----------------------|----------------------------------------------------------------|
| aces2  | 141 pass / 0 fail / 2 error | of 143; the errors are `pccd_oo_rfo_{ccd,lccd}_h2o` (`@Pccd_form_uhtau: Insufficient memory`), which fail the same way in a pure upstream build with ifx 2024.0 |
| aces3  | 37 pass / 0 fail      | `ccsd_f12` needs `--jobs 1` or `--timeout-scale 3` locally -- at `--jobs 4` it exceeds its 1200s budget on shared cores (549s alone) |
| pyaces | 51 pass / 0 fail      | 51 runnable of 125 classified cases (59 EXCLUDED, 15 SCRIPT)   |

Git tag `zero-point-20260812` is the historical zero point (aces2 121 pass/1 fail/3 error, aces3
37/0, pyaces 37 pass/16 fail_clean/4 crash/1 no_result of 58 candidates); everything above is
measured at master HEAD and supersedes it.

## The manual

```
python3 -m pip install --user -r docs/requirements.txt
docs/build.sh                    # HTML in docs/_build/html/index.html
```

The keyword reference, the JOBARC record catalog and the ACES III jobflow table are generated from
the source code each time the manual is built.
