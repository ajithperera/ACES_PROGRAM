# ACES_PROGRAM

Unified repository for the ACES quantum chemistry program system:

- `aces2/` — ACES II, the serial legacy code (also the trusted correctness reference for ACES III)
- `aces3/` — ACES III, the parallel SIAL/SIP code
- `pyaces/` — pyaces (acespy), the Python/F2PY interface to ACES II

Each subdirectory's history was preserved via `git subtree` from its original standalone repository
(`ajithperera/ACES-II`, `ajithperera/acesiii`, `ajithperera/acespy`), which remain on GitHub as
historical archives.

## Building and testing

Each subsystem keeps its own build system (aces2: `Makefiles/GNUmakefile` via `gmake`; aces3:
autotools; pyaces: `pyaces/build.sh`, f2py) -- nothing here replaces those. `configure`/`build.sh`/
`run_tests.sh` at this top level are a thin orchestrator that records the shared toolchain (Intel
compilers + Intel MKL, via `module load intel/2025.1.0 openmpi/5.0.7 python/3.12`) and drives all
three consistently:

```
./configure                      # writes build_config.sh (module set, default -np)
./build.sh                       # builds aces2, then aces3, then pyaces (or pass a subset)
./run_tests.sh                   # runs all three fast suites, prints one rollup table
```

On HiPerGator, run these via sbatch rather than on an interactive shell:
`sbatch build_logs/submit_unified_build.sh`, `sbatch build_logs/submit_unified_tests.sh`.

**Baseline**: current numbers, aces2/pyaces measured 2026-09-04 and aces3 2026-09-09, on the local workstation
(`build_logs/unified_test_*.log`) -- treat these as the reference point, not something to re-derive
from scratch each session:

| suite  | result                | notes                                                          |
|--------|-----------------------|----------------------------------------------------------------|
| aces2  | 125 pass / 0 fail / 0 error | case 020 runs `REF=UHF`; the RHF CCSD(TQ) path is still broken |
| aces3  | 37 pass / 0 fail      | `ccsd_f12` needs `--jobs 1` or `--timeout-scale 3` locally -- at `--jobs 4` it exceeds its 1200s budget on shared cores (549s alone) |
| pyaces | 51 pass / 0 fail      | 51 runnable of 125 classified cases (59 EXCLUDED, 15 SCRIPT)   |

Git tag `zero-point-20260812` is the historical zero point (aces2 121 pass/1 fail/3 error, aces3
37/0, pyaces 37 pass/16 fail_clean/4 crash/1 no_result of 58 candidates); everything above is
measured at master HEAD and supersedes it.
