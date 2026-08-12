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

**Baseline**: git tag `zero-point-20260812` marks all three suites confirmed solid (aces2 121
pass/1 fail/3 error, aces3 37/0, pyaces 37 pass/16 fail_clean/4 crash/1 no_result) -- treat this as
the reference point, not a number to re-derive from scratch each session.
