# ACES III fast regression suite

Small (2-4 core), automated pass/fail suite covering ACES III's major
capabilities, replacing one-at-a-time-by-hand comparison against the old
PBS-batch `tests/` archive. Run with:

```
module load intel/2025.1.0 openmpi/5.0.7
python3 tests/fast/run_tests.py
```

Each case in `cases/<name>/` supplies a `ZMAT`; `run_tests.py` runs `xaces3`
in a scratch dir and compares the resulting JOBARC record (`TOTENERG` or
`SCFENEG`, whichever `test_results` lists for that case) against a locked-in
reference value via `xtest_compare`, gated by a tolerance (default `1e-6`
a.u. RMS, see `TOLERANCE_OVERRIDES`).

## g-tensor cases (`gten_*`)

`gten_no`, `gten_cn`, `gten_h2op` exist to exercise the CCSD g-tensor
property calculation itself, not just SCF convergence -- `SCFENEG` alone
doesn't touch the actual g-tensor tables. Each carries a
`reference_gtensor.txt` (an extracted, known-good copy of the g-tensor
tables `print_rel_info.F` writes to `summary.out`, unit 66). `run_tests.py`
parses both the reference and the fresh run's `summary.out` and diffs every
tensor element (default tolerance `1e-3` ppm RMS across all tables, see
`GTENSOR_TOLERANCE`/`GTENSOR_TOLERANCE_OVERRIDES`).

**Known caveat, not a bug:** `gten_no`'s paramagnetic g-tensor component is
not invariant, reflecting the orbitally-degenerate <sup>2</sup>&Pi; ground
state under single-reference UHF -- not a code defect. NO's true ground
state is genuinely multi-determinantal (a symmetry-adapted combination of
real p<sub>x</sub>/p<sub>y</sub> orbitals with real orbital angular
momentum), which single-reference UHF cannot represent; forcing SCF to
always land the unpaired electron in one fixed real orbital would make the
number reproducible but not physically correct. A real fix would require
multi-reference/multi-determinantal treatment of degenerate states -- new
method development, not test-suite hygiene. `gten_no` is therefore checked
against `SCFENEG` only (not the g-tensor tables) until/unless that method
work happens.

## F12 cases (`ccsd_f12`, `ccsdt_f12`, and any `*_f12` property case)

`STGF12=ON` deliberately runs the primary/orbital SCF over a fully
decontracted basis (every primitive treated as its own basis function --
confirmed directly by diffing the generated `MOL` file for water/aug-cc-pVDZ:
41 contracted basis functions with `STGF12=OFF` vs. 57 with it on, matching
each shell's primitive count 1:1). This is expected/correct behavior for
ACES III's STG-F12 method, not a bug -- so SCF and CCSD energies for F12
cases will *not* match the pre-2026 `/blue/bartlett/perera/TESTS` archive,
which predates (or otherwise doesn't reflect) this decontraction. **Do not
use that archive as ground truth for any `STGF12=ON` case.** `test_results`
entries for F12 cases are self-consistent baselines from a verified clean
run of the current code instead (same convention as `steom_so`), not
archive-matched. Note the decontraction is F12-specific plumbing, not the
general-purpose `UNCONTRACT` ASV keyword -- tried in isolation, `UNCONTRACT=ON`
does not reproduce it (confirmed 2026-08-03).
