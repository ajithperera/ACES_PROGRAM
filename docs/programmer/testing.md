# Test suites

Each part has a fast regression suite; `ctest --test-dir build` runs them
all (see {doc}`../getting_started/index`).

## ACES II

`aces2/tests/fast/run_tests.py` runs every case in
`aces2/tests/fast/cases/` in its own temporary directory, reads one
`JOBARC` record with `xa2proc jareq d <RECORD> 1`, and compares it with the
reference in `aces2/tests/fast/test_results`:

```text
013a zmat.013a
TOTENERG 1
-7.627254832660E+01
```

(case name, record, reference value). The default tolerance is 10<sup>-6</sup>
hartree; `TOLERANCE_OVERRIDES` in the runner can loosen it for one case.
The result is PASS, FAIL (wrong value), or ERROR (the job did not finish or
the record is missing).

```bash
module load intel/2025.1.0 openmpi/5.0.7 python/3.12
cd aces2
python3 tests/fast/run_tests.py --jobs 4            # all cases
python3 tests/fast/run_tests.py --case 013a         # one case
```

`xaces2` must be on `PATH`.

A case is a directory with either:

- a `ZMAT`: run with `xaces2`, with `basis/GENBAS` copied in; or
- a `SCRIPT`: a shell script for jobs with several steps (for example an
  optimization followed by frequencies, or running member programs by
  hand). The script gets `$GENBAS` and `$ECPDATA` (the files in
  `basis/`) and `$USERDIR`, a directory outside its working directory where
  it can keep files between steps.

Every suite, ACES II, ACES III and pyaces, reads the basis sets from the one
`GENBAS` in `ACES_PROGRAM/basis/`. Add a basis set there, never to a copy in
a test directory.

**Adding a case**: create `cases/<name>/ZMAT` (or `SCRIPT`), run it with a
build you trust, check the result independently, and add the three lines
to `test_results`. Keep cases small enough to finish in seconds. A case
whose job is expected to exit with an error after writing its result goes
in `EXPECT_NONZERO_EXIT`.

The `erd_grad_*` cases check `INTEGRALS=ERD` gradients against the default
route. Their `SCRIPT` writes a one-cycle optimization (`OPT_MAXCYC=1`)
to `ZMAT` and sources `tests/fast/erdgrad_check.sh`. That runs the job
twice, once with `xvmol` + `xvdint` and once with `xerdints` +
`xerdgrad`. It fails if the two `GRADIENT` records differ by more than
10<sup>-8</sup> hartree/bohr. Then the runner checks the ERD energy as usual.

Other ACES II suites with their own runners: `aces2/tests/ringcc`,
`ringcl`, `rpa`, `jcp2020_eom_rccd` (published EOM ring-CC results), and
`pccd_slow` (larger pCCD jobs, run through Slurm).

`aces2/tests/legacy_archive/` (not tracked by git) holds the original ACES
II test collection: the `zmat.NNN` and `script.NNN` inputs, their old
outputs, and the `KEEP/` files. The `fast` cases with numeric names
(`001a`, `013a`, `s_034`, ...) were imported from it.

## ACES III

`aces3/tests/fast/run_tests.py` runs each case with `xaces3` under MPI
(`--np` ranks each, `--jobs` cases at once) and compares the energies with
the references in `aces3/tests/fast/test_results` using the
`xtest_compare` tool, with an RMS tolerance of 10<sup>-6</sup>. Run it inside a **single-node** job; the
cases keep their scratch files in node-local directories, and the runner
stops if Slurm gave it more than one node.

## pyaces

`pyaces/pyaces_suite_runner.py` runs each case's `ZMAT` through the
`aces2py` module and compares the result with the reference value in its
case table. Besides PASS it reports MISMATCH (finished, wrong value),
FAIL_CLEAN (stopped with a normal ACES error exit), CRASH, NO_RESULT, and
TIMEOUT.

## Rollup

`cmake/rollup.sh` reads the three suite logs and prints one table of
pass/fail/error counts. Compare it with the counts from the last run
before your change (commit messages that change the counts record them);
any change in the counts is the signal to look closer.
