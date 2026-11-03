# aces3 suite: 16 Intel MPI ranks on 4 cores

On a 4-core allocation with Intel MPI, `./run_tests.sh` runs the aces3
suite so oversubscribed that cases take many times longer than their
budget.

- **Status:** open (worked around by running aces3 by hand)
- **Area:** `run_tests.sh`, `aces3/tests/fast/run_tests.py`
- **Found with:** local workstation, ifort 2024.0 + Intel MPI 2021.11,
  4 cores (`taskset -c 20-23`), `ACES_PROGRAM_NP=4`

## Symptom

`run_tests.sh` calls the aces3 runner with `--jobs $ACES_PROGRAM_NP --np
$ACES_PROGRAM_NP`, which is 4 cases at once with 4 ranks each: 16 MPI
processes. On 4 cores each process got about 23 % of a core. After about
2 h 40 min the suite was at case 27 of 37; `dea_eomcc_f12` alone had run
1 h 50 min and was still in its second CCSD iteration (its budget is
2700 s × `--timeout-scale 3` = 8100 s).

Intel MPI waits for messages by spinning, so a rank that waits still uses
its whole time slice, and oversubscription slows every case down, not only
the waiting ones. The same `--jobs 4 --np 4` on this machine with Open MPI
finished 37/37 (with `--timeout-scale 3`), and the ccsd PBS job asks for
16 CPUs, so neither showed the problem.

## Workaround

Run aces3 one case at a time and the other suites as usual:

```
./run_tests.sh aces2 ring pyaces
cd aces3 && python3 tests/fast/run_tests.py --jobs 1 --np 4 --timeout-scale 3
```

With `--jobs 1` the suite passed 37/37 in about 3 h (10884 s of case
time) on the same 4 cores.

## Possible fixes

1. Let `run_tests.sh` take the aces3 job count separately (for example
   `ACES3_JOBS`, default `$ACES_PROGRAM_NP`), or choose it from the number
   of usable cores (`nproc`) divided by `--np`.
2. With Intel MPI, ask for less spinning when oversubscribed
   (`I_MPI_WAIT_MODE`, `I_MPI_SPIN_COUNT`). Not tried.

## Open questions

- After `subprocess.run(..., timeout=...)` kills `mpirun`, does the runner
  also stop the `xaces3` ranks? Not checked; the run above was stopped by
  hand before any case reached its timeout.
