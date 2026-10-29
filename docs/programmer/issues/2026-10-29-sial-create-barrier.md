# Missing barrier after `create` in the spin-orbit EOM programs

The ACES III case `so_ipeomcc_f12` failed now and then under Intel MPI
with `Block not found`. The cause was a race in the SIAL program
`eom_so_ci_ip.sial`: it creates distributed arrays and writes to them with
`PUT` without a `sip_barrier` in between. The same pattern is in
`eom_so_ci_ea.sial` and `eom_so_ci.sial`. This page records the symptom,
the cause, the fix and how it was checked, found on 2026-10-29.

- **Status:** fixed (commit pending)
- **Area:** aces3, SIAL programs in `aces3/src/aces/aces_sial/spnorbt_cc_new/`
- **Found with:** ccsd server, ifx 2023.2 + Intel MPI 2021.10; also listed in upstream's tag `baseline-20261028`

## Symptom

On the ccsd server (ifx 2023.2, Intel MPI 2021.10, `--np 4`) the case
stopped after about 19 s in the last SIAL program of its jobflow,
`eom_so_ci_ip.sio`:

```
Task 1, array 106, block 1   Block not found
Message type -6
 Task            2  Current instruction counter is          562
  current line number is          906
Abort(1) on node 2 (rank 2 in comm 0): application called MPI_Abort(MPI_COMM_WORLD, 1) - process 2
```

The rank that reports the missing block and the rank that is still at
line 906–907 change from run to run. Before the fix it failed 5 times in
13 runs (2 of 3, then 1 of 5 with the current ZMAT and 2 of 5 with the
older `scf_aguess` ZMAT, so the 2026-10 SCF change in the ZMAT is not
the cause). On the local workstation (ifx 2024.0, Open MPI) it was not
seen to fail. Upstream's tag `baseline-20261028` lists it as known:
"aces3 so_ipeomcc_f12 intermittent under Intel MPI".

## Cause

`PROC SO_INT` in `eom_so_ci_ip.sial` starts like this:

```
      create SO_x_ao
      ...
      create SO_z_pq

      PARDO mu,p
      tSO_x_pmu(mu,p) = 0.0
      PUT SO_x_pmu(mu,p) = tSO_x_pmu(mu,p)
      ENDPARDO mu,p
```

`create` allocates, on each worker, only the blocks of the distributed
array that this worker owns. It does not synchronise the workers. `PUT
X(...) = ...` sends the block to its owner as a "copy data" message
(message type −6, `msg_copy_data` in `aces3/src/main/sumz.c`). The owner
looks the block up and aborts if it is not there:

```c
if (info[index].block == -1)
{
   printf("Task %d, array %d, block %d   Block not found\n", ...);
   printf("Message type %d\n", info[index].msgtype);
   F77_NAME(abort_job, ABORT_JOB)();
}
```

So when one worker reaches the `PARDO` and sends its first `PUT` before
the owning worker has run the `create` lines, the owner has no such block
yet and the job stops. The worker that is still at line 906–907 in the
log is the one that has not finished its `create` lines. Whether this
happens depends on timing, which is why Intel MPI showed it and Open MPI
on the local workstation did not.

The SIAL manual requires the barrier (see `sip_barrier` in {doc}`../sial`):
"Must be used after distributed arrays are created, before distributed
arrays are deleted". In the other SIAL programs a `create` is followed by
a barrier in 3333 places and not in 996.

## Fix

An `execute sip_barrier` after each group of `create` lines that is
followed directly by a `PARDO`. All nine are in `PROC SO_INT`:

| File | New line | Arrays created just before |
|---|---|---|
| `eom_so_ci_ip.sial` | 913 | `SO_{x,y,z}_ao`, `SO_{x,y,z}_pmu`, `SO_{x,y,z}_pq` |
| `eom_so_ci_ip.sial` | 954 | `Fso_pq_{x,y,z}_aux` |
| `eom_so_ci_ip.sial` | 1516 | `Fso_pq_{x,y,z}`, `somfex_{x,y,z}_{1,2,3}` |
| `eom_so_ci_ea.sial` | 878, 919, 1476 | the same three groups |
| `eom_so_ci.sial` | 1666, 1707, 2277 | the same three groups |

All files are in `aces3/src/aces/aces_sial/spnorbt_cc_new/`. `gmake -C
aces3/src/aces/aces_sial` recompiles the changed `.sial` files and copies
the `.sio` files to `aces3/bin/sio`; the full `./build.sh` does the same.

The fix adds one synchronisation per group; the cost is negligible next to
the integral work in `SO_INT`.

## Verification

| Where | Case | Before | After |
|---|---|---|---|
| ccsd (ifx 2023.2, Intel MPI) | `so_ipeomcc_f12` (uses `eom_so_ci_ip`) | 5 of 13 failed | 10 of 10 passed, RMS error 4.2e-10 each time |
| ccsd | `eom_ee` (uses `eom_so_ci`) | passed | 2 of 2 passed |
| local workstation (ifx 2024.0, Open MPI) | both cases | passed | passed |

With a failure rate of 5 in 13, ten passes in a row by chance has a
probability of about 1 %. No test case runs `eom_so_ci_ea.sio`; that file
got the same change by analogy and was only checked to compile.

To repeat the check (on ccsd through PBS, on the workstation on the
allocated cores):

```
cd aces3
for i in $(seq 10); do
  python3 tests/fast/run_tests.py --np 4 --case so_ipeomcc_f12 | grep -E '^so_ipeomcc|Block not found'
done
```

## Not covered

Other `create` lines without a barrier remain in these programs and in
the rest of `aces3/src/aces/aces_sial`. Most are followed by
`list_to_blocks`, a `DO` loop or other local work, which then ends in a
barrier before the array is read, and none has been seen to fail. This
lists, for one directory, each `create` group that is not followed by a
barrier and the line that comes next:

```
cd aces3/src/aces/aces_sial/spnorbt_cc_new
for f in *.sial; do
  awk -v F=$f '/^ *[Cc][Rr][Ee][Aa][Tt][Ee] /{c=1; s=NR; next}
    c && NF && !/^ *#/ { if ($0 !~ /barrier/) print F ":" s ": " $0; c=0 }' $f
done
```

A `create` group followed directly by `PARDO ... PUT` or `PREPARE` is the
pattern to fix.
