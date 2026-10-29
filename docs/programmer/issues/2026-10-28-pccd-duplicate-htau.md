# pccd: duplicate `Pccd_form_htau` in one archive

Two aces2 pCCD cases stopped with a bogus "Insufficient memory" on the
local workstation because the link bound a call to the wrong one of two
routines with the same name.

- **Status:** fixed in `0a27e717`
- **Area:** aces2, `aces2/pccd/pccd_read_htau.F` (and `.f`)
- **Found with:** local workstation, ifx 2024.0, Open MPI; aces2 fast
  suite after the upstream merge of 2026-10-28

## Symptom

`pccd_oo_rfo_ccd_h2o` and `pccd_oo_rfo_lccd_h2o` failed with

```
@Pccd_form_uhtau: Insufficient memory
```

although the job had plenty of memory. A build of plain upstream failed
the same way, so the merge was not the cause. On HiPerGator both cases
passed.

## How to reproduce

On a build from before `0a27e717`:

```
cd aces2
python3 tests/fast/run_tests.py --case pccd_oo_rfo_ccd_h2o
```

## Cause

`pccd/pccd_read_htau.F` defined a second, unused `Pccd_form_htau` with
three arguments, next to the real one in `pccd/pccd_form_htau.F`. Both end
up in `libpccd.a`. GNU ld reads an archive once, front to back, and
resolves a symbol with the next member that defines it, so the copy a
binary gets depends on member order. The call in `pccd_form_htau_2d.F`
was bound to the three-argument copy, which read array data as `Maxcor`
and reported that as too little memory. HiPerGator's member order happened
to pick the right copy.

## Fix

The unused routine in `pccd_read_htau.F` was renamed to `Pccd_read_htau`
(its `Insmem` message too), with a comment saying why.

## Verification

aces2 fast suite 143/143 on the local workstation after the fix (157/157
with the cases upstream added since). `aces2/Scripts/dup_symbols.py`
lists the other symbols defined twice in one archive; see
{doc}`../duplicate_symbols`.
