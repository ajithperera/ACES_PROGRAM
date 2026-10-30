# Issues: bug reports and fixes

One page per failure: what was seen, how to reproduce it, the cause (or
what is known of it), the fix and how it was checked. Open reports are
for whoever picks the problem up next; fixed ones record why the code
looks the way it does.

File names are `YYYY-MM-DD-short-name.md`, with the date the problem was
found. Start a page with a one-line summary and this header:

```
- **Status:** open | fixed (commit) | known, not fixed
- **Area:** subsystem and files
- **Found with:** machine, compiler, MPI, test case
```

then the sections *Symptom*, *How to reproduce*, *Cause*, *Fix*,
*Verification* and, where something is left, *Not covered* or *Open
questions*. Add the page to the table and to the toctree below, and update
the status when it changes.

| Found | Report | Area | Status |
|---|---|---|---|
| 2026-10-28 | [pccd: duplicate `Pccd_form_htau` in one archive](2026-10-28-pccd-duplicate-htau.md) | aces2 `pccd` | fixed `0a27e717` |
| 2026-10-29 | [Missing barrier after `create` in the spin-orbit EOM programs](2026-10-29-sial-create-barrier.md) | aces3 SIAL | fixed in `76353e87` |
| 2026-10-29 | [EOM ring-CC: Ne `Ag` doubles correction differs between compilers](2026-10-29-ring-ne-ag-correction.md) | aces2 `vee` | open |
| 2026-10-29 | [aces3 suite: 16 Intel MPI ranks on 4 cores](2026-10-29-aces3-suite-oversubscription.md) | `run_tests.sh` | open |

For the symbols that are defined twice in one aces2 archive (the cause of
the pccd failure) see also {doc}`../duplicate_symbols`.

```{toctree}
:maxdepth: 1

2026-10-28-pccd-duplicate-htau
2026-10-29-sial-create-barrier
2026-10-29-ring-ne-ag-correction
2026-10-29-aces3-suite-oversubscription
```
