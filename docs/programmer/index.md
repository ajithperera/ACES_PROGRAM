# Programmer manual

```{note}
This part is being written.
```

Planned chapters:

- **ACES II architecture**: the `xaces2` driver and member executables;
  the `JOBARC` record store (`GETREC`/`PUTREC`) with a catalog of records
  checked against the code; the list-based integral and amplitude files
  (`GETLST`/`PUTLST`, `MOINTS` and related files); memory management.
  Based on the ACES II programmer manual and `aces2/Info`.
- **Building**: the CMake superbuild, the ACES II makefiles and host table,
  the ACES III autotools build, and how pyaces links ACES II.
- **Adding to ACES II**: a new keyword (`joda/asv_ctl.c` and `include/flags.h`),
  a new member executable, and a new test case.
- **Coding conventions and traps**: fixed-form column 72, per-file include
  guards, implicit typing, and the `_DEBUG_LVL0` debug prints.
- **ACES III**: SIAL programming and the SIP runtime. Based on
  `aces3/doc/SIALProgrammerGuide.tex` and `SIPDesign.pdf`.
- **pyaces internals**: how the F2PY module is built and why some ACES II
  libraries are left out of its link.
- **Test suites**: layout, reference values, and how to add a case.
