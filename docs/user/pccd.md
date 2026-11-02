# Pair coupled-cluster doubles (pCCD)

`CALC=PCCD` runs pair coupled-cluster doubles (pCCD). `PCCD_REF` sets the
orbitals it uses: `HF` (the default), `NATURAL`, `LOCAL`, or `OPTIMIZED` for
orbital-optimized pCCD (OO-pCCD). `OO_METHOD` chooses the orbital
optimizer (`RFO`, `BFGS`, `FULL-NR`, ...), and `OO_MAXCYC` sets the maximum
number of orbital steps. `CALC=OO-CCD` and `CALC=OO-LCCD` optimize the
orbitals for CCD and linearized CCD in the same way. With `PROGRAM=PCCD`,
`xaces2` runs the pCCD programs (`xrunpccd`) after `xjoda`, instead of the
standard sequence.

OO-pCCD for water in STO-3G (test `pccd_oo_rfo_xaces2_h2o`, reference
−75.0056655037):

```{literalinclude} ../../aces2/tests/fast/cases/pccd_oo_rfo_xaces2_h2o/ZMAT
:language: text
```

A doubly occupied configuration interaction (DOCI) energy at the OO-pCCD
orbitals is computed by running `xrunpccd` and then `xoorot` with control
files. The script in test `pccd_oo_rfo_doci_h2o` shows the steps (reference
DOCI energy −75.005667504).
