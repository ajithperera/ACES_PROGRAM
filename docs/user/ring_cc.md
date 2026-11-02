# Ring coupled-cluster and RPA

`CALC=RCCD` and `CALC=DRCCD` give ring CCD and direct ring CCD (the latter is
equivalent to direct RPA for the correlation energy). `xaces2` runs them with
their own programs, `xringcc` for the amplitudes and, for excitation
energies, `xringcl` (Λ) and `xvee`. The RHF ring-CCD energy combines the
singlet and triplet channels with equal weights (Eq. 49 of the paper cited
below) and equals the UHF value for a closed shell. Ring-CCD for water
(`aces2/tests/ringcc`, reference −76.1821730567, the same with `REF=UHF`):

```{literalinclude} ../../aces2/tests/ringcc/cases/rccd_rhf/ZMAT
:language: text
```

Excitation energies follow with `EXCITE=RPA` or `EXCITE=DRPA` on an SCF
reference, or with the EOM ring methods (`EXCITE=EOM(S)-RCCD`,
`EOM(S)-DRCCD` and their spin-flip forms) together with `EOMREF=RCCD` or
`EOMREF=DRCCD`. `RPA_SOLVE=DAVIDSON` replaces full diagonalization with an
iterative solver. The inputs in `aces2/tests/jcp2020_eom_rccd` reproduce
the calculations of V. Rishi, A. Perera and R. J. Bartlett,
*J. Chem. Phys.* **153**, 234101 (2020).
