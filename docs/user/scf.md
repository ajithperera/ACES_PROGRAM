# SCF

`REFERENCE` selects RHF (the default), UHF, ROHF, or a two-determinant
reference (`TWODET`) for open-shell singlets. Convergence is controlled by
`SCF_CONV` (10<sup>-N</sup> on the density, default 7), `SCF_MAXCYC`, and
the extrapolation keywords `SCF_EXTRAP`, `SCF_EXPSTAR` and
`SCF_EXPORDE`. `OCCUPATION` fixes the number of electrons in each irrep (see
{ref}`zmat-namelist`).

`HFSTABILITY=ON` checks whether the SCF solution is a minimum. `FOLLOW`
also moves to the lower solution along the unstable direction and reruns the
SCF there. Test `083` breaks the symmetry of a UHF wavefunction for
stretched HF this way, then runs CCSD(T) with the 1s orbital frozen
(reference total energy −100.153602773):

```{literalinclude} ../../aces2/tests/fast/cases/083/ZMAT
:language: text
```

Kohn–Sham DFT is selected with `SCF_TYPE=KS` and a `*INTGRT` namelist
that names the functional. LDA, GGA, global hybrid (B3LYP, QTP17) and
range-separated hybrid functionals (CAM-B3LYP, CAM-QTP-00/01/02, LC-QTP)
are available. {doc}`dft` describes the input and lists the test cases.
