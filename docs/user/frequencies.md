# Vibrational frequencies

`VIB=EXACT` computes the force constants analytically, so far for SCF
wavefunctions (test `004a`, RHF/6-31G\*\*, reference SCF energy
−76.0236149977):

```{literalinclude} ../../aces2/tests/fast/cases/004a/ZMAT
:language: text
```

`VIB=FINDIF` builds the force constants by finite differences of analytic
gradients, or of energies when gradients are not available, in
symmetry-adapted mass-weighted Cartesian displacements. It works with any
method. Test `001f` (SCF) and test `047d` (CCSD) are examples. `FD_IRREPS`
limits the displacements to the listed irreps.

(isotopes)=
## Isotopes

Frequencies use the mass of the most abundant isotope of each element. To
use other masses, put a file named `ISOMASS` in the work directory before
the run, with one mass (in amu) per line for each real atom, in Z-matrix
order (dummy atoms are skipped). For D<sub>2</sub>O with the water input
below, `ISOMASS` is

```text
2.01410178
15.99491462
2.01410178
```

and for H<sub>2</sub><sup>18</sup>O it is `1.00782503`, `17.99915961`,
`1.00782503`.

```text
Water SCF frequencies by finite differences
H
O 1 R
H 2 R 1 A

R=0.947
A=105.5

*ACES2(CALC=SCF,BASIS=6-31G*,VIB=FINDIF,SCF_CONV=10)
```

The output lists the masses (`@VIB1-I, Non-standard isotopic masses
used.`) before the normal-mode analysis. Harmonic frequencies in cm⁻¹ from
these runs:

| Mode | H<sub>2</sub>O | D<sub>2</sub>O | H<sub>2</sub><sup>18</sup>O |
|---|---|---|---|
| A<sub>1</sub> bend | 1826.05 | 1336.22 | 1818.51 |
| A<sub>1</sub> stretch | 4064.59 | 2930.63 | 4055.93 |
| B<sub>2</sub> stretch | 4182.31 | 3065.74 | 4165.05 |

The force constants do not depend on the masses, so a finished frequency
job can be reused: copy its directory, add `ISOMASS`, and run `xjoda`
alone. It reads the force constants from `JOBARC` and repeats only the
normal-mode analysis. For D<sub>2</sub>O this gives 1336.21, 2930.63 and
3065.74 cm⁻¹, within 0.02 cm⁻¹ of the full run. The infrared intensities
from such a rerun are not right (B<sub>2</sub>: 32.0 instead of
36.6 km/mol); run the whole job again when you need intensities.
