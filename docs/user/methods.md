# Methods and calculation types

A handful of keywords decide what an ACES II job does:

| What | Keyword | Examples |
|---|---|---|
| Reference wavefunction | `REFERENCE` | `RHF`, `UHF`, `ROHF`, `TWODET` |
| Correlation method | `CALCLEVEL` | `MBPT(2)`, `CCSD`, `CCSD(T)`, `PCCD` |
| Excited, ionized, attached states | `EXCITE`, `IP_CALC`, `EA_CALC`, `DIP_CALC`, `DEA_CALC` | `EXCITE=EOMEE` |
| Properties | `PROPS`, `ESTATE_PROP` | `PROPS=FIRST_ORDER` |
| Geometry | `*` on Z-matrix coordinates, `GEOM_OPT`, `OPT_METHOD` | `OPT_METHOD=EVFTS` |
| Frequencies | `VIBRATION` | `VIB=EXACT`, `VIB=FINDIF` |

Every example below is an input from the ACES II test suite, included
directly from `aces2/tests/`, with the reference energy the suite checks
against. They are small enough to run in seconds to minutes.

## SCF

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

Kohn–Sham DFT is selected with `SCF_TYPE=KS` and the `*VSCF` and
`*INTGRT` namelists. The test cases `103.ks`, `104.ks` and `105.ks` show
the input.

## Correlated energies

`CALCLEVEL` (shortened to `CALC`) selects the method. The main families:

Many-body perturbation theory
: `MBPT(2)`, `MBPT(3)`, `SDQ-MBPT(4)`, `MBPT(4)`; `OO-MP2` with orbital optimization

Coupled cluster
: `CCD`, `CCSD`, `CCSD(T)`, `CCSD[T]`, `ACCSD`, `ACCSD(T)`, `CC2`, `CC3`,
  the iterative triples approximations `CCSDT-1` to `CCSDT-4`, `CCSDT`,
  quadruples corrections and approximations up to `CCSDTQ`, and the linear
  methods `LCCD`, `LCCSD`, `LCCSDT`

Configuration interaction and QCI
: `CID`, `CISD`, `QCISD`, `QCISD(T)`, `QCISD(TQ)`

Pair coupled cluster and ring methods
: `PCCD` and its variants, `OO-CCD`, `OO-LCCD`; `rCCD`, `drCCD` (see below)

The {ref}`keyword-reference` lists every allowed value.

A CCSD energy for water (test `109`, reference −76.2697698466):

```{literalinclude} ../../aces2/tests/fast/cases/109/ZMAT
:language: text
```

`INCORE` keeps integral lists in memory to cut disk I/O, and
`ABCDTYPE=AOBASIS` computes the four-virtual-index contribution from AO
integrals, which avoids storing the largest MO integral list.

### Frozen orbitals and truncated virtual spaces

`DROPMO` removes orbitals from the correlation treatment, counted from the
lowest. `DROPMO=1` freezes the oxygen 1s, and `DROPMO=1-2/10` freezes
orbitals 1, 2 and 10. `TRUNC_ORBS=FNO` replaces the virtual space with
frozen natural orbitals, keeping the percentage given by `FNO_KEEP`
(test `124a.fno`, reference −76.2496011642):

```{literalinclude} ../../aces2/tests/fast/cases/124a.fno/ZMAT
:language: text
```

### Brueckner orbitals

`BRUECKNER=ON` repeats the CC calculation, rotating the orbitals until the
T1 amplitudes vanish to within `BRUCK_CONV`. Test `022` runs CCSD with
`BRUECKNER=ON`, which gives Brueckner CCD (B-CCD).

## Excited, ionized and electron-attached states

Excited states of the same charge come from `EXCITE`. `EOMEE` is EOM-CCSD on
a CCSD ground state; `TDA` (CIS), `RPA`, `CIS(D)`, `P-EOMEE`, `STEOM`, and
EOM methods built on ring-CCD ground states are also available. The number
of states in each irrep is given with `ESTATE_SYM`, one number per irrep,
separated by `/`. `ESTATE_PROP=EXPECTATION` also computes transition moments,
oscillator strengths, and one-electron properties of the excited states.

EOM-CCSD for water, one state in each of the four irreps (test `054a`,
reference ground-state energy −76.2725089278):

```{literalinclude} ../../aces2/tests/fast/cases/054a/ZMAT
:language: text
```

`IP_CALC` and `EA_CALC` compute ionization potentials and electron
affinities of a closed-shell parent state, for example with `IP_EOMCC` and
`EA_EOMCC`. `IP_SYM` and `EA_SYM` give the number of states in each irrep.
`DIP_CALC` and `DEA_CALC` give double ionization and double attachment
(TDA, EOM-CC or STEOM). The electron affinity of the sodium atom with
EA-EOM-CCSD on a UHF reference (test `042`, reference −161.973978215):

```{literalinclude} ../../aces2/tests/fast/cases/042/ZMAT
:language: text
```

## Pair coupled-cluster doubles (pCCD)

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

## Ring coupled-cluster and RPA

`CALC=RCCD` and `CALC=DRCCD` give ring CCD and direct ring CCD (the latter is
equivalent to direct RPA for the correlation energy). Ring-CCD for water
(`aces2/tests/ringcc`, reference −76.3039658461):

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

## Properties

`PROPS=FIRST_ORDER` (or `PROPS=1`) computes the dipole and quadrupole
moments, electric field gradients, spin densities and the Cowan–Griffin
relativistic correction from the relaxed density of the chosen method
(test `026`, CCSD, reference −76.2725454575):

```{literalinclude} ../../aces2/tests/fast/cases/026/ZMAT
:language: text
```

Other values of `PROPS`:

- `SECOND_ORDER`: static polarizabilities (CPHF, with `CALC=SCF`)
- `NMR`: NMR chemical shifts with GIAOs (SCF and MBPT(2), with
  `SPHERICAL=OFF`), as in test `037`
- `J_FC`, `J_SD`, `J_SO`, `JSC_ALL`: contributions to NMR spin–spin
  coupling constants (`CALC=CCSD`; the Fermi-contact and spin-orbit terms
  need `REF=UHF`; test `046`)
- `EOM_NLO`: static and frequency-dependent EOM-CCSD polarizabilities
  (tests `s_075b`, `s_075c`)

## Geometry optimization

Put `*` after the names of the Z-matrix coordinates to optimize. `xaces2`
then repeats energy and gradient steps until the geometry converges.
Gradients are analytic whenever the method has them, and otherwise come
from finite differences of energies (`GRAD_CALC=NUMERICAL` forces this).
A CCSD optimization of water (test `013a`, reference −76.2725483266):

```{literalinclude} ../../aces2/tests/fast/cases/013a/ZMAT
:language: text
```

Related keywords:

- `GEOM_OPT=FULL` optimizes every coordinate, starred or not; `RIC` uses
  redundant internal coordinates, which also works for Cartesian input
  (test `120a`).
- `OPT_METHOD` chooses the step algorithm: Newton–Raphson (`NR`), rational
  function (`RFA`), or eigenvector following for transition states
  (`EVFTS`).
- `CONVERGENCE`, `OPT_MAXCYC` and `MAX_STEP` control convergence.
- A starting Hessian can be read from `FCMINT`. Test `126c` searches for the
  HCN ⇌ HNC transition state with `OPT_METHOD=EVFTS`, starting from a
  Hessian placed with a `%FCMINT=` file directive (reference −93.1153442031):

```{literalinclude} ../../aces2/tests/fast/cases/126c/ZMAT
:language: text
```

`RESTART` (on by default) lets an interrupted optimization or
finite-difference job continue from its last saved step; see
[](#file-directives) for `SAVEDIR`.

## Vibrational frequencies

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
limits the displacements to the listed irreps, and an `ISOMASS` file
replaces the atomic masses.
