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

### Optimizations from Cartesian input

With `COORD=CARTESIAN` there is no Z-matrix to star, and the whole
structure is optimized in one of two coordinate systems:

`GEOM_OPT=RIC` (or `FULL`)
: Redundant internal coordinates. `xjoda` finds the bonds from covalent
  radii (with hydrogen bonds and links between fragments added). It builds
  every bond, bond angle, and dihedral from them; a nearly linear bend
  (within 5° of 180°) gets two perpendicular bending coordinates. Every
  atom with exactly three bonds also gets an out-of-plane coordinate (an
  improper dihedral), so that the pyramidalization of a planar center
  (H<sub>2</sub>CO, amides, NH<sub>3</sub> near planarity) is described; an empty
  file named `NO_RIC_OOPS` in the work directory leaves them out. Steps are
  taken in these coordinates and transformed back to Cartesians by
  iteration (Pulay et al., *J. Chem. Phys.* **96**, 2856 (1992)). A
  `CONNECT` file in the work directory replaces the automatic bonds.

`GEOM_OPT=CART`
: Plain Cartesian coordinates, with translations and rotations projected
  out of the Hessian. Only `OPT_METHOD=RFA` (minima) and `EVFTS`
  (transition states) are allowed.

Water with redundant internals (test `120c`):

```{literalinclude} ../../aces2/tests/fast/cases/120c/ZMAT
:language: text
```

**Symmetry.** The energy and gradient calculations use the point group of
the starting structure, as for Z-matrix input. That group is found with
a loose test for degenerate moments of inertia, but the atoms must still
match under the symmetry operations to 10<sup>−6</sup> bohr. Give
symmetric structures to about seven decimals: NH<sub>3</sub> with
0.866025 for √3/2 is found as C<sub>3</sub> only. The optimizer itself
treats every coordinate as independent. After every step it projects the
new Cartesian geometry onto the totally symmetric structures of the full
point group, so the point group never drops. A C1 starting structure is
optimized with `SYMMETRY=NONE`. If a structure that started with symmetry
reaches a higher point group (a bent molecule that becomes linear, say),
the run stops with `Ascent in symmetry detected`. Restart from the printed
coordinates, or use `SYMMETRY=NONE`. Before 10/2026 these optimizations
always ran with `SYMMETRY=NONE`.

**Constraints.** `OPT_CONTROL=CONSTRAINED` keeps selected coordinates at
their starting values. The list goes in a `*CONSTRAINS` block after the
`*ACES2` namelist, ended by `*END`. Constrained optimizations always run
with `SYMMETRY=NONE` and `NOREORI=ON`, so the coordinates refer to the
input orientation.

- With `GEOM_OPT=RIC` give the number of bond constraints followed by one
  line `i j` for each bond, then the number of angle constraints and the
  lines `i j k` (the angle at atom `j`), then the number of dihedral
  constraints and the lines `i j k l`. Trailing sections can be left
  off: `*END` may follow the bond or the angle lines. A constraint must match a
  coordinate that `xjoda` generated: an angle is looked up with `j` as
  its vertex first, then as any angle of the same three atoms; a
  dihedral is looked up as `i j k l` or `l k j i` first. Constraints that
  do not match are reported with the list of generated coordinates, and
  the optimization continues without them. Atoms are numbered in input
  order. The constraints are applied by projecting the Hessian and
  gradient (Peng et al., *J. Comput. Chem.* **17**, 49 (1996)). The test
  input `aces2/test/zmat.130` freezes two bonds and two angles:

  ```text
  *CONSTRAINS
  2
  1 3
  5 4
  2
  1 2 3
  5 6 4
  *END
  ```

- With `GEOM_OPT=CART` give the number of frozen Cartesian components and
  then their positions in the coordinate list, where atom *n* has x, y, z
  at 3*n*−2, 3*n*−1, 3*n*. `aces2/test/zmat.131` freezes five components:

  ```text
  *CONSTRAINS
  5
  1 5 9 11 16
  *END
  ```

**Convergence and limits.** An optimization has converged when the norm
of the (projected) gradient is below 10<sup>−`CONVERGENCE`</sup> (the
value printed as `RMS force` is this norm) and the largest component of
the next step is below 10<sup>1−`CONVERGENCE`</sup> bohr or radian. The
step test catches soft coordinates (torsions, pyramidalization) where a
small gradient still leaves the structure far from the minimum; it is
skipped once the gradient is a hundred times below the threshold. It
was added in 10/2026 and costs a few cycles on large floppy molecules
(glycine 9 → 13) whose energies were up to 10<sup>−5</sup> hartree too
high before. The Hessian update uses the step actually taken, which
differs from the requested one by what the back transformation could
not reproduce.
`INIT_HESSIAN=SPECIAL` (the default) starts RIC and CART optimizations
from Lindh's model Hessian (Lindh et al., *Chem. Phys. Lett.* **241**,
423 (1995)). Its force constants depend on the distances and
periodic-table rows of the atoms involved. RIC uses it as a diagonal
with a floor of 10<sup>−3</sup> a.u.; CART uses the full Cartesian
form, summed over all atom pairs, triples and quadruples. Compared with
the guesses used before 10/2026 it saves about a quarter of the cycles
for RIC and about half for CART. The old guesses were a constant
diagonal (1.0 for bonds, 0.25 for angles, 0.1 for dihedrals) for RIC
and a unit matrix for CART. An empty file named `OLD_GUESS_HESS` in the
work directory brings them back. Z-matrix optimizations keep their
constant diagonal, which did as well as the model Hessian on the test
suite. The
coordinate set is made once, at the first geometry. If the bonding changes
during the run, `xjoda` prints a message and keeps the original set. If
the iterative back transformation does not converge, the iterate closest
to the requested step is used.

**Transition states from reactant and product.** With
`OPT_CONTROL=QST` or `LST` and `OPT_METHOD=EVFTS`, a transition state
search with `GEOM_OPT=RIC` is guided by the synchronous transit path
between a reactant and a product (Peng and Schlegel, *Isr. J. Chem.*
**33**, 449 (1993)). The input holds three structures with the atoms in
the same order: the guess for the transition state (the main geometry),
then the reactant and the product, each after a one-line label and a
blank line:

```text
H2CO -> H2 + CO
C   0.000000   0.000000   0.000000
O   0.000000   0.000000   1.160000
H   0.376000   0.000000  -1.034000
H   1.503000   0.000000  -0.547000

Reactant
C   0.000000   0.000000   0.000000
O   0.000000   0.000000   1.200000
H   0.940000   0.000000  -0.540000
H  -0.940000   0.000000  -0.540000

Product
C   0.000000   0.000000   0.000000
O   0.000000   0.000000   1.130000
H   0.370000   0.000000  -2.800000
H  -0.370000   0.000000  -2.800000

*ACES2(CALC=SCF,BASIS=3-21G,GEOM_OPT=RIC,OPT_CONTROL=QST
OPT_METHOD=EVFTS)
```

The coordinates are built from the bonds of all three structures, so the
bonds being made and broken are included (Peng et al., 1996). The
tangent to the path at the current point is a straight line from
reactant to product for `LST`, and the circle through reactant, current
point and product for `QST`. The first two steps climb along this
tangent to the maximum of the path, and steps 3 and 4 do the same if
the estimated climb is longer than 0.05 au. Every later step is an
eigenvector-following step uphill along the Hessian eigenvector that
overlaps the tangent by more than 0.8, or else along the lowest one.
The starting Hessian is the model Hessian; its updates build up the
negative eigenvalue. For this input QST converges in 11 cycles and
LST in 16, to the transition state that an exact-Hessian search finds
from the same guess. Plain `EVFTS` from the guess wanders off. Before
10/2026 this path was broken for Cartesian input.

Without a guess for the transition state (QST2), give only two
structures: the main geometry is the reactant, and the product follows
after a label line and a blank line. `xjoda` then starts from the
midpoint of the linear synchronous transit path (Halgren and Lipscomb,
*Chem. Phys. Lett.* **49**, 225 (1977)): the product is rotated onto the
reactant, the interatomic distances are interpolated halfway, and the
Cartesians are fitted to them by least squares. The search continues as
above. For H2CO → trans-HCOH, QST2 takes 8 cycles (LST 9) to the
same transition state as the three-structure input. For vinyl alcohol →
acetaldehyde it finds the nonplanar 1,3-hydrogen shift transition state
in 22 cycles (one imaginary frequency, 2513i cm⁻¹). The planar
three-structure guess for that reaction leads to a different saddle
point, where H stays on O and the mode is out of plane (1491i cm⁻¹).

The midpoint has every symmetry that the reactant and product share
with the atoms in the given order, and the search keeps that symmetry.
When the transition state is less symmetric, QST2 cannot reach it. An
example is H2CO → H2 + CO, where the midpoint is C2v but the transition
state is planar Cs: QST2 ends at a C2v point with three imaginary
frequencies. Use three structures in such cases.

For checking the coordinate code, an empty file named `RIC_FDCHECK` in the
work directory makes `xjoda` compare every B-matrix element with central
differences of the coordinate values and print the largest error for each
kind of coordinate. At an exactly linear bend the angle has a kink, so a
large error there is expected.

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
