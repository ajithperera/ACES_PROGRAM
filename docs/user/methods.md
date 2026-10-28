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

`INTEGRALS=ERD` computes the integrals with `xerdints` and the gradients
with `xerdgrad` (Flocke's OED/ERD) in place of `xvmol` and `xvdint`. It
covers SCF, MBPT and CC gradients for RHF, UHF, ROHF and QRHF references,
with ECPs and Cartesian or spherical basis sets. Test `erd_opt_h2o_ccsd` is `013a` with
this one keyword added and ends at the same energy, −76.2725483266:

```{literalinclude} ../../aces2/tests/fast/cases/erd_opt_h2o_ccsd/ZMAT
:language: text
```

The `erd_grad_*` tests compare single gradients from the two routes
(UHF MBPT(2), ROHF CCSD, QRHF CCSD, spherical frozen-core MBPT(2) and an
SBKJC ECP); they agree to 10<sup>−12</sup> hartree/bohr. Analytic second
derivatives (`VIB=EXACT`), NMR and other properties still use `xvdint`.

The options below were checked in 10/2026 on SCF/3-21G minima (HOOH as a
Z-matrix, ethanol from Cartesians) and transition states (HCN ⇌ HNC as a
Z-matrix, H<sub>2</sub>CO → trans-HCOH from Cartesians). The cycle
counts quoted are from those runs.

**Coordinates.** Which coordinates are optimized depends on the input:

| Input | Keyword | Coordinates |
|---|---|---|
| Z-matrix with `*` | (none; `GEOM_OPT=PARTIAL` is set for you) | the starred Z-matrix variables |
| Z-matrix or Cartesian | `GEOM_OPT=RIC` (or `FULL`, the same) | redundant internal coordinates of the whole molecule (test `120a`) |
| Cartesian (`COORD=CARTESIAN`) | `GEOM_OPT=CART` | Cartesian coordinates |

RIC is the most robust choice for floppy molecules: the HOOH minimum
takes 43 cycles in Z-matrix coordinates and 5 for ethanol in RIC.
[](#optimizations-from-cartesian-input) has the details of RIC and CART.

**Step algorithm (`OPT_METHOD`).**

| Value | Search | Notes |
|---|---|---|
| `AUTO` (default) | minimum | becomes `MANR` |
| `NR` | minimum | Newton–Raphson; no protection against a wrong-signed Hessian (failed for the HOOH Z-matrix) |
| `MANR` | minimum | Newton–Raphson with a Morse adjustment of bond stretches; identical to `NR` for RIC in the tests |
| `RFA` | minimum | rational function step; safe when the Hessian has negative eigenvalues |
| `EVFTS` | transition state | eigenvector following (Baker), uphill along one mode and downhill along the rest |
| `MAEVFTS` | transition state | `EVFTS` with the Morse adjustment |

With `GEOM_OPT=CART` every minimum search is changed to `RFA` (with a
message). `NEGEVAL` chooses what NR and MANR do when the Hessian has a
negative eigenvalue (`ABORT`, `ABSVAL` or switch to `RFA`).
`OPT_METHOD=IGTS` and `QSD` were removed in 10/2026; they did not work.

**Hessian.** The step uses an approximate Hessian, set up by
`INIT_HESSIAN` and corrected every cycle by `HESS_UPDATE`:

- `INIT_HESSIAN=SPECIAL` (default): a model Hessian (Lindh's for RIC
  and CART, a constant diagonal for Z-matrices). `FCMINT` reads it from
  the `FCMINT` file; the file must exist. `EXACT` computes it
  analytically (SCF only; correlated methods get the SCF Hessian).
  `MOPAC` has not been tested.
- `EVAL_HESS=n` recomputes the exact Hessian every *n* cycles; `0` (the
  default with `INIT_HESSIAN=EXACT`) computes it only at the start.
  Giving `EVAL_HESS` implies `INIT_HESSIAN=EXACT`; `INIT_HESSIAN=SPECIAL`
  with `EVAL_HESS>0` stops with an error. The exact Hessian pays off for
  hard cases: the HOOH Z-matrix converges in 7 cycles with `EVAL_HESS=1`,
  8 with `EVAL_HESS=2`.
- `HESS_UPDATE`: `AUTO` (default) is `BFGS` for minima and `BOFILL` for
  transition states. `POWELL` and `PSB` are the same update (Powell's
  symmetric Broyden) and give identical results; `MS` is
  Murtagh–Sargent; `BOFILL` mixes Powell and MS; `NONE` keeps the
  starting Hessian. BFGS keeps the Hessian positive definite and is the
  only one that reliably converged the HOOH Z-matrix; for RIC and CART
  minima all of them converged in 5–10 cycles. For transition states use
  the default (`BOFILL`) or `MS`, which can build up a negative
  eigenvalue; `NONE` did not converge the HCN search. Before 10/2026
  `NONE` silently meant `BFGS`, and `MS` crashed on its first update.

**Step size.** `MAX_STEP` caps the step (millibohr, default 300).
`STP_SIZ_CTL=TRUST_RADIUS` (default) adjusts a trust radius from the
ratio of actual to predicted energy change (Ayala and Schlegel,
*J. Chem. Phys.* **107**, 375 (1997)); `NORM` scales the whole step to
`MAX_STEP` and `MAXIMUM` limits its largest component.

**Transition states.** `OPT_METHOD=EVFTS` follows the lowest Hessian
eigenvector uphill (`EIGENVECTOR=n` picks the *n*-th lowest instead),
and the Hessian is then transformed to curvilinear coordinates
automatically. The search needs a Hessian with the right negative
eigenvalue, so start from an exact Hessian (`INIT_HESSIAN=EXACT`), a
`FCMINT` file, or the reactant and product
([](#transition-states-from-reactant-and-product)). From a model Hessian
it can slide down to a minimum: the H<sub>2</sub>CO → HCOH guess ends at
H<sub>2</sub>CO with RIC and CART, but finds the transition state in 4
cycles with `INIT_HESSIAN=EXACT` (3 with `EVAL_HESS=1`; test
`ric_evfts_exact_hcoh`, reference −113.0500519970). When a search
converges at a point whose Hessian has no negative eigenvalue, the
output says `WARNING: this transition state search converged to a point
... It is NOT a transition state`. Test `126c` searches for the HCN ⇌ HNC
transition state, starting from a Hessian placed with a `%FCMINT=` file
directive (reference −93.1153442031):

```{literalinclude} ../../aces2/tests/fast/cases/126c/ZMAT
:language: text
```

Check any transition state with a frequency calculation (`VIB=EXACT`):
it must have exactly one imaginary frequency.

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
  gradient (Peng et al., *J. Comput. Chem.* **17**, 49 (1996)). Test
  `ric_constrained_130` (reference −167.9845273870) freezes the C–O and
  N–H bonds and, by the lookup rule above, the H–C–O angle at C and the
  H–N–H angle at N:

  ```{literalinclude} ../../aces2/tests/fast/cases/ric_constrained_130/ZMAT
  :language: text
  ```

- With `GEOM_OPT=CART` give the number of frozen Cartesian components and
  then their positions in the coordinate list, where atom *n* has x, y, z
  at 3*n*−2, 3*n*−1, 3*n*. Test `cart_constrained_131` (reference
  −167.9849004110) freezes five components:

  ```{literalinclude} ../../aces2/tests/fast/cases/cart_constrained_131/ZMAT
  :language: text
  ```

A minimum search in Cartesian coordinates, `GEOM_OPT=CART` with the RFA
step, is test `cart_rfa_hooh` (HOOH SCF/6-31G\*, reference
−150.7647867620; the RIC search from the same start ends 10<sup>−7</sup>
hartree away):

```{literalinclude} ../../aces2/tests/fast/cases/cart_rfa_hooh/ZMAT
:language: text
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

### Transition states from reactant and product

With
`OPT_CONTROL=QST` or `LST` and `OPT_METHOD=EVFTS`, a transition state
search with `GEOM_OPT=RIC` is guided by the synchronous transit path
between a reactant and a product (Peng and Schlegel, *Isr. J. Chem.*
**33**, 449 (1993)). The input holds three structures with the atoms in
the same order: the guess for the transition state (the main geometry),
then the reactant and the product, each after a one-line label and a
blank line (test `ric_qst3_h2co_h2_co`):

```{literalinclude} ../../aces2/tests/fast/cases/ric_qst3_h2co_h2_co/ZMAT
:language: text
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
from the same guess (−113.0500312190; the two agree to 10<sup>−8</sup>
hartree). Plain `EVFTS` from the guess wanders off. Before
10/2026 this path was broken for Cartesian input.

Without a guess for the transition state (QST2), give only two
structures: the main geometry is the reactant, and the product follows
after a label line and a blank line. `xjoda` then starts from the
midpoint of the linear synchronous transit path (Halgren and Lipscomb,
*Chem. Phys. Lett.* **49**, 225 (1977)): the product is rotated onto the
reactant, the interatomic distances are interpolated halfway, and the
Cartesians are fitted to them by least squares. The search continues as
above. For H2CO → trans-HCOH, QST2 takes 8 cycles (LST 9) to the
same transition state as a three-structure input for this reaction
and the exact-Hessian search (test `ric_qst2_h2co_hcoh`, reference
−113.0500519800):

```{literalinclude} ../../aces2/tests/fast/cases/ric_qst2_h2co_hcoh/ZMAT
:language: text
```

For vinyl alcohol →
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

### Reaction paths (IRC)

The intrinsic reaction coordinate, the steepest-descent path in
mass-weighted Cartesians from a transition state down to the minima on
either side, is followed with the Gonzalez–Schlegel method
(*J. Chem. Phys.* **90**, 2154 (1989); *J. Phys. Chem.* **94**, 5523
(1990)). Each point is found by minimizing the energy on a hypersphere
of radius *s*/2 around a pivot point half a step ahead, using analytic
gradients and an updated Hessian. The search is run with
`xa2proc PES_scan`, not `xaces2`, in a directory that holds `ZMAT` and
`GENBAS`. The `ZMAT` is a frequency calculation (`VIB=EXACT`, so SCF
for now) at the transition state, followed by an `*IRC` namelist:

```text
HCOH TS IRC
C -0.70275405  0.07016197 0.0
O  0.60245353 -0.06351953 0.0
H  0.09745863  1.02683686 0.0
H -1.28726224 -0.84915200 0.0

*ACES2(CALC=SCF,BASIS=3-21G,COORD=CARTESIAN,VIB=EXACT,SCF_CONV=10)

*IRC
STEP_SIZE=0.1
IRC_STEPS=40
IRC_SEARCH=FORWARD
```

(`PES_scan` also has single points on a grid of normal-mode
displacements and a direct-dynamics mode, chosen with `SCAN_TYPE` in a
`*PES_SCAN` namelist; these have not been tested.)

| `*IRC` keyword | Default | Meaning |
|---|---|---|
| `STEP_SIZE` | 0.1 | step *s* along the path, sqrt(amu)·bohr |
| `IRC_SEARCH` | `FORWARD` | `FORWARD` or `REVERSE`: the two directions along the transition vector |
| `IRC_STEPS` | 20 | number of points |
| `IRC_SEARCHES` | 20 | constrained optimization cycles allowed per point |
| `DO_NEXT` | `.FALSE.` | go on to the next point when one does not converge |
| `SADDLE` | `.TRUE.` | start from a transition state (`.FALSE.`: from any point, along the gradient; not tested) |

Which direction `FORWARD` takes depends on the sign of the transition
vector, so run both. `xa2proc` first runs the frequency calculation,
then one gradient calculation per constrained cycle. The output of
calculation *N* goes to `out.N`; if a script `runaces2a` (reference
point) or `runaces2b` (the other points) is on the `PATH` it is run
instead, so a site can stage the calculations elsewhere. The file
`IRC.map` collects, for every point, the path length *s*, the energy
and the Cartesian geometry in Å. The search stops after `IRC_STEPS`
points, when the path stops advancing (a minimum), or when the gradient
falls below 10<sup>−`CONVERGENCE`</sup>.

For the SCF/3-21G H<sub>2</sub>CO ⇌ trans-HCOH transition state
(−113.0500520, imaginary frequency 2707i cm⁻¹) with `STEP_SIZE=0.1`,
`REVERSE` reaches H<sub>2</sub>CO in 24 points (80 gradient
calculations; −113.2217692 at *s* = 2.30, optimized minimum −113.2218200)
and `FORWARD` reaches trans-HCOH in 20 points (71 calculations;
−113.1462658 at *s* = 1.90, optimized minimum −113.1462870). The
energy falls monotonically along both branches. The search stops when a
new point lies less than 0.1 *s* from the previous one. Test
`irc_hcoh_reverse` runs the first three `REVERSE` points of this path
(`IRC_STEPS=3`) and checks that the energy falls and that the third
point is at −113.0622774557, as in the full run.

Limitations: the points reuse the MOL file of the transition state, so
they keep its point group; a transition state whose imaginary mode is
not totally symmetric needs `SYMMETRY=NONE` in the frequency run.
`INTEGRALS=ERD` is not supported (only `xvmol` reads the geometry of
each point). The 0.3 step that was the default before 10/2026 is too
long for paths that move hydrogens: the first constrained search
walked back to the transition state.

Before 10/2026 `PES_scan` could not run: it called the
`runaces2a/b` scripts of a 2007 cluster, dropped the imaginary mode
when it read the frequencies, used a zero Hessian, looped forever in
the λ search, computed every energy and gradient with the basis
functions (and, in `xvdint`, the whole geometry) left at the transition
state, and let the points drift off the symmetry elements.

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
limits the displacements to the listed irreps.

(isotopes)=
### Isotopes

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
