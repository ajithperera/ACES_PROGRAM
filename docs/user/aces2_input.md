# ACES II input

An ACES II job runs in a directory that contains a `ZMAT` file (the input)
and a `GENBAS` file (the basis set library). Optional files such as
`ECPDATA`, `OLDMOS` or `GUESS` change how the job starts. This chapter
describes them all.

## The ZMAT file

`ZMAT` is read by `xjoda`, the first program `xaces2` runs. Lines should be
at most 80 characters; the parser does not check, and longer lines may be
cut off silently. The file has these parts, in order:

1. **Header** (optional): blank lines, comments (first non-blank
   character `#`), and file directives (first non-blank character `%`, see
   [](#file-directives)).
2. **Title**: the first line that is none of the above. Any text, up to 80
   characters.
3. **Geometry**: one line per atom, directly after the title, as a Z-matrix
   or as Cartesian coordinates. The first blank line ends it. A Z-matrix is
   followed by a blank line and then its parameter values.
4. **Namelists**: the `*ACES2` namelist, which every job needs, and any
   other namelists that some programs read (`*SIP` for ACES III, for
   example).
5. **Basis and ECP lines** (only with `BASIS=SPECIAL` or `ECP=ON`): one
   blank line after the namelist, then one line per atom (see
   [](#basis-by-atom)).

Anything after that is ignored. If nothing follows the last section, end the
file with a blank line.

Two small examples. A Z-matrix with dummy atoms:

```text
RHF CCSD property calculation of C6 in DZP
X
X 1 RX
C 2 R 1 A
C 2 R 1 A 3 T
C 2 R 1 A 4 T
C 2 R 1 A 5 T
C 2 R 1 A 6 T
C 2 R 1 A 7 T

RX=1.0
R=1.33
A=90.
T=60.

*ACES2(CALC=CCSD,PROPS=FIRST_ORDER,BASIS=DZP)

```

Cartesian coordinates:

```text
RHF MBPT(2) property calculation of formaldehyde in DZP
O 0.0  0.0          1.22
C 0.0  0.0          0.0
H 0.0  0.873489539  0.545816842
H 0.0 -0.873489539 -0.545816842

*ACES2(CALC=MBPT(2),PROPS=FIRST_ORDER,BASIS=DZP)

```

### Z-matrix (internal coordinates)

Each line places one atom relative to atoms defined before it. Atoms are
numbered by their line in the Z-matrix, dummy atoms included.

```text
A
B 1 AB
C 1 AC 2 CAB
D 3 CD 2 DCB 1 TAU
```

- Line 1: the atom symbol only.
- Line 2: symbol, the atom it is bonded to, and a distance name.
- Line 3: adds a second atom and a bond-angle name.
- Lines 4 onward: add a third atom and a dihedral-angle name.

The dihedral angle on line *I* is the rotation, seen along the axis from atom
`NCON(3I-1)` to atom `NCON(3I-2)`, that takes the projection of the bond to
atom *I* onto the projection of the bond to the third reference atom.
Clockwise is positive, and values must lie between -180° and 180°.

Rules:

- Every coordinate needs a **name**. Numbers are not allowed in the Z-matrix
  itself. Names are at most five characters.
- A name followed by `*` (for example `R*`) is **optimized** in a geometry
  optimization. Coordinates without `*` are held fixed.
- Coordinates that share a name are constrained to stay equal.
- After the Z-matrix, a blank line, then one `NAME=value` line for each
  distinct name. Text after the number is ignored, so it can hold a comment.
- Distances are in Ångström unless `UNITS=BOHR`. Angles are in degrees.
- Bond angles of exactly 0° or 180° are not allowed, because they make the
  dihedral undefined. Use a dummy atom for linear fragments.
- Two names have exact values built in, though you still have to give them
  a value: `TDA` is the tetrahedral angle, arccos(-1/3) = 109.4712...°, and
  `IHA` is the icosahedral angle, arccos(1/√5) = 63.4349...°.

**Dummy atoms** (symbol `X`) are points in space. They have no basis
functions and do not affect symmetry. Use them to avoid 180° angles and to
describe very symmetric molecules with no atom at the center.

**Ghost atoms** (symbol `GH`) have no nuclear charge but carry basis
functions, for example for counterpoise (BSSE) corrections. Ghost atoms
need `BASIS=SPECIAL` with a basis line for each atom, and only
single-point energies are supported with them.

#### The Z-matrix analyzer

`xjoda` checks the Z-matrix before a geometry optimization. It finds
coordinates that share a name but are not equivalent by symmetry, and the
reverse. It also compares the number of totally symmetric degrees of freedom
with the number of coordinates being optimized. It prints warnings for weak
Z-matrices and stops on very bad ones. Read its output and check that the
point group it reports is the one you intended.

A good rule is to include every coordinate that is symmetry-equivalent to
one you already have. For water, use both O–H distances and the H–O–H angle,
not one O–H distance plus the H–H distance. For a regular hexagon of carbon
atoms, a chain of five C–C bonds leaves the sixth bond out, so the gradient
cannot have the full symmetry and the first optimization step breaks it.
Build the ring around dummy atoms instead, as in the C6 example above. Do
not be afraid to use dummy atoms.

### Cartesian coordinates

One atom per line: the symbol and the x, y and z coordinates, free format.
`xjoda` recognizes Cartesian input automatically (a Z-matrix's first line has
only one word). `COORDINATES=CARTESIAN` forces it. Coordinates in bohr need
`UNITS=BOHR`.

### Orientation and symmetry

`xjoda` finds the full point group and runs the calculation in its largest
Abelian subgroup (D2h or one of its subgroups), unless `SYMMETRY` or
`SUBGROUP` asks for something else. Every calculation uses the
*computational orientation*, the standard orientation of that subgroup.
Orbital occupations by irreducible representation (`OCCUPATION` and similar
keywords) refer to that orientation.

For C2v, D2 and D2h the orientation is ambiguous. For water, for example,
the molecular plane may be xz or yz, which swaps b1 and b2. Run `xjoda` once
and check the orientation it chose before giving occupations by symmetry.
`SUBGRPAXIS` and `NOREORIENT` give finer control.

(zmat-namelist)=
## The `*ACES2` namelist

The namelist starts with `*ACES2` (or the older `*CRAPS`). Keywords can
follow in parentheses on the same line, or on the lines below, separated by
commas or line breaks:

```text
*ACES2(CALC=CCSD,BASIS=PVDZ,REF=UHF)
```

```text
*ACES2
# geometry optimization at the CCSD(T) level
CALC=CCSD(T),BASIS=DZP
SPHERICAL=OFF
!RESTART
```

The namelist ends at a blank line, at a line that starts with `*` (the next
namelist), at the `)` that matches the opening `(`, or at the end of the
file. Text after `#` is a comment. A line that is only a comment does not
end the namelist.

Keyword names are not case-sensitive and may be shortened to any prefix at
least as long as the minimum abbreviation shown in the
{ref}`keyword-reference` (for example `CALC` for `CALCLEVEL`). The value
depends on the keyword's type:

handle
: One name from a fixed list (`CALC=CCSD`, `REF=UHF`), case-insensitive.
  The position of the name in the list, counted from 0, is also accepted
  (`CALC=10` is `CALC=CCSD`), but names are clearer and the numbers can
  change between versions. An unknown name stops `xjoda` with
  "The value token is not recognized."

switch
: `ON` or `OFF` (also `YES`/`NO`, `TRUE`/`FALSE`, `1`/`0`). The keyword alone
  means `ON`, and a leading `!` means `OFF`: `SYM,!ECP` is the same as
  `SYMMETRY=ON,ECP=OFF`.

integer
: For example `PRINT=1`, `CHARGE=-1`.

tolerance
: An integer *N* standing for 10<sup>-N</sup>, for example `SCF_CONV=8`.

string
: Plain text, such as the basis name (`BASIS=PVDZ`), matched case-sensitively
  against `GENBAS`. Some string keywords hold arrays, described next.

**Arrays by irrep and spin** (`OCCUPATION`, `IP_SYM`, `EA_SYM`, `ESTATE_SYM`
and similar). The values for each irrep are separated by `-`, and the blocks
for each spin by `/`:

```text
OCCUPATION=3-1-1-0/2-1-1-0     # alpha / beta, 4 irreps
IP_SYM=1-0                     # one spin, 2 irreps
```

**Sets** (`DROPMO`, `FD_IRREPS`, `QRHF_GENERA`, `QRHF_ORBITA`, `QRHF_SPIN`).
Here `-` gives a range and `/` separates entries. A negative number can be
put in parentheses, `(-4)`, which is needed when it ends a range:

```text
DROPMO=1-3/10-12               # drop orbitals 1, 2, 3, 10, 11, 12
FD_IRREPS=1/3/4                # irreps 1, 3 and 4
QRHF_GEN=2/(-4),QRHF_ORB=3/1   # add an electron to irrep 2 (3rd lowest virtual),
                               # remove one from irrep 4 (highest occupied)
```

```{warning}
`DROPMO=1-31` drops orbitals 1 through 31. Very old versions of ACES II read
it as "orbitals 1 and 31".
```

**Memory.** `MEMORY_SIZE` takes a number with an optional unit: `W`
(integer words, the default) or `B` (bytes), with an optional `K`, `M` or `G`
prefix. Prefixes are powers of 1024. For example: `MEMORY_SIZE=2GB`.

**Environment variables.** A value written as `${NAME}` is replaced by the
value of the shell variable `NAME` (up to 80 characters). This makes it easy
to loop over methods or basis sets in a script:

```bash
cat > ZMAT <<'EOF'
H2 scan
H
H 1 R

R=0.7

*ACES2
CALC=${CALC}
BASIS=${BASIS}

EOF
for CALC in SCF CCSD; do
  for BASIS in PVDZ PVTZ; do
    export CALC BASIS
    rm -rf scr && mkdir scr && cp ZMAT GENBAS scr && (cd scr && xaces2 > ../$CALC.$BASIS.out)
  done
done
```

(basis-by-atom)=
## Basis sets

### Choosing a basis

There are three ways to give the basis set.

**One basis for all atoms.** `BASIS=name` uses the same basis on every
atom. `xjoda` looks up the entry `XX:name` in `GENBAS` for each element
`XX`. For example, `BASIS=CC-PVDZ` needs `C:CC-PVDZ`, `O:CC-PVDZ`,
`H:CC-PVDZ`, and so on.

```text
Formaldehyde
C   0.000000   0.000000   0.000000
O   0.000000   0.000000   1.205000
H   0.000000   0.940000  -0.580000
H   0.000000  -0.940000  -0.580000

*ACES2(CALC=SCF,BASIS=CC-PVDZ)

```

With `BASIS=SPECIAL` you name the basis yourself, after the namelist,
in one of two forms.

**One line per atom.** Leave one blank line after the namelist, then write
one line per atom in the order of the geometry, skipping dummy atoms. Each
line holds the full `GENBAS` name. Use this form when atoms of the same
element need different basis sets, or with ghost atoms.

```text
*ACES2(CALC=SCF,BASIS=SPECIAL)

C:CC-PVDZ
O:CC-PVTZ
H:3-21G
H:3-21G

```

**One line per element.** For larger molecules, give each element once.
Leave one blank line after the namelist, then a line with the number of
elements (an integer in columns 1–4), then one line per element, in any
order. `xjoda` gives every atom the basis of its element. This input is
the same as the one above:

```text
*ACES2(CALC=SCF,BASIS=SPECIAL)

3
O:CC-PVTZ
H:3-21G
C:CC-PVDZ

```

Dummy atoms are not counted and get no line, so this works with any
Z-matrix:

```text
NH3 with a dummy atom on the C3 axis
X
N 1 RN
H 2 NH 1 A
H 2 NH 1 A 3 T
H 2 NH 1 A 3 MT

RN=1.0
NH=1.012
A=112.0
T=120.0
MT=-120.0

*ACES2(CALC=SCF,BASIS=SPECIAL)

2
H:CC-PVDZ
N:CC-PVTZ

```

Rules for the per-element form:

- List every element in the molecule exactly once. A missing or repeated
  element is not reported, and the basis sets end up on the wrong atoms.
- If the number on the first line equals the number of lines in the
  geometry (dummy atoms included), `xjoda` reads the lines that follow as
  one per atom, in geometry order. This happens when every atom is a
  different element (CO, HCN, ...). For those molecules, write the lines
  in geometry order, or use the one-line-per-atom form.

`xjoda` does not check that the element in each line matches the atom it
is given to. To see which basis sets exist for an element, search `GENBAS`
for `XX:`.

### ECP names

With `ECP=ON` (which needs `BASIS=SPECIAL`), the ECP names from `ECPDATA`
follow the basis lines: one blank line, then the ECP names in either of the
same two forms. Atoms without an ECP get `XX:NONE`. Both forms below
describe CuH₂⁻ with the Hay–Wadt ECP on copper:

```text
CuH2- with an ECP on Cu
CU  0.000000   0.000000   0.000000
H   0.000000   0.000000   1.520000
H   0.000000   0.000000  -1.520000

*ACES2(CALC=SCF,BASIS=SPECIAL,ECP=ON,CHARGE=-1)

CU:HAY-WADT-VDZ
H:CC-PVDZ
H:CC-PVDZ

CU:HAY-WADT-VDZ
H:NONE
H:NONE

```

```text
*ACES2(CALC=SCF,BASIS=SPECIAL,ECP=ON,CHARGE=-1)

2
CU:HAY-WADT-VDZ
H:CC-PVDZ

2
H:NONE
CU:HAY-WADT-VDZ

```

The basis and ECP sections are independent: either can use either form,
and the elements need not be in the same order in both.

### GENBAS and ZMAT.BAS

`GENBAS` is the basis set library. On its first run, `xjoda` copies the
entries the job needs into `ZMAT.BAS` and reads that file instead of
`GENBAS` from then on. Delete `ZMAT.BAS` if you change the basis.

Each entry has this fixed format:

| Line | Format | Content |
|---|---|---|
| 1 | A80 | blank |
| 2 | A80 | name, `XX:BASIS` |
| 3 | A80 | comment |
| 4 | A80 | blank |
| 5 | I3 | number of shells, NS |
| 6 | NS(I5) | angular momentum of each shell |
| 7 | NS(I5) | number of contracted functions in each shell (NC) |
| 8 | NS(I5) | number of primitives in each shell (NE) |
| 9 | A80 | blank |
| 10 | NE(F14.7) | exponents of the shell |
| 11 | A80 | blank |
| 12 | NC(F10.7,1X) | contraction coefficients, one row per primitive |
| 13 | A80 | blank |

Lines 10–13 repeat for each shell. Shells must be in order of increasing
angular momentum (s first, then p, and so on).

### ECPDATA

Effective core potentials are read from `ECPDATA`. Each entry looks like this:

```text
*
CU:ECP-10-SK
# Stevens/Krauss ECP for Cu, 10 core electrons, LMAX = 2
*
NCORE = 10 LMAX = 2
d
-10.00000000 1 511.9951763
-72.55482820 2  93.2801074
-12.74502310 2  23.2206669
s-d
  3.00000000 0 173.1180854
 ...
*
```

`NCORE` is the number of core electrons the ECP replaces, and `LMAX` is the
highest angular momentum of the projectors. Each block (`d`, `s-d`, `p-d`,
...) is one angular term of the potential, given as Gaussian terms, one per
line: the coefficient, the power of *r*, and the Gaussian exponent. The form
follows L. R. Kahn, P. Baybutt and D. G. Truhlar, *J. Chem. Phys.* **65**,
3826 (1976).

## Starting orbitals and other input files

`OLDMOS`, `NEWMOS`
: `xvscf` writes the final SCF orbitals (in the symmetry-adapted basis) to
  `NEWMOS`. Copy it to `OLDMOS` to start a later job on the same molecule
  and basis from those orbitals, with `GUESS=READ_SO_MOS` and a `GUESS`
  file (below).

`AOBASMOS`, `OLDAOMOS`
: The same, but in the AO basis before symmetry adaptation
  (`GUESS=READ_AO_MOS`). This works even when the point group changes, as it
  does between finite-difference displacements.

`GUESS`
: Controls the occupation of the starting orbitals when they are read with
  `GUESS=READ_SO_MOS`. Every field must be given; there are no defaults. See
  the format below.

`FCMINT`
: A force constant matrix in internal coordinates, used as the starting
  Hessian of a geometry optimization. Optimizations and frequency jobs also
  write `FCM`, `FCMSCR` and `FCMFINAL`.

`ISOMASS`
: Atomic masses for vibrational frequencies (free format, one per real atom
  in Z-matrix order). If the file exists, its masses replace the standard
  ones.

`frequency`
: The frequencies at which dynamic polarizabilities are computed (free format).

### The GUESS file

This example is for a state of the water cation with four irreps (`NIRREP=4`)
and UHF orbitals (`NSPIN=2`). Text after the numbers on each line is ignored.

```text
H2O+ test.
  3  1  1  0        alpha occupation by irrep
  3  0  1  0        beta occupation by irrep
  0  0  0  0        alpha pairs to swap, irrep 1
  0  0  0  0        alpha, irrep 2
  0  0  0  0        alpha, irrep 3
  0  0  0  0        alpha, irrep 4
  0  0  0  0        beta pairs to swap, irrep 1
  0  0  0  0        beta, irrep 2
  0  0  0  0        beta, irrep 3
  0  0  0  0        beta, irrep 4
  0  0  0  0        alpha lock flags by irrep
  0  0  0  0        beta lock flags by irrep
  0  0  0  0        alpha print flags by irrep
  0  0  0  0        beta print flags by irrep
  0  0              stop flags
  1                 read orbitals from OLDMOS
  1                 copy alpha orbitals to beta
  0                 read GUESS only for the first SCF
```

In order:

1. A title line.
2. Occupation by irrep, one line per spin (format 8I3).
3. Orbital swaps, one line per irrep for each spin (format 4I3): two pairs
   of orbital numbers, so at most two swaps per irrep and spin.
4. Lock flags by irrep, one line per spin. A positive value keeps the
   occupation of that irrep following the starting orbitals ("minimum
   change", using C<sub>old</sub><sup>T</sup> S C) during the SCF.
5. Print flags by irrep, one line per spin. A positive value prints the
   starting orbitals of that irrep.
6. Two stop flags. A positive first value stops after the starting
   orbitals have been formed.
7. Read flag. A positive value reads the starting orbitals from `OLDMOS`.
8. Copy flag. A positive value reads only the alpha orbitals and copies them
   to beta, for example to start UHF from closed-shell RHF orbitals.
9. Reuse flag. A positive value makes every SCF in a multi-step job (an
   optimization, a stability analysis) read `GUESS` again. With 0, the file
   is read once and then deleted.

(file-directives)=
## Where files go: file directives

By default every file a job creates is written to the directory where
`xaces2` was started. A `%` line in the ZMAT header moves any file except
`ZMAT` itself:

```text
%MOINTS=/scratch/me/job1/MOINTS
%JOBARC=/home/me/safe/JOBARC
%JAINDX=/home/me/safe/JAINDX

Water SCF frequencies
...
```

The `SAVEDIR` directive names the directory used for restart archives of
geometry optimizations and finite-difference jobs, given as an absolute path
or relative to the working directory:

```text
%SAVEDIR=/home/me/restart/job1
```

The default is a directory named `SAVEDIR` in the working directory. Jobs
that run at the same time must not share a save directory, or they will
overwrite each other's files.

## Files a job creates

`JOBARC`, `JAINDX`
: The record store shared by all member programs: geometry, basis
  information, orbitals, energies, gradients and more. See the programmer
  manual.

`MOINTS`, `GAMLAM`, `MOABCD`, `DERINT`, `DERGAM`
: The lists of integrals and amplitudes used after the SCF. These are by far
  the largest files.

`MOL`, `IIII`, `IIJJ`, `IJIJ`, `IJKL`
: `MOL` holds the molecule and basis for the integral program; the others
  hold AO integrals. `IIII` grows quickly with basis size.

`HF2`, `HF2AA`, `HF2AB`, `HF2BB`
: Partly transformed integrals, normally deleted by `xintprc`.

`OPTARC`, `OPTARCBK`
: The history of a geometry optimization, and a backup of it.

`DIPOL`, `POLAR`, `DIPDER`, `POLDER`
: Dipole moments, polarizabilities, and their derivatives.

`TGUESS`, `LGUESS`
: Coupled-cluster T and Λ amplitudes kept for restarts.

`VPOUT`
: One-electron property integrals.

Use a fresh directory for each job, or delete these files between jobs.
Files left over from an earlier job (`JOBARC` in particular) are read by the
next one as if they belonged to it.
