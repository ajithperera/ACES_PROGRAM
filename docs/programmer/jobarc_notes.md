<!-- Descriptions of JOBARC records for the generated catalog (programmer/jobarc.md).
     '## NAME' starts a record (the 8-character record name, trailing blanks
     dropped). Which directories write and read each record comes from the code
     (docs/tools/gen_jobarc.py); do not repeat it here. Seeded 2026-10-23 from
     the ACES II Programmer's Manual (1990s), with obvious errors corrected. -->

## ANGMOMBF
Angular momentum quantum number *l* (s = 0, p = 1, ...) of each AO basis function, in the computational order. See also `ANMOMBF0`. Integer, NAO.

## ANMOMBF0
Angular momentum quantum number *l* of each AO basis function, in ZMAT order. See also `ANGMOMBF`. Integer, NAO.

## AO2SO
Transformation between the symmetry-adapted orbitals (columns) and the atomic orbitals (rows). Multiplied from the left, it takes an SO-basis object to the AO basis, so the name is the reverse of what it does. In general not unitary and often not square. The inverse is `AO2SOINV`. Floating point, NAO × NSO.

## AO2SOINV
Transformation between the atomic orbitals (columns) and the symmetry-adapted orbitals (rows). Multiplied from the left, it takes an AO-basis object to the SO basis. In general not unitary and often not square. Floating point, NSO × NAO.

## ATOMCHRG
Atomic number of every center in the Z-matrix; zero for dummy atoms. NATOMS.

## ATOMMASS
Atomic mass of every center in the Z-matrix; zero for dummy atoms. Floating point, NATOMS.

## BRUKTEST
Brueckner convergence flag read by the `xaces2` driver: 1 when the Brueckner orbitals have converged to `BRUCK_CONV`, 0 to run another iteration. Integer, 1.

## CENTERBF
Atomic center of each AO basis function; centers and functions both in the computational order. See also `CNTERBF0`. Integer, NAO.

## CMP2ZMAT
Transformation from the computationally ordered SO basis (columns) to the ZMAT-ordered AO basis (rows), acting from the left. In general not unitary and often not square. The inverse is `ZMAT2CMP`. Floating point, NAO × NSO.

## CNTERBF0
Atomic center of each AO basis function, in ZMAT order. See also `CENTERBF`. Integer, NAO.

## COORD
Cartesian coordinates of every center in the Z-matrix, in the computational orientation (bohr). Floating point, 3 × NATOMS.

## FULLAOSO
Like `AO2SOINV` (AO to SO, from the left) but always square: the SO basis here keeps functions such as spherical-harmonic contaminants that the computational basis drops. The inverse is `FULLSOAO`. Floating point, NAO × NAO.

## FULLSOAO
Like `AO2SO` (SO to AO, from the left) but always square, including the functions the computational basis drops. The inverse is `FULLAOSO`. Floating point, NAO × NAO.

## GRADIENT
Cartesian gradient, ordered ∂E/∂x, ∂E/∂y, ∂E/∂z atom by atom, atoms in the computational order. Floating point, 3 × NATOMS.

## HESSIANM
Cartesian Hessian (force constant) matrix, F(1x,1x), F(1x,1y), ..., atoms in the computational order. Floating point, 9 × NATOMS².

## IFLAGS
The first block of job-control flags set from the `*ACES2` namelist (see `joda/asv_ctl.c` and `include/flags.h`). Written by `xjoda`; later programs may change it. `IFLAGS2` holds the second block. Integer.

## IFLAGS2
The second block of job-control flags (see `IFLAGS`). Integer.

## LINEAR
1 if the molecule is linear, 0 otherwise. Integer, 1.

## MAP2ZMAT
Maps centers from the computational order to the ZMAT order: for computational center *i*, `MAP2ZMAT(i)` is its number in the ZMAT file. Integer, NATOMS.

## MODROPA
Alpha orbitals dropped from the correlation treatment (`DROPMO`), numbered by orbital energy from the lowest (1). Integer, `NUMDROPA`.

## MODROPB
Beta orbitals dropped from the correlation treatment. Integer, `NUMDROPB`.

## NAOBASFN
Number of AO basis functions. Because the SO basis can omit functions (spherical-harmonic contaminants), this can be larger than `NBASTOT`. Integer, 1.

## NAOBFORB
Number of AO basis functions on each symmetry-unique center ("orbit"). The AO analog of `NBASATOM`. Integer, NOrbits.

## NBASATOM
Number of symmetry-adapted functions on each symmetry-unique center. Integer, NOrbits.

## NATOMS
Number of centers in the Z-matrix, dummy atoms included. Integer, 1.

## NREALATM
Number of real (non-dummy) centers. Integer, 1.

## NUCREP
Nuclear repulsion energy (hartree). Floating point, 1.

## NBASTOT
Number of symmetry-adapted basis functions. See also `NAOBASFN`. Integer, 1.

## NMPROTON
Number of protons (sum of nuclear charges). Integer, 1.

## NOCCORB
Number of occupied orbitals in the correlated calculation, for each spin. Integer, 2.

## NVRTORB
Number of virtual orbitals in the correlated calculation, for each spin. Integer, 2.

## NUMBASIR
Number of SOs in each irrep. Integer, NIRREP.

## NUMDROPA
Number of alpha orbitals dropped from the correlation treatment. Integer, 1.

## NUMDROPB
Number of beta orbitals dropped from the correlation treatment. Integer, 1.

## OCCUPYA
Number of occupied alpha orbitals in each irrep. Integer, NIRREP.

## OCCUPYB
Number of occupied beta orbitals in each irrep. Integer, NIRREP.

## ORIENTMT
Rotation between the computational and canonical orientations. Floating point, 3 × 3.

## QRHFIRR
Irreps of the orbitals whose occupation is changed in a QRHF calculation. A negative value removes an electron from an orbital of that irrep, a positive one adds an electron. QRHFTOT.

## QRHFLOC
Position of each QRHF orbital within its irrep. For removal, 1 is the highest occupied orbital and larger values go down from the Fermi level; for addition, 1 is the lowest virtual and larger values go up. QRHFTOT.

## QRHFTOT
Number of orbitals whose occupation changes in a QRHF calculation. 1.

## S2SCF
Expectation value ⟨S²⟩ of the reference wavefunction. Floating point, 1.

## SCFENEG
SCF energy. Floating point, 1.

## SCFEVALA
Diagonal of the alpha Fock matrix in the MO basis (orbital energies), in the post-SCF orbital order. Floating point, NMO.

## SCFEVALB
As `SCFEVALA`, for beta spin.

## SCFEVCA0
Alpha MO coefficients in the SCF orbital order. Floating point, NAO × NAO.

## SCFEVCB0
As `SCFEVCA0`, for beta spin.

## SCFEVECA
Alpha MO coefficients in the post-SCF orbital order, as used by the correlated methods. Floating point, NMO × NAO.

## SCFEVECB
As `SCFEVECA`, for beta spin.

## SCFKICK
Flag read by the `xaces2` driver when `HFSTABILITY=FOLLOW`: 1 reruns the SCF from rotated orbitals after an instability that can be followed within the current symmetry and spin restrictions, 0 ends the job. Integer, 1.

## TOTENERG
Total energy of the current method. First written by `xvscf` and then updated by each program that adds to the energy, so after a job it holds the final energy. This is the record the test suites compare. Floating point, 1.

## UHFRHF
SCF type: 0 RHF, 1 UHF, 2 ROHF. First written by `xjoda`; later programs may change it. 1.

## ZMAT2CMP
Transformation from the ZMAT-ordered AO basis (columns) to the computationally ordered SO basis (rows), acting from the left. In general not unitary and often not square. The inverse is `CMP2ZMAT`. Floating point, NSO × NAO.
