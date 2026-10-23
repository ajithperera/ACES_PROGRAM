<!-- Keyword descriptions for the generated keyword reference.
     '# Category' starts a group, '## KEYWORD' (the full name as in
     aces2/joda/asv_ctl.c, without the '#') starts a description.
     Names, types, defaults and value lists come from the code; do not
     repeat them here.  Seeded 2026-10-23 from the ACES II 2.5.0 user
     manual (2006); text not yet reviewed against the current code. -->

# System: general

## CALCLEVEL

Selects the method for the energy (and its derivatives) after the SCF. See {doc}`../user/methods` for the method families and examples. `CALC=HFDFT` is obsolete; use `SCF_TYPE=KS` for Kohn–Sham DFT. Not every value is implemented for every reference, derivative order or property.

## MEMORY_SIZE

Sets the “maximum” memory to allocate. Most member executables will not allocate more than the amount specified, but some will attempt to allocate extra bits (on the order of a few megabytes) that are of little consequence compared to the main memory heap. There is special parsing logic for the corresponding value string. The tail of the string is checked for case-insensitive units and optional prefixes. Units can be W (integer words) or B (bytes), and prefixes can be K, M, or G (kilo, mega, and giga, respectively). The number part is read with the atof() C function. [NOTE: The prefixes scale the units by 2<sup>10</sup> (1024), not 1000, and the default unit is integer words (4 bytes for 32-bit binaries and 8 bytes for 64-bit binaries).]

## RESTART

Controls the creation, updating, and reading of a checkpoint directory. Currently, only coarse-grain restarts are available, which means geometry optimizations and finite difference calculations can be restarted from the last checkpoint. Fine-grain restarts are not available yet but will also checkpoint SCF, T, and Lambda iterations.

## GRAD_CALC

Specifies whether to calculate a gradient analytically or numerically. This should be used only to override analytical gradients since the program will attempt to use them whenever possible.

## DERIV_LEV

Specifies whether or not derivatives of the energy are to be calculated and if so then whether first or second. =NONE Derivatives not calculated, =FIRST First derivatives to be calculated, =SECOND Second derivatives to be calculated. This is automatically set.

# System: molecular control

## SYMMETRY

Specifies which subgroup (“computational” point group) of the full point group is to be used in the energy and/or gradient calculation. =OFF forces a no symmetry run (in C1), =ON runs the calculation in the largest self-adjoint subgroup ( D2h and its subgroups), and =FULL uses the full point group. Currently, ACES II does not support groups with degenerate representations, so the FULL option has no value unless joda is used to make input files for another program package.

## SUBGROUP

Specifies a lower computational point group symmetry to use (provided it is a subgroup of the full group). Acceptable values are: DEFAULT, C1, C2, CS, CI, C2V, C2H, D2, and D2H. SUBGROUP=C1 is equivalent to SYMMETRY=OFF. The DEFAULT option uses the highest order Abelian subgroup (including the full group).

## SUBGRPAXIS

This is a somewhat complicated keyword to use. The value can be X, Y, or Z. The use of the keyword is best described by example: Suppose one is running a calculation on water, and wishes to run it in the Cs point group with the “special” plane being the one that bisects the H-O-H bond angle. SUBGRPAXIS specifies which Cartesian direction in the C2v frame becomes the special direction in the Cs frame. Normally, joda will orient water in the yz plane. In this example, the y axis in the C2v frame should be the z axis in the Cs frame, so SUBGRPAXIS=Z. When the true Abelian subgroup is either C2v or D2h, the orientation is not well defined, and it might be necessary to run xjoda twice. If SUBGROUP=0 in the first pass, then the reference orientation for the true Abelian subgroup can be determined and the appropriate value of SUBGRPAXIS can be selected.

## NOREORI

Forces the program not to change the internal orientation of the molecule specified by the input Cartesian coordinates. To function correctly, this keyword also requires turning off symmetry (SYMMETRY=OFF).

# System: debug

## PRINT

Controls the amount of printing in all member executables except joda. A value of 1 will produce a modest amount of additional output over the default value, which includes some useful information such as SCF eigenvectors.

## JODA_PRINT

Controls the amount of debug printing performed by joda. The higher the value, the more information is printed. Values of 25 or higher generally do not produce anything of interest to the general user, and values greater than 999 will dump the core vector to disk.

# I/O subsystem

## FILE_RECSIZ

Sets the length of the physical records used in direct access file I/O. This value should always be chosen as a multiple of 512 bytes and is parsed like the MEMORY SIZE keyword. The default value (-1) determines a size based on the value of MEMORY SIZE.

## CACHE_RECS

Sets the number of physical records (of length FILE RECSIZE) held in the I/O cache. The maximum number of records is 128, and the default value (-1) determines a number based on the value of MEMORY SIZE.

## INCORE

This keyword can be used to reduce disk I/O (usually by a significant amount). Acceptable values are: NONE, NOABCD, NOABCI, NOABIJ, ALL, and T. =NONE and =ALL load no and all lists into core memory, respectively. =NOABCD, =NOABCI, and =NOABIJ load all lists that are at most the size of <Ab|Ci>,<Ab|Ij>, and<Ij|Ka> molecular integral lists, respectively. These options are implemented in vcc, lambda, dens, anti, bcktrn, and vee. The =T option only applies to vcc, which will load T1 and T2 amplitudes into memory. All of the other executables will treat =T as =NONE (including lambda, although this should change soon).

# Chemical system

## COORDINATES

Specifies whether the molecular system is defined in INTERNAL (Z matrix) or CARTESIAN (xyz) coordinates.

## UNITS

Specifies the units used for molecular geometries. The value can be ANGSTROM or BOHR.

## CHARGE

Sets the molecular charge.

## MULTIPLICTY

Sets the spin multiplicity (2 S + 1).

# Basis set

## BASIS

Defines the name of the basis set to use for all atoms in the molecule. If BASIS=SPECIAL, then each atom (in coordinate order) must be specified after the *ACES2 namelist.

## SPHERICAL

Controls whether (=ON) spherical harmonic (5d, 7f, 9g, ...) or (=OFF) Cartesian (6d, 10f, 15g, ...) basis functions are used.

## LINDEP_TOL

Sets the tolerance for linear dependencies in the basis set. The basis set is considered linearly dependent and eigenvectors of the overlap matrix are neglected if the associated eigenvalues are less than 10<sup>-N</sup> .

## GENBAS_1

This keyword applies to first-row elements (H and HE) and specifies the number of contracted Gaussian functions per shell. There is usually no need to use this keyword, but it can be useful for using a subset of the functions in a particular entry in the GENBAS file, particularly for generally contracted basis sets. For example, if entry H:BASIS in the GENBAS file contains 7 contracted s functions, 4 p functions, and one d function, then setting GENBAS 1=730 would eliminate the last p function and the d function. The default for this keyword is to use the unaltered entry in GENBAS.

## GENBAS_2

This keyword performs the same function as GENBAS 1 but applies only to second-row atoms.

## GENBAS_3

This keyword performs the same function as GENBAS 1 and GENBAS 2, but applies only to third-row atoms.

## CONTRACTION

Specifies the contraction scheme used by the integral and integral derivative programs. This can be set to SEGMENTED or GENERAL. [NOTE: Even for truly segmented basis sets, both integral and integral derivative programs run significantly faster in GENERAL mode.]

## ECP

Controls whether effective core potentials (a kind of pseudopotential) are used (=ON) or not (=OFF). ECP=ON requires BASIS=SPECIAL and specification of the ECP data sets (in a file named ECPDATA).

# Integrals

## INTEGRALS

Specifies which integral package to use. This is not a very robust keyword and should always be set to VMOL unless otherwise directed to change it. Other values include SEW ARD and GAMESS, but these options will not work unless the binaries were linked to the appropriate libraries, which require third-party license agreements.

=ERD computes the integrals with the OED/ERD library (`xerdints`) instead of `xvmol` + `xvmol2ja`. It writes the same files and JOBARC records, for Cartesian and spherical basis sets, ECPs and DKH. As with VMOL, it cannot be combined with DIRECT or SYMMETRY=FULL. Only the xaces2 driver uses it so far; the other drivers (runpccd, aux, pyaces) still always run VMOL.

## INTGRL_TOL

Sets the tolerance for storing the two-electron integrals such that integrals having an absolute value greater than 10<sup>-N</sup> will be stored on disk.

## DIRECT

Controls whether two-particle AO integrals are calculated on-the-fly or calculated once and stored in a file. This cannot be used with VMOL integrals.

# Reference

## REFERENCE

Specifies the type of SCF reference to use. Supported references include spin-restricted Hartree-Fock (=RHF), spin-unrestricted Hartree-Fock (=UHF), and spin-restricted open-shell Hartree-Fock (=ROHF). These apply to Kohn-Sham DFT references as well (e.g., RHF is spin-restricted Kohn-Sham). REF=TWODET (two-determinant reference for open-shell singlet CC calculations) is available, but the orbitals are obtained from a closed-shell state and are not optimum open-shell singlet (or low spin triplet) SCF orbitals. The orbital occupancy for subsequent CCSD calculations is specified by the various QRHF keywords.

## NONHF

Signals the correlation energy code if the reference wavefunction satisfies the HartreeFock equations ( Fai = 0). Usually there is no need to set this keyword since standard non-HF reference functions (QRHF and ROHF) set this flag internally; however, if you input a set of orbitals to the correlation energy code directly, then it might be necessary to set NONHF=ON.

## ORBITALS

Specifies semicanonical orbitals ( Fi⁄=j = Fa⁄=b = 0) in non-HF calculations. “Semicanonical” orbitals are obtained by diagonalizing the occupied-occupied ( ij) and virtualvirtual (ab) blocks of the spin-orbital Fock matrix and can be advantageously exploited in certain post-SCF calculations. (This is particularly true for ROHF-MBPT and noniterative triple excitation corrections.) There is no default value for this parameter, and considerable logic is used to determine the orbital type in post-SCF non-HF calculations (if the keyword is not defined). It is strongly recommended that this keyword not be used by anyone who is not thoroughly familiar with non-HF CC/MBPT methods, since the logic used to set the default value is sound. =STANDARD uses the orbitals obtained in the reference function calculation without modification, and =SEMICANONICAL forces a transformation to semicanonical orbitals.

## PERT_ORB

Specifies whether the gradient formulation assumes that the perturbed orbitals are not those in which the Fock matrix is diagonal (=STANDARD). =CANONICAL means that the perturbed orbitals are assumed to be canonical. This keyword must be set to CANONICAL in gradient calculations with methods that include triple excitations (i.e., MBPT(4), CCSD[T], CCSD(T), and QCISD(T)).

# SCF: general

## SCF_TYPE

Switches the SCF code between the standard Hartree-Fock program (=HF) and the Kohn-Sham DFT program (=KS). If SCF TYPE=KS, then the KS-DFT program generates the SCF reference and requires separate namelists (viz. *VSCF and *INTGRT).

## FOCK

Specifies the algorithm for constructing the Fock matrix in SCF calculations. =PK uses the PK-supermatrix approach while =AO constructs the matrix directly from the AO integrals. In general, PK is faster but results in considerable use of disk space when out-of-core algorithms are required.

## CHECK_SYM

Specifies the action taken when the density matrix does not transform as the totally symmetric irreducible representation of the full molecular point group. =NORMAL terminates the run if the molecule is nonlinear, while =OVERRIDE allows the job to continue but prints a warning message. Nonsymmetric density matrices result from calculations on electronically degenerate states or from broken-symmetry SCF solutions. Often, such calculations do not give meaningful results and inexperienced users are encouraged to use CHECK SYM=NORMAL in their calculations on nonlinear molecules. For Π, ∆, Φ, ... states of linear molecules, however, meaningful calculations can still be performed even though the density matrix is not symmetric.

## SCF_PRINT

[This keyword is currently not used but exists for future compatibility.]

# SCF: orbital control

## GUESS

Specifies where the initial SCF eigenvectors are read from. The HF-SCF executable checks multiple places for pre-existing orbitals, but this keyword can override that logic. =MOREAD→ JOBARC file =CORE→ core Hamiltonian =READ SO MOS→ OLDMOS file (from NEWMOS) =READ AO MOS→ OLDAOMOS file (from AOBASMOS) Other options include =NDDO, =W ALT PRJDEN, =MIN BASIS, and =HUCKEL, although some of these are still experimental and not fully supported.

## OCCUPATION

Specifies the orbital occupancy of the reference wavefunction in terms of the occupation numbers per irrep per spin. The occupancy is specified by NIRREP (or 2*NIRREP) integers that define the number of occupied orbitals of each symmetry type. (NIRREP is the number of irreducible representations in the computational point group.) If there are no orbitals of a particular symmetry type, then a zero must be entered. If the reference wavefunction is for an open-shell system, then two strings of NIRREP occupation numbers separated by a forward-slash (/) define the α and β sets of orbitals. An example of the use of the OCCUPATION keyword for the water molecule would be OCCUPATION=3-1-1-0. For the 2A1 water cation (an open-shell system), the keyword would be specified by OCCUPATION=3-1-1-0/2-1-1-0. It should be noted that the VMOL integral program orders the irreducible representations in a strange way, which most users do not perceive to be a logical order. Hence, it is advised to run a single-point SCF calculation to determine the initial number and ordering of the irreducible representations. The occupation keyword may be omitted, in which case an initial orbital occupancy is determined by diagonalizing the core Hamiltonian. In many cases, SCF calculations that use the core Hamiltonian guess will converge to the lowest energy SCF solution, but this is not guaranteed.

## LOCK_ORBOCC

Controls orbital occupancies among symmetry blocks in the SCF iterations. =ON locks the occupation to the OCCUPATION value array (or the initial guess if OCCUPATION is undefined); =OFF permits the occupation to change. [NOTE: If the OCCUPATION keyword is defined, then LOCK ORBOCC is switched on automatically.]

## NEWVRT

Signals if the SCF virtual orbitals are to be replaced by a set determined from an N-1 potential. This is an orthogonal transformation within the virtual space. As long as appropriate fij and fab are included, the energies of standard single-reference methods are unchanged. However, for TD-CC methods, this keyword mixes one of the occupied orbitals with the virtual space and the results are changed. It is anticipated that TDCCSD calculations will be improved by the use of this keyword, but this has yet to be demonstrated. NEWVRT is useful for interpreting results of EOM-EE calculations since excitations tend to be more easily identified as involving one occupied and one virtual orbital. =OFF does not rotate the virtual space, and =ON rotates the virtual space.

# SCF: iteration control

## SCF_CONV

Sets the convergence criterion for the SCF equations. Equations are considered converged when the maximum change in density matrix elements is less than 10<sup>-N</sup> .

## SCF_MAXCYC

Sets the maximum number of SCF iterations.

## DAMP_TYP

Specifies what type of damping is used during the SCF iterations. The value can be NONE (no damping), DA VIDSON (Davidson’s empirical dynamical damping scheme), and STATIC (a fixed damping parameter). Note that RPP convergence extrapolation is not turned on until the damp factor has gone below DAMP TOL and the energy change is below 0.05 a.u. DAMP TYP=DA VIDSON is recommended even though NONE is the default.

## DAMP_TOL

If DAMP TYP=DA VIDSON, then the cutoff is defined to be N∗ 0.01, so the default cutoff is 0.1.

## DAMPSCF

If DAMP TYP=STATIC, then the static damping factor used in the SCF iterations is N∗ 0.01, so the default damping factor is 0.2.

## SCF_EXTRAP

(Formerly RPP.) Controls whether or not the reduced partitioning procedure is to be used to accelerate convergence of the SCF equations.

## SCF_EXPORDE

(Formerly RPP_ORDER.) Sets the number of density matrices to be used in the RPP convergence acceleration procedure.

## SCF_EXPSTAR

(Formerly RPP_LATEST.) Sets the latest iteration for initiation of the RPP convergence acceleration procedure. RPP is switched on when the error falls below a certain threshold, but in difficult cases, in which the iterations are oscillatory, it is necessary to force it on when the error is still large. In these cases, the RPP will begin on the iteration number specified by this parameter.

## LSHF_A1

Sets the doubles-singles level shifting parameter α, in which α = N∗0.01. This keyword is currently only meaningful in ROHF calculations.

## LSHF_B1

Sets the singles-virtuals level shifting parameter β, in which β = N∗0.01. This keyword is currently only meaningful in ROHF calculations.

# SCF reference adjustments

## DROPMO

Lists which molecular orbitals will be dropped from the post-SCF calculation. The orbitals are ordered by eigenvalue from the most stable (negative energy) to the most unstable (largest positive energy) regardless of irrep. The delimiter is a forward-slash (/) that separates single orbitals and orbital ranges (x-y inclusive). For example, the string “1-5/55/62-64” will cause vtran to drop orbitals 1, 2, 3, 4, 5, 55, 62, 63, and 64. For UHF calculations, the appropriate orbitals are deleted from both spin cases.

## HFSTABILITY

Checks the stability of RHF and UHF wavefunctions and optionally rotates the orbitals to a lower SCF solution. There are three possible options for this keyword. =OFF does nothing and =ON performs a stability check and returns the number of negative eigenvalues in the orbital rotation Hessian. =FOLLOW performs the stability check and then proceeds to rotate the SCF orbitals in the direction of a particular negative eigenvalue of the orbital rotation Hessian (see the explanation of the ROT EVEC keyword below), after which the SCF is rerun.

## ROT_EVEC

Defines which eigenvector of the orbital rotation Hessian will be used to rotate the original SCF orbitals. By default, it will use the eigenvector with the lowest eigenvalue of irrep 1 (the totally symmetric part of the block-factored Hessian). This choice often leads to the lowest energy SCF solution. For RHF stability checks, only those instabilities which correspond to RHF solutions will be considered. It is important to understand that following non-symmetric eigenvectors lowers the symmetry of the wavefunction and that following RHF →UHF stabilities leads to a UHF solution. To converge the SCF roots associated with such instabilities, one must run the calculation in reduced symmetry and as a closed-shell UHF case, respectively. If not set to 0, the program follows the vector with the N th lowest eigenvalue having the proper symmetry (totally symmetric) and spin (RHF →RHF or UHF→UHF) properties.

## BRUECKNER

Controls whether Brueckner orbitals are to be determined for the specified CC method.

## BRUCK_CONV

Sets the convergence criterion in Brueckner iterations. The calculation is considered converged when the largest single excitation amplitude falls below 10<sup>-N</sup> .

## QRHF_GENERA

The presence of this keyword turns on QRHF logic, and the value of each element specifies which irrep is created (positive) or annihilated (negative). In the QRHF scheme, an RHF (closed-shell) SCF calculation is performed with the OCCUPATION keyword or the CHARGE and MULTIPLICITY keywords. The correlated calculation is then performed on an open-shell reference generated from this closed-shell state by removing, adding, or exciting electrons. Any number of electrons may be added and all of the electrons may be removed from the SCF reference. The array elements associated with this keyword must be separated by forward-slashes (/) and be positive or parenthesized negative. The absolute value of each element specifies the symmetry block involved in the addition or removal of electrons. The numerical ordering of the symmetry blocks is consistent with then OCCUPATION keyword (the symmetries of the integrals). By default, the electrons of each irrep are added to the lowest unoccupied α orbital in the symmetry block and removed from the highest occupied β orbital. Different orbitals and spins can be specified with the QRHF ORBITAL and QRHF SPIN keywords. [NOTE: Gradients and property calculations are currently available only for cases involving addition or removal of electrons. Mixed cases involving both processes are not available. Gradients and properties are available for open-shell singlet CCSD wavefunctions for the case that the open-shell orbitals have different symmetries.] ???

## QRHF_ORBITA

By default, in QRHF calculations, electrons are removed from the highest occupied β orbital in a symmetry block (symmetry block HOMO), while electrons are added to the lowest unoccupied α orbital within a symmetry block (symmetry block LUMO). The purpose of the QRHF ORBITAL keyword is to allow additional flexibility in choosing which orbitals will have their occupation numbers altered. The value of this keyword gives the index with respect to the default orbital for the orbital which will be populated or depopulated. For calculations involving more than one removal or addition of electrons, values are separated by forward-slashes and correspond one-for-one to the QRHF GENERAL value array. For example, specifying QRHF GENERAL=2/(-4) and QRHF ORBITAL=3/2 means that an electron will be added to the third lowest virtual in symmetry block 2 and another will be removed from the second highest occupied orbital in symmetry block 4. Examples in Sections 9.1.4 and ?? (pages 71 and ??) further illustrate the QRHF input options.

## QRHF_SPIN

Specifies the spins of QRHF occupations by overriding the remove-fromβ and add-to-α behavior. An element value of 1 means α spin and 2 means β spin. [NOTE: This option allows one to construct low-spin determinants, which generally are unsuitable for single-reference coupled-cluster calculations. An important exception is the open-shell singlet coupled-cluster method (REFERENCE=TWODET).]

## UNO_REF

Controls the use of the new UNO reference state. The ACES II calculation is initialized by an open-shell SCF calculation specified in the usual way (either UHF or RHF). The resulting α and β density matrices from the SCF calculation are then added into a spatial density matrix and this is diagonalized to yield unrestricted natural orbitals (Pulay’s UNOs). These orbitals are used to create a new reference determinant. The UNOs can be used in conjunction with QRHF and NEWVRT keywords to create more complicated vacuum states for post-Hartree-Fock calculations. [NOTE: The UNO keywords together with MAKERHF can be used to create a reference determinant that is used primarily in connection to certain types of multireference calculations run in mrcc.]

## UNO_CHARGE

Sets the charge of the final reference state.

## UNO_MULT

Sets the spin-multiplicity of the final reference state. The orbitals are occupied in order of their natural occupation number. Often this keyword is used to create a closed shellreference (vital if the mrcc module is used). In this case, the program can proceed as an RHF calculation.

## MAKERHF

MAKERHF=ON instructs the post-SCF logic to proceed as an RHF calculation and must be specified in an mrcc calculation.

# Post-SCF file options

## SINGLE_STORE

Controls storing permutations of commonly resorted two-particle quantities like MO integrals and T amplitudes. If a calculation is running out of disk space, then setting SINGLE STORE=ON might allow the calculation to finish.

## ABCDTYPE

Specifies the way that VVVV molecular orbital integrals (having four virtual MOs) are handled in post-SCF calculations. =STANDARD uses a technique that minimizes CPU time but makes liberal use of disk storage (particularly during the integral processing program intprc). =MULTIPASS avoids creating intermediate sort files during intprc by reading the HF2 files (created by vtran) multiple times (note that this option requires more CPU processing time and the setting of HF2 FILE=SA VE). =AOBASIS uses an AO-based algorithm to evaluate all terms involving the VVVV MO integrals. Again, use of this option results in considerably longer CPU times but significantly reduces the amount of disk storage (more so for UHF and ROHF references than RHF).

## AO_LADDERS

Specifies the algorithm used by ACES II when ABCDTYPE=AOBASIS. =MULTIPASS uses an approach where the AO integral file is read a number of times to maximize vectorization and is usually the optimal strategy on vector supercomputers. =SINGLEPASS determines the contributions with only a single pass through the AO integral files at the cost of significantly reduced vectorization.

## SAVE_INTS

Controls deletion of the AO integrals file(s) after the AO to MO integral transformation. If SA VE INTS=ON, then these files are not deleted. Any method that requires AO integrals during a post-SCF calculation will automatically switch on this keyword.

## VTRAN

Specifies what type of integral transformation is to be performed in the transformation program vtran. =FULL/PARTIAL allows the program to choose the appropriate type of transformation, while =FULL forces a full integral transformation, and =PARTIAL skips the VVVV integrals. This keyword is set automatically based on other relevant keywords.

## HF2_FILE

Specifies whether the HF2 file series (including HF2AA, HF2BB, and HF2AB) is deleted after the first stage of integral processing. The default is to delete these files; however, when ABCDTYPE=MULTIPASS, these files must not be deleted and the program sets HF2 FILE=SA VE.

## ABCDFULL

This is a debug aid and is concerned with the storage of VVVV integrals and effective Hamiltonian elements. It is only relevant to RHF calculations, and it should never be necessary for users to set this keyword. Possible values are: =UNKNOWN (the program will choose an appropriate value), =ON (full storage), and =OFF (reduced storage).

## HBARABCD

Controls the formation (and storage) of the VVVV block of the effective Hamiltonian. It should never be necessary for users to set this keyword. Possible values are: =UNKNOWN (the program will choose an appropriate value), =ON (full storage), and =OFF (reduced storage).

## HBARABCI

Controls the formation (and storage) of the VVVO block of the effective Hamiltonian. It should never be necessary for users to set this keyword. Possible values are: =UNKNOWN (the program will choose an appropriate value), =ON (full storage), and =OFF (reduced storage).

## GAMMA_ABCD

Controls the evaluation of the two-particle density elements when all four indices correspond to virtual orbitals. This works in conjunction with the ABCDTYPE keyword. When ABCDTYPE=AOBASIS, the ABCD integrals are not stored on disk and GAMMA ABCD is set to DIRECT. When ABCD integrals are stored to disk, GAMMA ABCD must be set to DISK.

# Post-SCF calculations

## CC_CONV

Sets the convergence criterion for the CC and Lambda equations. Equations are considered converged when the maximum change in amplitudes is less than 10<sup>-N</sup> .

## CC_MAXCYC

Sets the maximum number of CC or Lambda iterations.

## CC_EXTRAPOL

(Formerly RLE.) Specifies the type of convergence acceleration used during the CC iterations. =STANDARD uses the RLE method of Bartlett and Purvis with periodic extrapolation of the solution vector, =DIIS uses the DIIS approach of Pulay, =NOJACOBI uses the RLE method with continuous extrapolation, and =OFF uses no convergence acceleration. In general, DIIS works well for a wide range of molecular systems while RLE can be better for certain cases. NOJACOBI might offer advantages in cases where the reduced subspace becomes singular too rapidly. NOJACOBI requires some additional disk storage, which might be disadvantageous for very large calculations. =OFF is generally a bad idea for CC calculations but might be preferred by some CI calculations.

## CC_EXPORDER

(Formerly ORDER_RLE.) Sets the maximum number of iterates to include in the R matrix used by RLE and DIIS. The maximum value allowed is 25.

## XFORM_TOL

Sets the tolerance for storing transformed integrals. Integrals less than 10<sup>-N</sup> are neglected and not stored on disk.

## NTOP_TAMP

Sets the N largest T amplitudes to print for each spin case and excitation level.

## TAMP_SUM

Sets how often the largest T amplitudes are printed. =0 only prints amplitudes at the beginning and end of the run, =1 prints amplitudes after every iteration, =2 prints amplitudes after every other iteration, etc.

## DENSITY

Specifies whether the relaxed density matrix is computed for correlated wavefunctions. =RESPONSE skips the orbital relaxation terms that contribute to the density matrix. This keyword should only be used by advanced users and ACES II developers who are testing new theories.

## DOHBAR

Controls what action is taken by the linear response program. =ON calculates and saves the full effective Hamiltonian. =OFF solves the “lambda” linear response equations.

# Excited states: general

## EXCITE

Selects the method for excited states of the same charge: `TDA` (CIS), `RPA`, `DRPA`, `EOMEE` (EOM-CCSD), `CIS(D)`, `P-EOMEE` (partitioned EOM-CCSD), `STEOM`, and EOM variants built on ring-CCD ground states (`EOM(S)-RCCD`, `EOM(S)-DRCCD` and their spin-flip and doubles forms, used with `EOMREF`). The number of states per irrep is set with `ESTATE_SYM`. See {doc}`../user/methods`.

## ESTATE_TOL

Sets the tolerance used in converging EOM-CC excited state calculations. The iterative diagonalization continues until the RMS residual falls below 10<sup>-N</sup> .

## ESTATE_MAXC

Sets the maximum number of iterative diagonalization steps for each root in excited state calculations.

## EOM_MAXCYC

Sets the number of iterations in EOM excited state calculations. If it has the value N and N ROOT roots have been requested for a given symmetry, then the program will allow up to N∗N ROOT iterations to find all of the requested roots for that symmetry. This keyword is only used by the iterative solution of the HF Stability analysis and has nothing to do with EOM.

## PROGRAM

## EOMREF

## IMEM_SIZE

## EOM_PRJCT

## NT3EOMEE

## ACC_SYM

# Excited states: properties

## ESTATE_PROP

Specifies whether any excited state one-electron properties are to be calculated, but it only applies to EOM-CC calculations. Proper use of this keyword might require some fairly advanced knowledge of quantum chemistry and the available options are discussed here. The values are: =OFF→ no properties or transition moments are calculated =EXPECTATION→ transition moments and dipole strengths are calculated along with selected one-electron properties which are evaluated as expectation values =UNRELAXED→ selected one-electron properties are calculated in an approximation that neglects relaxation of molecular orbitals =RESPONSE→ selected one-electron properties are calculated as analytic first derivatives of the energy Except for EOMCC calculations on two-electron systems (which are exact), properties obtained by the three approaches will not be equivalent. The default value for this keyword is slightly complicated. For TDA calculations, the default is EXPECTATION since the evaluation of transition moments involves only a negligible amount of additional computation relative to the evaluation of the excitation energies. For EOMCC, the default is OFF since evaluation of any transition moments or properties requires approximately twice the computational time. ESTATE PROP=RESPONSE is not available for EOMCC calculations. Transition moments and dipole strengths are evaluated by default for all values of ESTATE PROP other than OFF.

# Excited states: affinities

## EA_CALC

Specifies the method for calculating electron affinities of a closed-shell parent state. The values are: =NONE, =MBPT(2) (strict second-order perturbation theory), =SO DYSON (second-order Green’s Function, iterating the MBPT(2) formula), =OVGF, =P EOMEA (partitioned EA-EOM, the 2p1h-2p1h block is treated as diagonal), =EA EOMCC (Hbar is diagonalized over 1p and 2p1h configurations), =COMBO, =OS CCSD (a stateselective multireference CC method, in which both orbitals and cluster amplitudes are optimized for the final state of interest). This last method requires additional input in the mrcc gen namelist. Analytical gradients are available for EA-EOMCC and PEOMEA, and require additional input though the mrcc gen namelist. Properties and transition properties of the EA states can be requested by additional input in mrcc namelists ( *ea calc section).

## EA_SYM

Specifies the number of states to be calculated by the EA CALC calculation. A string (e.g., 4-2-2-0) has to be provided that indicates the number of doublet states in each symmetry block. The example string requests 4 doublet states of A1 symmetry, 2 doublet states in symmetry blocks 2 and 3, and 0 in block 4. The number of EA states to be calculated can also be specified by listing an energy range using ea low and ea high keywords in the mrcc gen namelist. This can be particularly relevant in vibrational frequency calculations in which the symmetry changes at different geometries.

## DEA_CALC

Specifies the method for calculating double electron affinities after calculation of the closed-shell parent state. This type of calculation can be very useful to describe certain multireference situations like the C atom or the N 2O2 molecule. The values are: =NONE, =TDA (the analogue of CIS, in which a CI calculation is performed over 2-particle states), =EOMCC, =STEOM (Similarity Transformed Equation of Motion), =OS CCSD, =SS STEOM. DEA-STEOM is only meaningful if EA CALC=EA EOMCC. OS-CCSD and SS-STEOM are multireference CC methods, which are still in an experimental stage. Analytical gradients are available for TDA and STEOM and require additional input though the mrcc gen namelist. Properties and transition properties of the DEA states can be requested by additional input in mrcc namelists (*dea calc section).

## DEA_SYM

Specifies the number of states to be calculated by a DEA calculation. A string (e.g., 4-2-2-0/2-1-1-1) has to be provided that indicates the number of singlet states in each symmetry block, followed by the number of triplet states in each symmetry block. The example string requests 4 singlet states of A1 symmetry, 2 singlet states in symmetry blocks 2 and 3, and 0 in block 4. In addition, 2 triplet states will be calculated in block 1 and 1 triplet in blocks 2, 3, and 4 each. The triplet vector and the forward-slash are optional.

# Excited states: electronic (absorption)

## EE_SYM

Specifies the number of states to be calculated by an EE calculation. A string (e.g., 4-2-2-0) has to be provided that indicates the number of singlet states in each symmetry block. The example string requests 4 singlet states of A1 symmetry, 2 states in symmetry blocks 2 and 3, and 0 in block 4. The number of EE states to calculate can also be specified by listing an energy range using ee low and ee high keywords in the mrcc gen namelist. This can be particularly relevant in vibrational frequency calculations in which the symmetry changes at different geometries. RHF is limited to singlets, but UHF can do singlets and triplets using the N-N/N-N notation like DEA SYM.

## EE_SEARCH

Specifies the character of the states to calculate. If EE SYM is not specified, the program attempts to determine all states of the given character; otherwise, it uses the symmetry constraints imposed by EE SYM. The values are: =LOWEST (find the lowest energy roots regardless of character), =CORE (find excitations from core orbitals using guess vectors from a projected TDA matrix), =LUMO (find excitations to the LUMO), and =HOMO (find excitations to the HOMO).

## ESTATE_SYM

Specifies the number of excited states which are to be determined in each irreducible representation of the computational subgroup. The program attempts to find all of the lowest roots, but this is not guaranteed since the eigenvalue problem is not solved by direct matrix diagonalization, rather by an iterative (modified Davidson) algorithm. For excited state gradient calculations (TDA only), only one root can be specified, so only one non-zero entry in the string is allowed and that entry must be set to one. The format used for this keyword is identical to that used in the OCCUPATION keyword. For example, for a computational subgroup having four symmetry species, the string ESTATE SYM=3-1-0-2 specifies that 6 total roots should be searched for, three in the first block, one in the second block, and two in the fourth block. This keyword is the same as EE SYM, but EE SYM is read when PROGRAM=ACES3 and ESTATE SYM is read when PROGRAM=ACES2.

# Excited states: ionizations

## IP_CALC

Specifies the method for calculating ionization potentials of the closed-shell parent state. The values are: =NONE, =MBPT(2) (strict second-order perturbation theory), =SO DYSON (second-order Green’s Function, iterating the MBPT(2) formula), =OVGF, =P EOMIP (partitioned IP-EOM, the 2h1p-2h1p block is treated as diagonal), =IP EOMCC (H-bar is diagonalized over 1h and 2h1p configurations), =COMBO, and =OS CCSD (a state-selective multireference CC method, in which both orbitals and cluster amplitudes are optimized for the final ionized state of interest). This last method requires additional input in the mrcc gen namelist. Analytical gradients are available for IP-EOMCC and P-EOMIP and require additional input though the mrcc gen namelist. Properties and transition properties of the IP states can be requested by additional input in mrcc namelists ( *ip calc section).

## IP_SYM

Specifies the number of states to calculate by an IP calculation. A string (e.g., 4-2-20) has to be provided that indicates the number of doublet states in each symmetry block. The example string requests 4 doublet states of A1 symmetry, 2 doublet states in symmetry blocks 2 and 3, and 0 in block 4. The number of IP states to be calculated can also be specified by listing an energy range using ip low and ip high keywords in the mrcc gen namelist. This can be particularly relevant in vibrational frequency calculations, in which the symmetry changes at different geometries.

## IP_SEARCH

Specifies the character of the IP states to zoom in on. If IP SYM is not specified, the program attempts to determine all IP’s of the given character; otherwise, it uses the symmetry constraints imposed by IP SYM. The values are: =V ALENCE (all valence IP’s), =LOWEST (only the ground state of the cation), =COREIP (only the core ionization potentials—IPs larger than about 100 eV), =SHAKEUP (shake-up states), and =KOOPMANS (all principle IP’s having predominantly 1h character). KOOPMANS type can be very difficult for inner valence ionization potentials, and SHAKEUP is not implemented yet.

## DIP_CALC

Specifies the method for calculating double IP states of the closed-shell parent state. This type of calculation can be very useful to describe certain multireference situations like biradicals. The values are: =NONE, =TDA (the analogue of CIS, in which a CI calculation is performed over 2-hole states), =EOMCC (H-bar is diagonalized over 2h and 3h1p calculations), =STEOM (Similarity Transformed Equation of Motion), =OS CCSD, and =SS STEOM. To run a meaningful DIP-STEOM calculation one must also set IP CALC=IP EOMCC. OS-CCSD and SS-STEOM are multireference CC methods, which are still in an experimental stage. Analytical gradients are available for TDA and STEOM and require additional input though the mrcc gen namelist. Properties and transition properties of the DIP states can be requested by additional input in mrcc namelists ( *dip calc section).

## DIP_SYM

Specifies the number of states to be calculated by a DIP calculation. A string (e.g., 4-2-2-0/2-1-1-1) has to be provided that indicates the number of singlet states in each symmetry block, followed by the number of triplet states in each symmetry block. The example string requests 4 singlet states of A1 symmetry, 2 singlet states in symmetry blocks 2 and 3, and 0 in block 4. In addition 2 triplet states will be calculated in block 1 and 1 in blocks 2, 3, and 4 each. The triplet vector and the forward-slash are optional.

# Excited states: gradients

## ZETA_TYPE

Specifies the algorithm used to solve linear equations (zeta equations) in excited state gradient theory. =POPLE uses Pople’s orthogonal subspace approach, and =DIIS uses Pulay’s DIIS approach.

## ZETA_CONV

Sets the convergence criterion for the iterative solution of the Zeta equations and Zvector equations. The solutions are considered to be converged when the error falls below 10<sup>-N</sup> .

## ZETA_MAXCYC

Sets the maximum number of cycles allowed for the solution of the Zeta equations.

## RESRAMAN

# Properties

## PROPS

Specifies whether properties other than the energy, geometry, and vibrational frequencies are calculated. Acceptable values are: =OFF→ no properties are calculated =FIRST ORDER→ dipole moment, quadrupole moment, electrical field gradients, spin densities, etc. and the approximate perturbative relativistic correction to the energy by Cowan-Griffin =SECOND ORDER→ frequency-independent polarizabilities, commonly known as the CPHF polarizabilities, provided CALC=SCF =EOM NLO→ both frequency-dependent and -independent polarizabilities at the CCSD level using EOM-CC (see EOMPROP) =J SO→ paramagnetic and diamagnetic spin-orbit contributions to the NMR spin coupling constant =J SD→ spin-dipole contribution to the NMR spin coupling constant =J FC→ Fermi-contact contribution to the NMR spin coupling constant =JSC ALL→ all four contributions =NMR→ NMR chemical shifts (limited to SCF and MBPT(2) with SPHERICAL=OFF) Note that options J SO, J SD, J FC, and JSC ALL require CALC=CCSD and that options J SO, J FC, and JSC ALL require REF=UHF . (see EOMPROP)

## EOMPROP

Specifies the method of calculating EOM-CCSD second-order properties (polarizabilities and spin-spin coupling constants). =CILIKE uses a CI-like formula, which is not rigorously size-extensive, =LINEAR removes unlinked diagrams from the CI-like formula, =QUADRATIC uses a size-extensive quadratic formula, and =COMBO computes all approximations. This keyword only has meaning when PROP =EOM NLO, =J SO, =J FC, =J SD, or =JSC ALL.

## TDHF

Controls whether a time-dependent Hartree-Fock calculation of nonlinear optical properties is to be performed. This keyword can only be used for closed-shell SCF calculations with no dropped MOs. The nonlinear properties which are to be calculated are specified by a namelist which is put at the end of the ZMAT file. A description of this namelist and an example can be found in Section 9.1.10 (page 77).

## CPHF_CONVER

Sets the convergence criterion for the iterative solution of the CPHF and Z-vector equations. The solutions are considered to be converged when the error falls below 10-N .

## CPHF_MAXCYC

Sets the maximum number of cycles allowed for the solution of the CPHF and Z-vector equations.

## TREAT_PERT

This keyword is used for certain types of correlated second-derivative calculations [presently only GIAO NMR shift calculations] and directs ACES II to either treat all perturbations at once or treat them sequentially. The latter approach results in less demand for physical disk space at the cost of increased cpu time. Available options are SIMULTANEOUS and SEQUENTIAL.

## XFIELD

Sets the X-component of an external electric field. The value must be specified as an integer and the field used by the program will be N∗ 10-6. This allows field strengths |E| > 10-6 to be used.

## YFIELD

Sets the Y-component of an external electric field. See above.

## ZFIELD

Sets the Z-component of an external electric field. See above.

## PRP_INTS

Specifies the types of property integrals that are computed. Setting this to FULL or PARTIAL computes the full set or sub set, respectively. This works in conjunction with the PROPS keyword and the defaults are set automatically.

# Geometry optimization: general

## GEOM_OPT

Specifies the scope of coordinate optimizations. This keyword is automatically set to PARTIAL if a geometry optimization is implied with asterisks in the internal coordinate Z matrix. Setting GEOM OPT=FULL will optimize all coordinates regardless of asterisks in the Z matrix (COORD=INTERNAL) and will optimize a Cartesian input geometry using Redundant Internal Coordinates (RICs).

## CURVILINEAR

Specifies whether or not the Hessian matrix is transformed (nonlinearly) to curvilinear internal coordinates. =OFF turns the transformation off if analytical force constants are not available, while it is always performed if CUR VILINEAR=ON. =NO unconditionally turns the transformation off. This keyword is set automatically.

# Geometry optimization: stepping algorithm

## OPT_METHOD

Specifies the geometry optimization strategy. =AUTO, =NR uses a straightforward Newton-Raphson minimum energy search, =RF A uses the Rational Function Approximation minimum energy search (and can be used when the initial structure is in a region where the number of negative Hessian eigenvalues is nonzero), =MANR uses a Morse-adjusted Newton-Raphson minimum energy search (and is very efficient if the Hessian is available), =EVFTS uses Cerjan-Miller eigenvector following for finding a transition state (can be started in a region where the Hessian index, the number of negative Hessian eigenvalues, is not equal to one), =MAEVFTS uses a Morse-adjusted eigenvector following for finding a transition state.

## MAX_STEP

Sets the largest step size in millibohr.

## STP_SIZ_CTL

Controls how the step size is scaled. =TRUST RADIUS uses the dynamic scaling by Fletcher, =NORM uses the absolute step length, and =MAXIMUM uses the largest individual step in the internal coordinate space.

## EIGENVECTOR

Sets which eigenvector of the totally symmetric part of the block-factored Hessian is to be followed uphill in a transition state search. Eigenvectors are indexed by their eigenvalues (the lowest eigenvalue is 1, the next lowest is 2, etc.). The default should always be used if you are not looking for a specific transition state which you know corresponds to motion along a different mode. The value of EIGENVECTOR has no meaning if OPT METHOD is not set to EVFTS or MAEVFTS.

## NEGEVAL

Specifies the behavior if negative eigenvalues are encountered in the totally symmetric Hessian during an NR or MANR search. If NEGEV AL=ABORT, then the job will terminate with an error message. If NEGEV AL=ABSV AL, then the program will just switch the eigenvalue to its absolute value and keep plugging away (this is strongly discouraged). If NEGEV AL=RF A, then OPTMETHOD is switched to RF A internally and the optimization is continued.

# Geometry optimization: iteration control

## CONVERGENCE

Sets the convergence criterion for geometry optimizations. The optimization terminates when the RMS gradient is below 10<sup>-N</sup> Hartree/Bohr.

## OPT_MAXCYC

Sets the maximum number of optimization cycles.

## INIT_HESSIAN

Specifies the initial Hessian for geometry optimizations. =SPECIAL generates an approximate Hessian internally, =FCMINT imports an approximate Hessian from the FCMINT file, =MOPAC generates the FCMINT file by running mopac, and =EXACT uses analytical second derivatives to generate an exact Hessian (limited to SCF).

## HESS_UPDATE

Specifies the algorithm used to update the Hessian. =NONE uses the same Hessian throughout the optimization, =POWELL uses the Powell update, =BFGS uses the Broyden-Fletcher-Goldfarb-Shanno update, =MS uses the Murtagh-Sargent update, =BOFILL uses the Bofill mixture of Powell and Murtagh-Saregent, and =PSB uses the Powell symmetric Broyden update. =POWELL and =BOFILL are the defaults for minimum energy and transition state searches, respectively.

## EVAL_HESS

Sets the cycle interval for recomputing the Hessian. For correlated calculations, the Hessian is evaluated only at the SCF level .

# Geometry optimization: integral derivatives

## TRANS_INV

Specifies whether or not translational invariance is exploited in derivative calculations. =USE uses translational invariance, while =IGNORE does not.

# Frequencies and other 2nd-order properties

## VIBRATION

Specifies the method for calculating vibrational frequencies. =EXACT (SCF only) performs normal mode analysis on an analytic force constant matrix and computes rotationally projected frequencies and infrared intensities. =FINDIF signals ACES II to compute the force constant matrix by finite difference of analytically computed gradients or energies using symmetry-adapted mass-weighted Cartesian coordinates.

## RAMAN

Controls raman intensities and depolarization ratios of numerical vibrational frequency calculations (VIB=FINDIF). The Raman intensities are limited to the SCF and CCSD level (via EOM-CCSD), and as a result RAMAN keyword must be used in conjunction with CALC=SCF or CALC=CCSD. We recommend setting EOMPROP=QUADRATIC for raman intensities at the CCSD level. [NOTE: This method requires an extra input file named frequency, which is described in Section 6 (page 22).] (for numerical gradients and Hessians)

# Finite displacements

## FD_STEPSIZE

Sets the step length (in 10<sup>-4</sup> amu1/2 Bohr) used in generating the force constant matrix by finite difference of Cartesian gradients. (The default is 0.005 amu 1/2 Bohr.)

## FD_IRREPS

Lists the symmetry types to be evaluated in a VIBRATION=FINDIF calculation. The numbers of the irreducible representations for which vibrational analysis is to be performed are separated by forward-slashes. For example, FD IRREP=1/3/4 computes the frequencies of modes transforming as the first, third, and fourth irreducible representations. If a symmetry is specified for which there are no vibrational modes, the program will terminate. The labels of the irreducible representations for this keyword are not usually the same as those used in the symmetry-adapted integrals. Moreover, for some point groups (like those of linear molecules), the two sets of labels refer to different subgroups. There is still no straightforward way to determine what they will be without starting a calculation. One joda run will list the relevant irreducible representations. If all vibrational frequencies are desired, this keyword need not be included. The default is to calculate all irreps.

## FD_PROJECT

Controls whether or not rotational degrees of freedom are projected out of symmetryadapted coordinates. =ON uses rotationally projected coordinates, while =OFF retains the rotational degrees of freedom. At a stationary point on the potential surface, both options will give equivalent harmonic force fields, but =OFF should be used at non-stationary points.

## FD_USEGROUP

Specifies the point group to be used in generating the symmetry-adapted vibrational coordinates. =FULL specifies the full molecular point group, =COMP specifies the Abelian subgroup used in the electronic structure calculation.

## POINTS

Specifies either single- (=SINGLE) or double-sided (=DOUBLE) numerical differentiation in the finite difference evaluation of the Hessian. Two-sided numerical differentiation is considerably more accurate than the single-sided method, and its use is strongly recommended for production work.

# External interfaces

## EXTERNAL

Specifies what type of external file interface for third-party programs to create. Current interfaces include =HYPERCHEM and =MOLDEN. a2proc help This is not a keyword per se but is used by ACES II to generate the external interfaces, and it can be used to generate such interfaces after a calculation has finished.

# Experimental, obsolete, and unused

## SOLVENT

Sets the dielectric constant used to determine the orbitals. A cavity size may be specified as well by creating a file named “radius”, which is read by the SCF code. This contains the cavity radius in Å. If the file is not present, the program uses a value calculated from 0.5 Å plus half of the longest internuclear distance.

## TURBOMOLE

## POLYRATE

## CCR12

## KUCHARSKI

## KS_POT

Obsoleted by the *VSCF namelist.

## FUNCTIONAL

Obsoleted by the *INTGRT namelist.

## EOMXFIELD

## EOMYFIELD

## EOMZFIELD

## RDO

## FILE_STRIPE

## RESET_FLAGS

## PSI

## INSERTF

## GLOBAL_MEM

# Added since 2006

## TRUNC_ORBS

Truncates the virtual space before the correlated calculation. `FNO` replaces the virtual orbitals by frozen natural orbitals from a preliminary MBPT(2) calculation and keeps the fraction given by `FNO_KEEP`; `OVOS` uses optimized virtual orbitals. See test `124a.fno`.

## FNO_KEEP

With `TRUNC_ORBS=FNO`, the percentage of virtual orbitals kept (for example `FNO_KEEP=60`).

## DAMP_CC

Damping of the coupled-cluster amplitude update, in units of 0.01: the default 100 is a full (undamped) update, smaller values mix in the previous amplitudes. Can help difficult CC convergence. Used by `xvcc`, the pCCD and ring-CC codes.

## DKH_ORDER

Order of the Douglas–Kroll–Hess scalar-relativistic Hamiltonian; 0 (the default) turns it off.

## UNCONTRACT

With `DKH_ORDER` > 0, builds the DKH one-electron operator in the uncontracted basis instead of the contracted one.

## FFIELD

Field strength, in units of 10<sup>-6</sup> a.u., for finite-field calculations of polarizabilities (see `ORDER_ALPHA`). A non-zero value makes `xaces2` run the finite-field sequence; a positive value requires `SYMMETRY=OFF`.

## ORDER_ALPHA

For finite-field calculations (`FFIELD`), the multipole operator whose response is computed: `DIPOLE`, `QUADRUPOLE` or `OCTUPOLE`.

## OO_METHOD

Orbital optimizer for orbital-optimized pCCD and related methods (`PCCD_REF=OPTIMIZED`, `CALC=OO-CCD`, `OO-LCCD`): `RFO` (rational function), `BFGS`, `FULL-NR` (Newton–Raphson), `BB-SD`, `AMSGRAD`, `ROTGRAD`. See {doc}`../user/methods`.

## PCCD_REF

Orbitals used by pCCD: `HF` (canonical SCF), `KS`, `BRUECKNER`, `NATURAL`, `LOCAL`, or `OPTIMIZED` for orbital-optimized pCCD.

## OO_MAXCYC

Maximum number of orbital-optimization steps in OO-pCCD and related methods.

## OO_DIIS

DIIS acceleration in the pCCD amplitude and orbital iterations (on by default).

## RPA_SOLVE

How RPA-type eigenvalue problems are solved: `FULL_DIAG` (full diagonalization, the default) or `DAVIDSON` (iterative, for the lowest roots only).

## EE_ROOT_FIND

Eigenvalue solver for EOM-CC excitation energies in `xvee`: `DAVIDSON` (the default), `BLOCK_DAVIDSON`, `LANCZOS`, or `FULL_DIAG`.

## IP_ROOT_FIND

Eigenvalue solver for EOM-CC ionization potentials in `xvip`: `DAVIDSON` (the default) or `FULL_DIAG`.

## EOM_SOLVE

For EOM-IP in `xvip`, whether the right (`RIGHT`, the default) or left (`LEFT`) eigenvectors are solved for.

