











      subroutine numintAG(dograd)

      implicit none







c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler























cYAU - ACES3 stuff . . . we hope - #include <aces.par>






      integer           iuhf
      common /iuhf_com/ iuhf
      save   /iuhf_com/




c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end





c This contains flags that are set in the INTGRT namelist.  See the
c file initintgrt.F for a description of each of them.
c
c The following are exceptions:
c    int_ks           : .true. if we are doing Kohn-Sham
c    int_ks_finaliter : .true. if this is the final iteration of KS
c    int_ks_exch      : which potential to use to calculate exchange
c    int_ks_corr      : which potential to use to calculate correlation
c    int_kspot        : which hybrid functional to use to calculate potential
c                       (if equal to fun_special, use int_ks_exch and
c                       int_ks_corr)
c    int_dft_fun      : which functional to use with any SCF density
c                       Added and modified by Stan Ivanov
c                       (if fun_special the functional is user-defined
c                        if fun_hyb_name then use hybrid functional) 
c    int_printlev     : 0 if we are doing a dft calculation, 1 if we
c                       are doing the final iteration of a KS calculation,
c                       2 if we are doing a KS iteration.
c These are set in the calling routines, NOT in the namelist.
c                     : Additions by S. Ivanov
c     num_acc_ks      : .true.  if numerical accelerator is used for KS 
c                        Default is .true.
c     ks_exact_ex     : .true. if exact LOCAL exchange is used for KS
c                        Deafult is .false.
c     int_tdks        : .true. if time-dependent KS calculation is
c                        requested
c                        Default is .false.
c     int_ks_scf      : .true. if the actual KS SCF energy is being
c                        calculated and printed out. Default is .false.
c
      integer int_numradpts,int_radtyp,int_partpoly,int_radscal,
     &    int_parttyp,int_fuzzyiter,int_defenegrid,int_defenetype,
     &    int_defpotgrid,int_defpottype,int_kspot,
     &    int_ks_exch,int_ks_corr,int_dft_fun,

     &    int_printlev,
     &    int_printscf,int_printint,int_printsize,int_printatom,
     &    int_printmos,int_printocc,
     &    potradpts, numauxbas,int_ksmem,int_overlp

      logical int_ks,num_acc_ks,ks_exact_ex,int_tdks,int_ks_scf,
     &        int_ks_finaliter

      double precision
     &    int_radlimit,coef_pot_nonlocal

      common /intgrtflags/  int_numradpts,int_radtyp,int_partpoly,
     &    int_radscal,int_parttyp,int_fuzzyiter,int_defenegrid,
     &    int_defenetype,int_defpotgrid,int_defpottype,int_kspot,
     &    int_ks_exch,int_ks_corr,int_dft_fun,

     &    int_printlev,
     &    int_printscf,int_printint,int_printsize,int_printatom,
     &    int_ks_finaliter,int_printmos,int_printocc,
     &    potradpts, numauxbas,int_ksmem,int_overlp
c

c  prakash added int_ksmem to the common block
      common /intgrtflagsd/ int_radlimit,coef_pot_nonlocal
      common /intgrtflagsl/ int_ks,num_acc_ks,ks_exact_ex,int_tdks,
     &                      int_ks_scf

      save /intgrtflags/
      save /intgrtflagsl/
      save /intgrtflagsd/

c The following are parameters used in the namelist

      integer int_prt_never,int_prt_dft,int_prt_ks,int_prt_always
      parameter (int_prt_never   =1)
      parameter (int_prt_dft     =2)
      parameter (int_prt_ks      =3)
      parameter (int_prt_always  =4)

      integer int_radtyp_handy,int_radtyp_gl
      parameter (int_radtyp_handy=1)
      parameter (int_radtyp_gl   =2)

      integer int_partpoly_equal,int_partpoly_bsrad,
     &    int_partpoly_dynamic
      parameter (int_partpoly_equal  =1)
      parameter (int_partpoly_bsrad  =2)
      parameter (int_partpoly_dynamic=3)

      integer int_radscal_none,int_radscal_slater
      parameter (int_radscal_none  =1)
      parameter (int_radscal_slater=2)

      integer int_parttyp_rigid,int_parttyp_fuzzy
      parameter (int_parttyp_rigid=1)
      parameter (int_parttyp_fuzzy=2)

      integer int_gridtype_leb
      parameter (int_gridtype_leb=1)


c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end


c Two different "stacks" of memory will be available in the course of an
c Aces3 calculation.  Integer memory will be stored in icore.  Real memory
c will be stored in dcore.
c
c Memory allocation is done in one of two ways, dynamic or nondynamic.  This
c is controlled by the parameter dynmem.  If dynmem is 1, dynamic memory
c is used.  Otherwise, nondynamic memory is used.
c
c Two additional parameters nondynimem and nondyndmem control how much
c nondynamic memory is allocated initially in the icore and dcore stacks.
c If dynmem is 1, both of these parameters SHOULD be 1.

c Dynamic memory allocation is done in a fairly straighforward way.  The
c program runs through everything twice.  The first time is just to determine
c how much dynamic memory is required.  The second time is use to actually
c allocate the memory and use it.  In some cases, it may be difficult to
c determine how much memory to use in advance in each of the stacks.  To
c aid this, there is also the option of only running through a program a
c single time.  When this happens, both stacks are allocated initially and
c used throughout the program.  If the memory requirement for either stack
c is exceeded, the program crashes.
c
c Two parameters which are used in either case are memknown and maxicore.
c Both MUST be set in the calling program.
c
c   memknown : 0 if a dummy run is being done to determine memory usage
c              1 if the memory usage is known from a previous dummy run
c             -1 if no dummy run is done (only one run through the program)
c   maxicore : the maximum amount of icore to allocate (this MAY be
c              adjusted IF a dummy run detects that more is needed, but in
c              the case of a single run, this will be rigidly applied)
c
c This include file is required whenever calls to any of the memory functions
c (setptr and relptr) are called or when part of the execution depends on
c whether the first or second run is being performed.
c
c ***NOTE***
c It is important that icore be aligned on a floating point boundary.  One
c way which seems to insure this is to have icore the FIRST element in the
c common block.  So, make sure that no integer is ever inserted before icore
c in the common block declaration.

      integer dynmem,nondynimem,nondyndmem
      parameter (dynmem=1)
      parameter (nondynimem=1)
      parameter (nondyndmem=1)

      double precision   dcore(nondyndmem)
      common /dcore_com/ dcore
      save   /dcore_com/

      integer    kscore(nondynimem)
      common /kscore_mem / kscore

      integer            memknown,maxicore
      common /kscore_com/ memknown,maxicore
      save   /kscore_com/





c This common block contains the molecule and basis set information
c read in from the JOBARC file.  Since much of this information is
c used in a large number of modules, and since most of the information
c is relatively small compared to the other things held in memory,
c a large percentage of the data stored in the JOBARC file is stored
c here, even though some modules will not use all of it.

c   maxangshell - The maximum number of angular momentum shells.  Since this
c                 is used VERY infrequently, set it high enough to never
c                 cause a problem.
c   spinc(2)    - The characters 'A' and 'B' (useful for alpha/beta labels)
c   natoms      - The number of atoms in the Z-matrix (including X/GH).  After
c                 remove is called, natoms becomes equivalent to nrealatm.
c   natomsx     - The number of atoms in the Z-matrix (including X/GH).  This
c                 does not change.
c   nrealatm    - The number of atoms in the Z-matrix (including GH).
c   naobasfn    - The number of AOs in the basis
c   nbastot     - The number of symmetry adapted orbitals in the basis (the AO
c                 basis may be larger than the SO basis if spherical orbitals
c                 are used since harmonic contaminants are deleted)
c   linear      - 1 if the molecule is linear
c   orientmt    - 3x3 matrix which relates the computational and canonical
c                 orientations
c   nucrep      - Nuclear repulsion energy in a.u.
c   nmproton    - Number of protons in the molecule.
c
c   compptgp    - Point group
c   fullptgp    -
c   compordr    - Order of the point group
c   fullordr    -
c   compnirr    - Number of irreps in the point group
c   fullnirr    -
c   compnorb    - Number of unique atoms (orbits) in the point group
c   fullnorb    -
c   c1symmet    - 1 if the molecule is C1 symmetry
c   nirrep      - The same as compnirr (since nirrep is used so commonly,
c                 this is included for conveniance)
c                 ***NOTE*** nirrep is read in twice and is stored in /sym/
c                            so it is not actually included here
c   totprim     - Total number of primitive functions in the molecule
c   maxshlprim  - Largest number of primitives in a single shell
c   maxshlao    - Largest number of AOs in a single shell
c   maxshlorb   - Largest number of primitive orbitals (primitive functions
c                 times the number of AOs) in a single shell
c   maxangmom   - Largest angular momentum for any atom
c   maxshell    - Larges number of angular momentum shells for any atom
c   noccorb(2)  - The number of alpha and beta occupied orbitals
c   nvrtorb(2)  - The number of alpha and beta virtual orbitals

c The parameter maxorbit is needed because of how dynamic memory is used.
c Two runs of the program are needed.  The first to calculate memory usage,
c the second to use it.  In order to calculate totprim, we have to know the
c orbit population vector (the number of each type of atom).  BUT, this is
c stored in dynamic memory since we do not know how long this vector is.
c In the future, joda or vmol will write this information to JOBARC, and
c this problem will disappear.  In the meantime, we have to introduce a
c genuine limit on the size of the molecule.  It may have no more than
c maxorbit sets of unique atoms.  This limit is ONLY used in the subroutine
c basis, so it probably will disappear when the information in the MOL file
c is put in JOBARC.
c    maxorbit   - the number of symmetry unique atoms

c The following are pointers to real arrays
c
c   zatommass(natoms)  - Atomic mass of all atoms (X=0.0, GH=100.0)
c   zcoord(3,natoms)   - Coordinates of all atoms (computational orientation)
c   zalpha(totprim)    - The alpha for each primitive function
c   zprimcoef(totprim,naobasfn)
c                      - The primitive to AO coefficients
c
c The following are pointers to integer arrays
c
c   patomchrg(natoms)  - Atomic number of all atoms (X=0, GH=110)
c   pfullclss(fullordr)- Class type vector
c   pcompclss(compordr)-
c   pfullpopv(natoms)  - Number of atoms in each orbit
c   pcomppopv(natoms)  -
c   pfullmemb(natoms)  - Atoms sorted by point group orbits
c   pcompmemb(natoms)  -
c   pnprimatom(natoms) - Number of primitive functions for each atom
c   pnshellatom(natoms)- Number of different angular momentum shells for each
c                        atom (takes on values of 1,4,9,16, etc.)
c   pnangatom(natoms)  - The number of different angular momentum for each
c                        atom (takes on values of 1,2,3,4, etc.)
c   pnaoatom(natoms)   - Number of AOs for each atom
c   pnshellprim(maxshell,natoms)
c                      - The number of primitive functions in each shell
c                        of each atom
c   pnshellao(maxshell,natoms)
c                      - The number of AOs in each shell of each atom
c   pprimoff(maxshell,natoms)
c   paooff(maxshell,natoms)
c                      - The primcoef matrix is a block diagonal matrix.
c                        Each shell of each atom has a block.  If you have
c                        a list of all primitive functions, pprimoff(ishell,
c                        iatom) tells the location of the first primitive
c                        function in the block (ishell,iatom) and paooff
c                        contains similar information for the AOs.
c
c ***NOTE***  Because joda stores pfullpopv/pcomppopv as size natoms, we
c             do to, but they should be of size fullnorb/compnorb.  The
c             first ones have real values.  The remaining ones are 0.

      double precision orientmt(3,3),nucrep
      integer natoms,nrealatm,naobasfn,nbastot,linear,compnirr,
     &    fullnirr,compnorb,fullnorb,compordr,fullordr,nmproton,
     &    c1symmet,totprim,maxshlprim,maxshlorb,maxshell,noccorb(2),
     &    nvrtorb(2),maxshlao,maxangmom,natomsx
      integer patomchrg,zatommass,zcoord,pfullclss,pcompclss,
     &    pfullpopv,pcomppopv,pfullmemb,pcompmemb,pnprimatom,
     &    pnshellatom,pnaoatom,pnshellprim,pnshellao,
     &    zalpha,zprimcoef,pprimoff,paooff,pnangatom

      common /mol_com/ orientmt,nucrep,
     &    natoms,nrealatm,naobasfn,nbastot,linear,compnirr,
     &    fullnirr,compnorb,fullnorb,compordr,fullordr,nmproton,
     &    c1symmet,totprim,maxshlprim,maxshlorb,maxshell,noccorb,
     &    nvrtorb,maxshlao,maxangmom,natomsx,
     &    patomchrg,zatommass,zcoord,pfullclss,pcompclss,
     &    pfullpopv,pcomppopv,pfullmemb,pcompmemb,pnprimatom,
     &    pnshellatom,pnaoatom,pnshellprim,pnshellao,zalpha,
     &    zprimcoef,pprimoff,paooff,pnangatom
      save   /mol_com/

      character*4 compptgp,fullptgp
      character*1 spinc(2)
      common /molc_com/ compptgp,fullptgp,spinc
      save   /molc_com/




c Analytic values determined in the anlytc subroutine

c All HF values are determined with the HF density
c All correlated values are determined using the correlated relaxed density

c ehfkin   : Hartree-Fock kinetic energy
c eckin    : correlated
c ehfnatr  : Hartree-Fock nuclear-electron attraction energy
c ecnatr   : correlated
c ehfcoul  : Hartree-Fock coulomb energy
c eccoul   : correlated
c ehfx     : Hartree-Fock exchange energy
c ecx      : correlated (using the correlated relaxed density in the
c            Hartree-Fock exchange energy expression)

c ehar     : Hartree energy using HF orbitals
c echar    : Hartree energy using natural orbitals and occupations
c escf     : SCF energy
c ecor     : energy determined using the relaxed density in the HF
c              energy expression

      double precision
     &    ehfkin,eckin,ehfnatr,ecnatr,ehfcoul,eccoul,ehfx,ecx,
     &    ehar,echar,escf,ecor
      common /exch/ ehfkin,eckin,ehfnatr,ecnatr,ehfcoul,eccoul,
     &    ehfx,ecx,ehar,echar,escf,ecor
      save /exch/




c This contains information about each of the possible grids for
c performing the numerical integration.

c###########################################################################
c MISC
c###########################################################################
c maxgrdatm : The largest atomic number for which the Slater and Bragg-Slater
c             atomic size is known.
c atmrad    : Atomic size using Slater's' rules for the radial integration
c xbsl      : The Bragg-Slater radii (one for each atom)
c             ***NOTE*** This is a genuine constraint.  Only atoms smaller
c             then this (currently 86) may be calculated.
c numangfct : number of angular momentum functions
c minpt     : The number of points used in the integration over
c             interatomic paths to find the minimum density point
c             between two atoms
c pangfct   : a pointer to the array of x,y,z angular momentum for each
c             angular momentum function

      integer maxgrdatm
      parameter (maxgrdatm=86)
      integer minpt
      parameter (minpt=100)

      integer numangfct

      double precision
     &    atmrad(maxgrdatm),xbsl(maxgrdatm),
     &     TA(maxgrdatm),multiEX(maxgrdatm)
      common /grid/  numangfct
      save /grid/

      common /gridd/ atmrad,xbsl,TA,multiEX
      save /gridd/

      integer pangfct
      common /gridp/ pangfct
      save /gridp/

c###########################################################################
c RADGRD file
c###########################################################################
c maxanggrid: The maximum number of different angular grids which can be
c             used in any given calculation.
c             ***NOTE*** This is a genuine constraint, but it must be used
c             since we must be able to keep a record of which angular grids
c             are used (before we have any allocated memory) since we have
c             to know how many grids are used in order to determine how much
c             memory to allocate.  This is set high enough it should never
c             be a problem.
c gridlist  : A list of all grids used (see the comment on maxanggrid).  It
c             is of dimension (maxanggrid,3) to keep track of the type and
c             subtype, and the number of times each grid is used.  The type
c             refers to how the grid is arrived at.
c                1 : Lebedev
c                2 : ?
c             The subtype refers to the degree of the grid of this type.
c numgrid   : The number of different angular grids used in the calculation.
c maxangpts : The maximum number of points in any of the angular grids used.
c maxanggrd : The grid with the maximum number of angular points.
c numradpts : The number of different radial points.
c ntotrad   : The total number of points in all angular grids at all radial
c             points (i.e. the entire integration grid)
c
c iradint   : determines if the Handy method (1) or Gauss-Legendre (2)
c             radial integration is used
c autosiz   : A flag which sets whether the polyhedra are (1) equally
c             sized, (2) sized according to Bragg-Slater radii or
c             (3) automatically sized according to the minimums in
c             density.
c slater    : A flag which determines whether (0) Slater's' rules are
c             used to determine the atomic size and scale the radial
c             integration or (1) no scaling is used.
c rigid     : A flag which determines whether rigid (0) or fuzzy
c             partitioning is used.
c nitr      : The number of iterations of the equations which create the
c             'fuzzy' boundary.

      integer maxanggrid
      parameter (maxanggrid=1000)

      integer gridlist(maxanggrid,3),numgrid,maxangpts,numradpts,
     &    iradint,autosiz,slater,rigid,nitr,maxanggrd,ntotrad

      common /radgrd/  gridlist,numgrid,maxangpts,numradpts,
     &    iradint,autosiz,slater,rigid,nitr,maxanggrd,ntotrad
      save /radgrd/

c Memory pointers
c
c pradgrid(numradpts) : The angular grid to use at each radial point.
c pgrdangpts(numgrid) : The number of angular points in each grid.
c zgridxyz(3,maxangpts,numgrid)
c                     : The x,y,z coordinate of each angular point in each grid
c zgridwt(maxangpts,numgrid)
c                     : The weight at each point.
c pintegaxis(natoms,3): Contains information about how much of each axis to
c                       integrate over.  If integaxis(iatom,i) is set to i,
c                       integrate only over the positive half of the i^th
c                       axis.  Otherwise, integrate over the entire axis.

      integer pgrdangpts,pradgrid,zgridxyz,zgridwt,pintegaxis

      common /radgrdp/ pgrdangpts,pradgrid,zgridxyz,zgridwt,
     &    pintegaxis
      save /radgrdp/

c###########################################################################
c Old stuff
c###########################################################################

c polist    : Contains an ordered list of unique atoms
c zatmvc    : The x, y, and z distance between each pair of atoms.
c zrij      : The distance between each pair of atoms.
c zatmpth   : The cartesian coordinates for the path integration between
c              all atom pairs
c zptdis    : The distance from atom i to a point along the path between
c              atoms i and j
c zprsqrd   : The distance squared from each atom to a point along all the
c              paths between all the atoms
c zpthpt    : The cartesian coordinates with respect to each atom for
c              the points along all the paths between all the atoms
c zbslrd    : The Bragg-Slater radii.
c zaij      : Surface shifting parameter dependent on the distance between
c               pairs of atoms.

      integer polist,ixx,iyy,izz,zatmvc,zrij,zatmpth,zptdis,
     &    zprsqrd,zpthpt,zbslrd,zaij
      common /gridold/ polist,ixx,iyy,izz,zatmvc,zrij,zatmpth,
     &    zptdis,zprsqrd,zpthpt,zbslrd,zaij
      save /gridold/




c This commonblock contains values for the functionals involved in the
c numerical integration either for DFT, plotting or fitting

c totele    : The total number of electrons
c etotekin  : The total kinetic energy
c etotenatr : The total nuclear-electron attraction energy
c etottf    : The total Thomas Fermi kinetic energy
c etotw     : The total Weizacker kinetic energy
c xldax     : The total LDA exchange energy
c becke     : The total Becke exchange energy
c lda       : The total LDA correlation energy
c xlyp      : The total LYP correlation energy
c icntr     : The integration center
c idns      : a flag =0 for SCF orbitals and =1 for natural orbitals

      integer icntr,idns
      double precision
     &    totele,etotkin,etotnatr,etottf,etotw,xldax,becke,
     &    lda,xlyp

      common /int_com/  icntr,idns
      save   /int_com/
      common /intr_com/ totele,etotkin,etotnatr,etottf,etotw,xldax,
     &                  becke,lda,xlyp
      save   /intr_com/

c array pointers

c zpcoeff(2): alpha/beta MO to primitive function transformation matrix
c zxocc     : alpha/beta orbital occupation

      integer
     &    zpcoeff(2),zxocc
      common /molecp/
     &    zpcoeff,zxocc
      save /molecp/




c The following values are used in determining the values of the various
c functionals:
c
c    roa    : The alpha density
c    rob    :
c    ro     : roa+rob
c    roinv  : 1/ro
c    rom    : roa-rob
c    rs     : (3 / 4 pi ro)**(1/3)
c    zeta   : (roa-rob)/ro
c
c    gradx  : d(roa)/dx + d(rob)/dx
c    grady  :
c    gradz  :
c    gradxm : d(roa)/dx - d(rob)/dx
c    gradym :
c    gradzm :
c    gro2   : gradx**2 + grady**2 + gradz**2
c    gro    : sqrt(gro2)
c    gro2a  : [d(roa)/dx]**2 + [d(roa)/dy]**2 + [d(roa)/dz]**2
c    gro2b  :
c
c    hesxx  : d2(roa)/dx2 + d2(rob)/dx2
c    hesyy  :
c    heszz  :
c    hesxy  :
c    hesxz  :
c    hesyz  :
c    hesxxm : d2(roa)/dx2 - d2(rob)/dx2
c    hesyym :
c    heszzm :
c    hesxym :
c    hesxzm :
c    hesyzm :
c    xlap   : hesxx + hesyy + heszz
c    xlapm  : hesxxm + hesyym + heszzm
c    trm1   :
c    trm1m  :
c    trm2   :
c    trm2m  :
c    grdaa  : gradroa.gradroa
c    grdbb  :
c    grdab  :
c    thresh : Parameter which tells what is the cutoff value
c           :  for is for the density, functionals and potetentials.
c           : The flag is cutoff and the default value is 10**(-12).
      double precision
     &    roa,rob,ro,rom,rs,zeta,roinv,
     &    gradx,grady,gradz,gradxm,gradym,gradzm,gro,gro2,gro2a,gro2b,
     &    hesxx,hesyy,heszz,hesxy,hesxz,hesyz,
     &    hesxxm,hesyym,heszzm,hesxym,hesxzm,hesyzm,
     &    xlap,xlapm,trm1,trm1m,trm2,trm2m,
     &    grdaa,grdbb,grdab,thresh

c There are several types of exchange and correlation energy functionals
c as well as hybrid schemes.  There are also corresponding potentials
c for use in Kohn-Sham.  Define the number of each type of energy functional
c and parameters to make calling them clearer.
c
c We want to have freedom to sellect potential and functional for
c better flexibility. Usually, the choice of functional determines 
c the poential  to be used in the  KS calculations.
c 
c
c fun_num_exch  : the number of exchange energy functionals (or potentials)
c fun_exch      : the total exchange energy contribution from the current
c                 atom
c tot_exch      : the total exchange energy contribution from all atoms
c coef_exch     : the coefficient of each exchange functional
c ene_exch      : the total energy with this functional
c nam_exch      : the name of the functional
c coef_pot_exch : the coefficient of each exchange component in the KS potential
c coef_pot_corr : the coefficient of each correlation component in the KS potential

      integer fun_num_exch,fun_num_corr,fun_num_hyb
      parameter (fun_num_exch = 5)
      parameter (fun_num_corr = 6)
      parameter (fun_num_hyb  = 1)

      double precision
     &    fun_exch(fun_num_exch),fun_corr(fun_num_corr),
     &    fun_hyb(fun_num_hyb),
     &    tot_exch(fun_num_exch),tot_corr(fun_num_corr),
     &    tot_hyb(fun_num_hyb),
     &    ene_exch(fun_num_exch),ene_corr(fun_num_corr),
     &    ene_hyb(fun_num_hyb),
     &    coef_exch(fun_num_exch),coef_corr(fun_num_corr),
     &    coef_pot_exch(fun_num_exch),coef_pot_corr(fun_num_corr)

      double precision vxc_ksalpha,vxc_ksbeta

      integer fun_exch_none,fun_corr_none,fun_special,
     &    fun_dft_none,
     &    fun_exch_lda,fun_exch_becke,fun_exch_pbe,fun_exch_pw91,
     &    fun_exch_hf,
     &    fun_corr_vwn,fun_corr_lyp,fun_corr_pbe,
     &    fun_corr_pw91, fun_corr_wl,fun_corr_wi,
     &    fun_hyb_b3lyp

      character*50
     &    nam_exch(fun_num_exch),nam_corr(fun_num_corr),
     &    nam_hyb(fun_num_hyb)
      character*80 nam_func,nam_kspot
      character*10
     &    abb_exch(fun_num_exch),abb_corr(fun_num_corr),
     &    abb_hyb(fun_num_hyb)

      parameter (fun_exch_none = 0)
      parameter (fun_exch_lda  = 1)
      parameter (fun_exch_becke= 2)
      parameter (fun_exch_pbe  = 3)
      parameter (fun_exch_pw91 = 4)
      parameter (fun_exch_hf   = 5)

      parameter (fun_corr_none = 0)
      parameter (fun_corr_vwn  = 1)
      parameter (fun_corr_lyp  = 2)
      parameter (fun_corr_pbe  = 3)
      parameter (fun_corr_pw91 = 4)
      parameter (fun_corr_wl   = 5)
      parameter (fun_corr_wi   = 6)

      parameter (fun_special   = 0)

      parameter (fun_dft_none  = -1)
      parameter (fun_hyb_b3lyp = 1)

      common /dftfunc/ fun_exch,fun_corr,fun_hyb,tot_exch,tot_corr,
     &    tot_hyb,ene_exch,ene_corr,ene_hyb,coef_exch,coef_corr,
     &    roa,rob,ro,rom,rs,zeta,roinv,
     &    gradx,grady,gradz,gradxm,gradym,gradzm,gro,gro2,gro2a,gro2b,
     &    hesxx,hesyy,heszz,hesxy,hesxz,hesyz,
     &    hesxxm,hesyym,heszzm,hesxym,hesxzm,hesyzm,
     &    xlap,xlapm,trm1,trm1m,trm2,trm2m,
     &    grdaa,grdbb,grdab,thresh,vxc_ksalpha,vxc_ksbeta,
     &    coef_pot_exch,coef_pot_corr

      common /dftfuncc/ nam_exch,nam_corr,nam_hyb,nam_func,
     &    abb_exch,abb_corr,abb_hyb,nam_kspot

      save /dftfunc/
      save /dftfuncc/

c Parameters for hybrid methods:
      double precision
     &    b3lypa,b3lypb,b3lypc
      parameter (b3lypa=0.20d0)
      parameter (b3lypb=0.72d0)
      parameter (b3lypc=0.81d0)



      integer ncount,ncnt,nnn,ifctr,zirrtmp,zradpt,zrwt,
     &    zwtintr,zsize,zrelsiz,zrtmp,zrint,setptr,iradpt,iangpt,
     &    grid,zvalmo,zcdnt,zrsqrd,zvalprim,zxnat,zgradprim,zgradmo,
     &    znull,pnull,zvalao,zvalgradao
      integer
     &   pmapatom,zgradS,zgradient,equ

      integer ppjj,imap2z,tgradient,NA
      integer basorder,zgradfinal,pchrgRatom

      integer zvalgrad2ao,zgrad2prim
      double precision
     &    ftotele,ftottf,ftotw,
     &    factor,totwt,ddot

       double precision
     &    exch_energy,corr_energy,exch_corr_energy,energy

      character*2 atmnam,atomsymb

      logical evalpt,init,print_post_ks
      logical dograd
      integer ihfdftgrad
      integer i, j,zdismap,zmapsymatom
      integer cosym,pk,xtm
      integer im2z,zgrad,zgrad2E,zcma,zcmB,zcc,zdd,zgradE

      double precision
     &    times(20)

      integer istat, fputc, mod
      character*1 twirl(8)

      integer zw,zdp,zdz,zp,zdmmuij,zdmmuji

      data twirl /'|','/','-','\\','|','/','-','\\'/


      call callstack_push('NUMINT')

c local pointers

      znull  =setptr(1,1,    1)
      pnull  =setptr(1,0, 1)

      zcdnt     = setptr(1,1, natoms*3)
      zrsqrd    = setptr(1,1, natoms)
      zvalprim  = setptr(1,1, totprim)
      zvalmo    = setptr(1,1, nbastot*2)
      zgradprim = setptr(1,1, totprim*3)
      zgradmo   = setptr(1,1, nbastot*3*2)

c zwtintr   : An intermediate in calculating the weight due to each atom.
c zradpt    : The array of radial points

      zxnat   = setptr(1,1, nbastot)
      zirrtmp = setptr(1,1, natoms)
      zradpt  = setptr(1,1, int_numradpts)
      zrwt    = setptr(1,1, int_numradpts*maxangpts)
      zwtintr = setptr(1,1, natoms)
      zsize   = setptr(1,1, natoms*natoms)
      zrelsiz = setptr(1,1, natoms*natoms)
      zrtmp   = setptr(1,1, natoms)
      zrint   = setptr(1,1, natoms*maxangpts)

c Pointers 
         zvalao     = setptr(1,1, naobasfn)
         zvalgradao = setptr(1,1, 3*naobasfn)

         zcmA=setptr(1,1,naobasfn*nbastot*iintfp)
         zcmB=setptr(1,1,naobasfn*nbastot*iintfp)
         zcc=setptr(1,1,naobasfn*naobasfn*(iuhf+1))
         zdd=setptr(1,1,naobasfn*nbastot*(iuhf+1))
         zgradE=setptr(1,1,natoms*3*natoms)
         zgradfinal=setptr(1,1,natoms*3*natoms)
         zgrad2E=setptr(1,1,natoms*3)
         zgrad  =setptr(1,1,natoms*3)
         im2z=setptr(1,0,natoms)
         pmapatom=setptr(1,0,natoms)
         zgradS=setptr(1,1,natoms*3)

         tgradient=setptr(1,1,natoms*3)

         zdismap=setptr(1,1,natoms*natoms)
         zmapsymatom=setptr(1,0,natoms*natoms*natoms)
         pchrgRatom=setptr(1,0,natoms)
         ppjj=setptr(1,0,natoms*natoms)
         imap2z=setptr(1,0,natoms)
         zw       =setptr(1,1,3*natoms)
         zdp      =setptr(1,1,3*natoms*natoms)
         zdz      =setptr(1,1,3*natoms)
         zp       =setptr(1,1,natoms)
         zdmmuij   =setptr(1,1,3)
         zdmmuji  =setptr(1,1,3)

         zgrad2prim=setptr(1,1,totprim*9)
         zvalgrad2ao=setptr(1,1,naobasfn*9) 

         basorder=setptr(1,0,naobasfn)
         equ=setptr(1,0,natoms*natoms) 
         zgradient=setptr(1,1   ,natoms*3)
        
          cosym=setptr(1,0,3)
          pk=setptr(1,0,natoms)
          xtm=setptr(1,0,3) 


       call getrec(1,'JOBARC','CNTERBF0',naobasfn,
     &              kscore(basorder))
       call getrec(1,'JOBARC','CCOEMOA',naobasfn*nbastot*iintfp,
     &              dcore(zcmA))
       call getrec(-1,'JOBARC','CCOEMOB',naobasfn*nbastot*iintfp,
     &              dcore(zcmB))
       call mkdenpra(dcore(zcmA),dcore(zcc),dcore(zxocc),dcore(zdd),
     &               dcore(zcmB),kscore(basorder))


c      call mkpra(dcore(zprimcoef))
       call getrec(20,'JOBARC','HFDFTGRA',1,ihfdftgrad)
       if(ihfdftgrad.eq.1) then   
        do i=1,fun_num_exch
            coef_pot_exch(i)=coef_exch(I)
         end do
         do i=1,fun_num_corr
            coef_pot_corr(i)=coef_corr(i) 
         end do
       end if
      coef_pot_exch(fun_exch_hf)=0.d0
      
      call intpath(kscore(pnshellatom),kscore(pnshellprim),
     &             dcore(zalpha),dcore(zprsqrd),dcore(zvalprim),
     &             dcore(zpcoeff(1)),dcore(zvalmo),
     &             dcore(zpthpt),dcore(zsize),dcore(zptdis),
     &             dcore(zrelsiz),dcore(zbslrd),kscore(pangfct))



       call dzero(dcore(zgradE),natoms*3*natoms)
       call dzero(dcore(zgradS),natoms*3)
       call izero(kscore(zmapsymatom),natoms*natoms*natoms)
       call izero(kscore(ppjj),natoms*natoms)

c Set up ordered list of unique atoms
      ncount=compnorb
      if (c1symmet.eq.1) ncount=natoms
c       ncount=natoms
      do i=1,20
         times(i)=0.d0
      end do
      ftotele = 0.d0
      ftottf  = 0.d0
      ftotw   = 0.d0
        if( natoms .ne. ncount) then
      call equi(kscore(pcompmemb),ncount,kscore(pcomppopv),kscore(equ))
      call dismap(kscore(pcomppopv),kscore(equ),ncount,dcore(zdismap),
     &            dcore(zcoord))
      call symatomap(dcore(zdismap),kscore(pcomppopv),kscore(equ),
     &               ncount,kscore(zmapsymatom),kscore(ppjj),
     &             dcore(zcoord),kscore(cosym),kscore(pk),
     &             kscore(xtm))
       end if
      do ncnt=1,ncount
         nnn=ncnt

c          if( ncnt .eq. 8 ) then
c        Determine center and symmetry factor
c           (number of symmetry equivalent atoms)
         icntr=kscore(polist+ncnt-1)
c         ifctr=kscore(pcomppopv+ncnt-1)
c         factor=dble(ifctr)
         if (c1symmet.eq.1) then
            icntr=ncnt
            ifctr=1
            factor=1.d0
         end if
c        icntr=ncnt
        factor=1.d0

c        Set up grid for integrating only symmetry unique octants
         call octag(kscore(patomchrg),dcore(zatmvc),
     &            dcore(zrij),     kscore(pradgrid),
     &            dcore(zrelsiz),  dcore(zrtmp),
     &            dcore(zrint),    dcore(zradpt),
     &            dcore(zrwt),     dcore(zbslrd),
     &            dcore(zgridxyz), kscore(pgrdangpts))

         init=.true.
         do iradpt=1,int_numradpts
            grid=kscore(pradgrid+iradpt-1)
            do iangpt=1,kscore(pgrdangpts+grid-1)
               evalpt=.true.
          call symoctag(kscore(patomchrg), dcore(zatmvc),
     &                     dcore(zrij),      dcore(zaij),
     &                     dcore(zcdnt),     dcore(zrsqrd),
     &                     dcore(zirrtmp),   dcore(zwtintr),
     &                     totwt,            kscore(pradgrid),
     &                     dcore(zradpt),    dcore(zrwt),
     &                     kscore(pintegaxis),dcore(zgridxyz),
     &                     dcore(zgridwt),   kscore(pgrdangpts),
     &                     iradpt,           iangpt,
     &                     grid,             evalpt,dcore(zw),
     &                    dcore(zdp),dcore(zdz),dcore(zp),
     &                    dcore(zdmmuji))

               if (evalpt) then
                  call integAG(kscore(pnshellatom),kscore(pnshellprim),
     &                       dcore(zcdnt),      dcore(zalpha),
     &                       dcore(zrsqrd),     dcore(zvalprim),
     &                       dcore(zpcoeff(1)), dcore(zvalmo),
     &                       totwt,             dcore(zgradprim),
     &                       dcore(zgradmo),    dcore(zxnat),
     &                       dcore(zxocc),      
     &                       factor,            ifctr,
     &                       dcore(zvalao),     dcore(zvalgradao),
     &                       dcore(zprimcoef),  dcore(zcoord),
     &                       kscore(pcompmemb),  kscore(pangfct),
     &                            init,
     &                       iradpt,            iangpt,
     &                       times,             dograd,
     &                       dcore(zcc),dcore(zgradE),kscore(pnaoatom),
     &                       dcore(zw),dcore(zgrad2E),
     &         dcore(zgrad2prim),dcore(zvalgrad2ao),kscore(basorder))
               end if
c           end do iangpt=1,kscore(pgrdangpts+grid-1)
            end do
c        end do iradpt=1,int_numradpts
         end do

c      end if
        if( natoms .ne. ncount) then
      call gradsymgen(ncnt,icntr,dcore(zgradE),dcore(zcoord),
     &     kscore(pcomppopv),ncount,kscore(equ),
     &  kscore(zmapsymatom))
      end if
c     end do ncnt=1,ncount
      end do

      call gradsymsum(ncount,dcore(zgradE),dcore(zgradfinal),
     &      dcore(zcoord),kscore(pcomppopv),kscore(equ),
     &      kscore(pmapatom),dcore(zgradS))
      CALL GETREC(20,'JOBARC','NATOMS',1,NA)

           call analgrad(dcore(zgradS),dcore(zgradient),
     &          dcore(zgrad),dcore(zgrad2E),kscore(im2z),
     &           kscore(patomchrg),kscore(imap2z),
     &           kscore(pchrgRatom),dcore(tgradient),NA)


      call relptr(1,1,   znull)
      call relptr(1,0,pnull)

      call callstack_pop
      return
      end

