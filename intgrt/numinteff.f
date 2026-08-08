











      subroutine numinteff(valao,valgradao,totwt,dograd,
     &             intnumradpts,max_angpts,ncount,zksvxc,kshf)

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



      integer intnumradpts,max_angpts,rshift

      integer ncount,ncnt,nnn,ifctr,zirrtmp,zradpt,zrwt,
     &    zwtintr,zrtmp,zrint,setptr,iradpt,iangpt,
     &    grid,zvalmo,zvalprim,zxnat,zgradprim,zgradmo,
     &    znull,pnull,zvalao,zvalgradao,zscrnn,zfourind,
     &    zorder
      integer
     &    zksdint, zdens, zgradient, zgradscr, zgradord,
     &    ihfdftgrad

      double precision
     &    ftotele,ftottf,ftotw,
     &    factor,ddot

       double precision
     &    exch_energy,corr_energy,exch_corr_energy,energy

      double precision
     &    zksvxc(naobasfn,naobasfn,iuhf+1)

      character*2 atmnam,atomsymb

      logical evalpt,init,print_post_ks
      logical dograd, kshf, hfdftgrad

      integer i, j, naonao
      integer zdensity,zgraddensity

      double precision
     &    times(20)

      double precision
     & valao(naobasfn,max_angpts,intnumradpts,ncount),
     & valgradao(naobasfn,max_angpts,intnumradpts,ncount),
     & totwt(ncount,intnumradpts,max_angpts)

      integer  zcmA,zcmB
     
      integer groupoint

      integer istat, fputc, mod
      character*1 twirl(8)
      data twirl /'|','/','-','\\','|','/','-','\\'/


      call callstack_push('NUMINTNEW') 
c local pointers
      znull  =setptr(1,1,    1)
      pnull  =setptr(1,0, 1)

c zvalprim  : the value of each primitive function at the current point
c zvalmo    : the value of each alpha/beta MO at the current point
c zvalao    : the value of each AO at the current point
c zvalgradao: the value of the x/y/z gradient component of each AO at
c           : the current point
c zgradprim : the value of the x/y/z gradient component of each primitive
c             function at the current point
c zgradmo   : the value of the alpha/beta MO x/y/z gradient component at the
c             current point
c       ncount=compnorb
c      if (c1symmet.eq.1) ncount=natoms
c prakash not going for unique atoms


      zvalmo    = setptr(1,1, nbastot*maxangpts*2)
      zgradmo   = setptr(1,1, nbastot*maxangpts*3*2)

      zcmA      =setptr(1, 1,  nbastot*naobasfn*iintfp)
      zcmB      =setptr(1, 1,  nbastot*naobasfn*iintfp) 
      zdensity  =setptr(1,   1, 1*maxangpts*2)
      zgraddensity=setptr(1, 1, 1*maxangpts*3*2)

c zwtintr   : An intermediate in calculating the weight due to each atom.
c zradpt    : The array of radial points

      zxnat   = setptr(1,1, nbastot*maxangpts)
      zirrtmp = setptr(1,1, natoms)
      zradpt  = setptr(1,1, int_numradpts)
      zrwt    = setptr(1,1, int_numradpts*maxangpts)
      zwtintr = setptr(1,1, natoms)
      zrtmp   = setptr(1,1, natoms)
      zrint   = setptr(1,1, natoms*maxangpts)
c Pointers 
      zfourind   = znull
      zscrnn     = znull
      zorder     = pnull
c      if (int_ks) then
c         zksvxc     = setptr(1,1,naobasfn*naobasfn*(iuhf+1))
c      else
         if (int_tdks) then
            zvalao   = setptr(1, 1, naobasfn)
            naonao   = rshift((nbastot*(nbastot+1)),1)
            zscrnn   = setptr(1, 1, naonao)
            zfourind = setptr(1, 1, naonao*(naonao+1)*(iuhf+1)/2)
            zorder   = setptr(1, 0, naobasfn)
         end if 
c      end if

c Kohn-Sham gradient pointers
      if (dograd) then
         zksdint  =setptr(1,1   ,naobasfn*naobasfn*3*(iuhf+1))
         zdens    =setptr(1,1   ,naobasfn*naobasfn*(iuhf+1))
         zgradscr =setptr(1,1   ,naobasfn*naobasfn*3)
         zgradient=setptr(1,1   ,natoms*3)
         zgradord =setptr(1,0,naobasfn*2)
      else
         zksdint  =znull
         zdens    =znull
         zgradscr =znull
         zgradient=znull
         zgradord =pnull
      end if

      if (memknown.eq.0) goto 999

       call getrec(1,'JOBARC','CCOEMOA',nbastot*naobasfn*iintfp,
     &               dcore(zcmA))

       call getrec(-1,'JOBARC','CCOEMOB',nbastot*naobasfn*iintfp,
     &               dcore(zcmB))
       
      if (int_ks) then
         call dzero(zksvxc,naobasfn*naobasfn*(iuhf+1))
         if (dograd) then
c         o zero the derivative-integral matrix
            call dzero(dcore(zksdint),naobasfn*naobasfn*3*(iuhf+1))
c         o set the HF exchange coefficient to zero since VDINT already
c           calculated its contribution
            coef_pot_exch(fun_exch_hf)=0.d0
         end if
      else
         if (int_tdks) then
            call dzero(dcore(zfourind),naonao*(naonao+1)*(iuhf+1)/2)
            call getrec(1,'JOBARC','KSSRTORD',naobasfn,kscore(zorder))
            call getrec(1,'JOBARC','CMP2ZMAT',naobasfn*nbastot*iintfp,
     &                  zksvxc)
         end if
      end if

 
c Loop over integration centers (nuclei) and set the radial grid using
c Slater's rules.

c Set up ordered list of unique atoms
c      ncount=compnorb
c      if (c1symmet.eq.1) ncount=natoms

      do i=1,20
         times(i)=0.d0
      end do
C
      call getrec(20,'JOBARC','HFDFTGRA',1,ihfdftgrad)
      if (ihfdftgrad .eq. 1 ) then
          hfdftgrad=.true.
      else
          hfdftgrad=.false.
      end if
C
      if(hfdftgrad.or.kshf) then
        do i=1,fun_num_exch
            coef_pot_exch(i)=coef_exch(I)
        end do
C
         do i=1,fun_num_corr
            coef_pot_corr(i)=coef_corr(i)
         end do
C
         coef_pot_exch(fun_exch_hf)=0.d0
      end if
C
      ftotele = 0.d0
      ftottf  = 0.d0
      ftotw   = 0.d0
      do ncnt=1,ncount
         nnn=ncnt

c        Determine center and symmetry factor
c           (number of symmetry equivalent atoms)


c praskash start
c  not going for unique atoms

         icntr=kscore(polist+ncnt-1)
         ifctr=kscore(pcomppopv+ncnt-1)
c         write(*,*) 'ifctr=',ifctr
         factor=dble(ifctr)
c         write(*,*) 'factor=',factor
         if (ncount.eq.natoms) then
            icntr=ncnt
            ifctr=1
            factor=1.d0
         end if


c prakash end
         init=.true.
         do iradpt=1,int_numradpts

            grid=kscore(pradgrid+iradpt-1)
            
               groupoint=kscore(pgrdangpts+grid-1)
         call findmo(dcore(zcmA),dcore(zcmB),
     &    dcore(zvalmo),
     &    dcore(zgradmo),dcore(zxnat),
     &    dcore(zxocc),ncount,iradpt,groupoint,
     &    dcore(zdensity),dcore(zgraddensity),
     &    valao,valgradao,intnumradpts,max_angpts,ncnt)

            do iangpt=1,kscore(pgrdangpts+grid-1)
               evalpt=.true.

               if (evalpt) then
       call integ_eff1(
     &              totwt,     
     &               zksvxc,
     &                factor,            ifctr,
     &             valao,     valgradao,
     &        dcore(zfourind),   dcore(zscrnn),
     &                init,ncount,
     &           iradpt,            iangpt,
     &     times,             dcore(zksdint),dograd,
     &     dcore(zdensity),dcore(zgraddensity),groupoint,
     &        intnumradpts,max_angpts,ncnt,kshf)
               end if
c           end do iangpt=1,kscore(pgrdangpts+grid-1)
            end do
c        end do iradpt=1,int_numradpts
         end do

c Print the value of the integrated density

         if (.not.int_ks) then
            if (int_printatom.gt.int_printlev) then
               atmnam=atomsymb(kscore(patomchrg+icntr-1))
               write(*,2000) icntr,atmnam,ifctr
               write(*,2001) totele
            end if
         end if

 2000    format(/'Atom',i3,' is ',a2,', it involves ',i3,
     &          ' equivalent atoms')
 2001    format('   The electron density                   =',
     &          f20.12,' electrons')

c        Add atomic result to total

         ftotele=ftotele+factor*totele
         do i=1,fun_num_exch-1
            tot_exch(i)=tot_exch(i)+factor*fun_exch(i)
         end do
         do i=1,fun_num_corr
            tot_corr(i)=tot_corr(i)+factor*fun_corr(i)
         end do

         ftottf=ftottf+factor*etottf
         ftotw=ftotw+factor*etotw

c     end do ncnt=1,ncount
      end do

c The energy corresponding to fun_exch_hf is ehfx

      tot_exch(fun_exch_hf)=ehfx

c Calculate hybrid values

      tot_hyb(fun_hyb_b3lyp)=
     &      (b3lypa*ehfx)+
     &      (1.d0-b3lypa)*tot_exch(fun_exch_lda)+
     &      b3lypb*(tot_exch(fun_exch_becke)-tot_exch(fun_exch_lda))+
     &      tot_corr(fun_corr_vwn)+
     &      b3lypc*(tot_corr(fun_corr_lyp)-tot_corr(fun_corr_vwn))

      if (int_ks) then
         if (dograd) then
            call dintproc(
     &              dcore(zksdint),
     &              dcore(zdens),
     &              dcore(zgradscr),
     &              dcore(zgradient),
     &              naobasfn,nbastot,
     &              dcore(zgradscr+naobasfn*naobasfn),
     &              kscore(zgradord),
     &              kscore(zgradord+naobasfn))
         else
            call putrec(1,'JOBARC','KSTOTELE',iintfp,ftotele)
c            call putrec(1,'JOBARC','KSPOTENT',
c     &          naobasfn*naobasfn*(iuhf+1)*iintfp,dcore(zksvxc))
         end if
      end if

c Print the value of the integrated density and functionals
      if (.not.int_ks) then

c         if (int_tdks) then
c We shall dump the integrals for needed for vtran
c        HERE!!!!!!!
c         end if

         write(*,2100) ftotele
c        write(*,2107) ftottf
c        write(*,2108) ftotw

         do i=1,fun_num_exch
            write(*,2121) nam_exch(i),tot_exch(i)
         end do
         do i=1,fun_num_corr
            write(*,2120) nam_corr(i),tot_corr(i)
         end do
         do i=1,fun_num_hyb
            write(*,2122) nam_hyb(i),tot_hyb(i)
         end do

c         print *, 'DEBUG: POTENTIAL/FUNCTIONAL COEFFICIENTS'
c         do i=1,fun_num_exch
c            write(*,2123) nam_exch(i),coef_pot_exch(i),coef_exch(i)
c         end do
c         do i=1,fun_num_corr
c            write(*,2123) nam_corr(i),coef_pot_corr(i),coef_corr(i)
c         end do

         write(*,2199)

         exch_energy=ddot(fun_num_exch,tot_exch,1,coef_pot_exch,1)
         corr_energy=ddot(fun_num_corr,tot_corr,1,coef_pot_corr,1)
         exch_corr_energy=exch_energy+corr_energy
         energy=ehar+exch_corr_energy

c Check if the density is from a KS calculation.
c If so, print the energies from the KS calculation.

         call getrec(1,'JOBARC','KSPRINT',iintfp,print_post_ks)

         if (int_dft_fun.gt.fun_dft_none) then
            if (int_kspot.eq.fun_hyb_b3lyp) write(*,6664)
            exch_energy=0.d0
            corr_energy=0.d0
            exch_corr_energy=0.d0
            energy=0.d0

            exch_energy=ddot(fun_num_exch,tot_exch,1,coef_exch,1)
            corr_energy=ddot(fun_num_corr,tot_corr,1,coef_corr,1)

            exch_corr_energy=exch_energy+corr_energy
            energy=ehar+exch_corr_energy 

            write(*,6676) ehar
            write(*,6677) exch_energy
            write(*,6678) corr_energy
            write(*,6675) exch_corr_energy
            write(*,6679) energy  
         end if

         call putrec(1,'JOBARC','TOTENERG',iintfp,energy)

c      o write out the HF exchange coefficient for use in VDINT
         call putrec(1,'JOBARC','PCHF',
     &       iintfp,coef_pot_exch(fun_exch_hf))

 2100 format(/'   Total density integrates to             : ',
     &    f20.12,' electrons')
 2107 format('The total TF kinetic energy                  =',
     &    f20.12,' a.u.')
 2108 format('The total Weizsacker kinetic energy          =',
     &    f20.12,' a.u.')

 2120 format('   Corr tot contrib : ',A20,' : ',f20.12)
 2121 format('   Exch tot contrib : ',A20,' : ',f20.12)
 2122 format('   Hybr tot contrib : ',A20,' : ',f20.12)
 2123 format(A20,' : ',f20.12,', ',f20.12)

c Print some total energies
 2320 format('   Corr total energy: ',A20,' : ',f20.12)
 2321 format('   Exch total energy: ',A20,' : ',f20.12)
 2322 format('   Hybr total energy: ',A20,' : ',f20.12)

 6664 format('Gaussian94-like B3LYP with its SCF Density (Full VWN)'
     &,/,
     &t1,'The original B3LYP corresponds to func=b3lyp, kspot=lda,vwn'/) 
 6665 format(' Exc  from   KS   Calculation:     ',f20.12,' a.u.')
 6666 format(' Ehar from   KS   Calculation:     ',f20.12,' a.u.')
 6667 format(' Ex   from   KS   Calculation:     ',f20.12,' a.u.')
 6668 format(' Ec   from   KS   Calculation:     ',f20.12,' a.u.')
 6669 format(' Final Result from KS Calculation: ',f20.12,' a.u.'/)

 6675 format(' Exc  with the SCF Dens and Func:  ',f20.12,' a.u.')
 6676 format(' Ehar with the SCF Dens and Func:  ',f20.12,' a.u.')
 6677 format(' Ex   with the SCF Dens and Func:  ',f20.12,' a.u.')
 6678 format(' Ec   with the SCF Dens and Func:  ',f20.12,' a.u.')
 6679 format(' Final Result with the SCF Dens and Func:',
     &f20.12,' a.u.'/)

      end if

 2199 format('     ')

  999 continue
      call relptr(1,1,   znull)
      call relptr(1,0,pnull)

      call callstack_pop
      return
      end

