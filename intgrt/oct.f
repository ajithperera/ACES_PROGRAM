










      subroutine oct(atomchrg,atmvc,rij,
     &    radgrid,relsiz,rtmp,rint,radpt,rwt,
     &    bslrd,gridxyz,grdangpts)

c This routine sets up the numerical integration grid

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


      integer
     &    atomchrg(natoms),radgrid(int_numradpts),
     &    grdangpts(numgrid)
      double precision
     &    atmvc(natoms,natoms,3),
     &    rij(natoms,natoms),
     &    relsiz(natoms,natoms),rtmp(natoms),
     &    rint(natoms,maxangpts),radpt(int_numradpts),
     &    rwt(int_numradpts,maxangpts),bslrd(natoms),
     &    gridxyz(3,maxangpts,numgrid)
      double precision pi
      
      double precision
     &    x23i(23),a23i(23),R,x50i(50),a50i(50),
     &    x45i(45),a45i(45),x40i(40),a40i(40),
     &    x35i(35),a35i(35),x30i(30),a30i(30)

      integer
     &    iradpt,iangpt,i,midpt,j,jcntr
      double precision
     &    xnradp1,rad,xiradpt,eps,
     &    xn,xinlngth,xinl,xmid,xi,zz,p1,p2,xj,p3,pp,zz1,
     &    rwttmp,xoulngth

      double precision
     &    FbyT,SbyT,alpha,rad1,w1,w2,w3,
     &     w4,w5,w6,w7,w8,w9,x,diradpt,a,one,two,three
         parameter(one=1.0D0)
         parameter(two=2.0d0)
         parameter(three=3.0d0)

       data(x30i(i),i=1,30)/0.000938253265447d0,0.006107726543114d0,
     & 0.016031080676759d0,0.030705186162534d0,0.050036197942358d0,
     & 0.073865597603489d0,0.101979339965684d0,0.134113324110492d0,
     & 0.169957816772565d0,0.209161640244516d0,0.251336426144445d0,
     & 0.296061051010481d0,0.342886289952091d0,0.391339685613067d0,
     & 0.440930608618850d0,0.491155473348428d0,0.541503065266541d0,
     & 0.591459931268120d0,0.640515781563337d0,0.688168850061483d0,
     & 0.733931159722301d0,0.777333639800269d0,0.817931043294885d0,
     & 0.855306615414549d0,0.889076468144764d0,0.918893624471120d0,
     & 0.944451717719704d0,0.965488410519963d0,0.981789011002416d0,
     & 0.993194470701601d0/

       data(a30i(i),i=1,30)/0.137676443941630d0,0.195885975601622d0,
     & 0.210268725136789d0,0.206561492289121d0,0.193856500144686d0,
     & 0.176584044933802d0,0.157256165869557d0,0.137415766886845d0,
     & 0.118046532595713d0,0.099777943557785d0,0.082999750605787d0,
     & 0.067931912187925d0,0.054670739318578d0,0.043221425989051d0,
     & 0.033522285387789d0,0.025463611447442d0,0.018902826458941d0,
     & 0.013676886808064d0,0.009612528011925d0, 0.006534705005718d0,
     &  0.004273451964695d0,0.002669310941103d0,0.001577437115953d0,
     & 0.000870468589335d0,0.000440241650354d0,0.000198433961056d0,
     & 0.000076224316174d0,0.000023066690603d0,0.000004686380402d0,
     & 0.000000416211557d0/

       data (x23i(i),i=1,23)/
     &   0.001505892474584,0.009949112846861,
     &   0.026212787562514,0.050215014094684,
     &   0.081660512821456,0.120092694277137,
     &   0.164915797166005,0.215411397226014,
     &   0.270754073298502,0.330027592716307,
     &   0.392241993524124,0.456351572232027,
     &   0.521273620526389,0.585907670629421,
     &   0.649154963911799,0.709937833943106,
     &   0.767218686284763, 0.820018260579743,
     &   0.867432877982759, 0.908650423041717,
     &   0.942964955771737,0.969790557585917,
     &   0.988681214124179/

        data(a23i(i),i=1,23)/
     &  0.193979975198118,0.262363963659648,
     & 0.267428126763970, 0.248662024728409,
     & 0.219828495609352,0.187495823926848,
     & 0.155234272915304,0.125066174707640,
     & 0.098100418872551,0.074860364849399,
     & 0.055477247414685,0.039817159733692,
     & 0.027571576765108,0.018325604666402,
     & 0.011611019656801,0.006947749930384,
     & 0.003875771461057,0.001978563506850,
     &  0.000898892530630,0.000347557889758,
     &  0.000105728680567,0.000021565953939,
     & 0.000001920578888/

       data(a35i(i),i=1,35)/0.112092876985421d0, 0.163605423965229d0,
     & 0.180057148429348d0,0.181471124457663d0,0.174971287793732d0,
     & 0.164053261027971d0,0.150722110311121d0,0.136234748428706d0,
     & 0.121420715568996d0,0.106841864276382d0,0.092880478278897d0,
     & 0.079792025456380d0,0.067738974224230d0,0.056813802853714d0,
     & 0.047055502143372d0,0.038461968941276d0,0.030999683895043d0,
     & 0.024611510338003d0,0.019223130211810d0,0.014748441657154d0,
     &  0.011094126050064d0,0.008163519610316d0,0.005859879132167d0,
     & 0.004089102843288d0,0.002761949813728d0,0.001795790936293d0,
     & 0.001115918869016d0,0.000656441886258d0,0.000360786261700d0,
     & 0.000181832868912d0,0.000081715613329d0,0.000031311727232d0,
     & 0.000009456588660d0,0.000001918357355d0,0.000000170197233d0/

       data(x35i(i),i=1,35)/0.000710288160721d0,0.004587102256450d0,
     & 0.012012783293299d0,0.022997219083556d0,0.037494445348316d0,
     & 0.055420399987502d0,0.076658733330740d0,0.101063934846327d0,
     & 0.128463607203279d0,0.158660460020927d0,0.191434244320567d0,
     & 0.226543722395174d0,0.263728714255856d0,0.302712236149867d0,
     &  0.343202733087210d0,0.384896399547925d0,0.427479577789020d0,
     &  0.470631220106360d0,0.514025399364001d0,0.557333850725481d0,
     & 0.600228526610458d0,0.642384146341758d0,0.683480721674270d0,
     & 0.723206039367624d0,0.761258082155899d0,0.797347369872117d0,
     & 0.831199203116472d0,0.862555792773512d0,0.891178260062636d0,
     & 0.916848494189018d0,0.939370859849387d0,0.958573761748160d0,
     & 0.974311130330776d0,0.986464197666822d0,0.994946681075095d0/

       data(x40i(i),i=1,40)/0.000556967498706d0,0.003573534297375d0,
     & 0.009340342752575d0,0.017870745434866d0,0.029141100831243d0,
     & 0.043103594816730d0,0.059690257935624d0,0.078814978080940d0,
     & 0.100374843362374d0,0.124251233458652d0,0.150310823319684d0,
     & 0.178406572113624d0,0.208378732027199d0,0.240055893384047d0,
     & 0.273256073107231d0,0.307787848247445d0,0.343451533113967d0,
     & 0.380040396522862d0,0.417341914330888d0,0.455139051495079d0,
     & 0.493211567242711d0,0.531337336470804d0,0.569293680168549d0,
     & 0.606858697439291d0,0.643812591571502d0,0.679938982557982d0,
     & 0.715026198481030d0,0.748868538263425d0,0.781267498428097d0,
     & 0.812032956714113d0,0.840984305669462d0,0.867951529701627d0,
     & 0.892776219567142d0,0.915312519064272d0,0.935428000181421d0,
     & 0.953004466542159d0,0.967938696011113d0, 0.980143176810235d0,
     & 0.989547125989027d0,0.996100198799301d0/

      data(a40i(i),i=1,40)/0.093465963789691d0,0.139195760758708d0,
     & 0.156165705182621d0,0.160459227209418d0,0.157833390206324d0,
     & 0.151126918069830d0,0.141978734962672d0, 0.131427129313022d0,
     & 0.120169213151007d0,0.108690050821019d0,0.097333650698496d0,
     & 0.086345110970585d0,0.075897263990631d0,0.066108452722889d0,
     & 0.057054968845518d0,0.048780135326969d0,0.041301197466255d0,
     & 0.034614730691338d0,0.028701008962378d0,0.023527618658139d0,
     & 0.019052504390898d0,0.015226570732274d0,0.011995923407948d0,
     & 0.009303806970138d0,0.007092278358123d0,0.005303644058252d0,
     & 0.003881680839129d0,0.002772655007931d0,0.001926151976128d0,
     & 0.001295726072580d0,0.000839379609339d0,0.000519879915638d0,
     & 0.000304923211442d0,0.000167154649196d0,0.000084054503519d0,
     & 0.000037701251993d0,0.000014423103218d0,0.000004350341946d0,
     & 0.000000881637013d0,0.000000078165783/
     
      data(a45i(i),i=1,45)/0.079412084797958d0,0.120227864969661d0,
     & 0.136973736570282d0,0.142889769940632d0,0.142740692715896d0,
     & 0.138885029284403d0,0.132690322698210d0,0.125028305231930d0,
     & 0.116489871765386d0,0.107492204798495d0,0.098337637463301d0,
     & 0.089248493333103d0,0.080388975566380d0,0.071879625850580d0,
     & 0.063807300655344d0,0.056232327633949d0,0.049193823523536d0,
     & 0.042713774751944d0,0.036800260739957d0,0.031450066423046d0,
     & 0.026650847463803d0,0.022382958549859d0,0.018621020470679d0,
     & 0.015335278537374d0,0.012492789246762d0,0.010058461351638d0,
     & 0.007995970071228d0,0.006268558014830d0,0.004839732807681d0,
     & 0.003673868936701d0,0.002736719656914d0,0.001995843695529d0,
     & 0.001420950803110d0,0.000984169817301d0,0.000660242741836d0,
     & 0.000426648341131d0,0.000263658862601d0, 0.000154333689953d0,
     & 0.000084453973862d0,0.000042402560342d0,0.000018993824772d0,
     & 0.000007258307489d0,0.000002187324000d0,0.000000442980480d0,
     & 0.000000039256135/

        data(x45i(i),i=1,45)/0.000448797230789d0,0.002863743176724d0,
     & 0.007472580719957d0,0.014288795399632d0,0.023299971981916d0,
     & 0.034477447370902d0,0.047779250073919d0,0.063151503776824d0,
     & 0.080529304021336d0,0.099837383985378d0, 0.120990693832148d0,
     & 0.143894949860442d0,0.168447181097710d0,0.194536287513316d0,
     & 0.222043617088493d0,0.250843565126355d0,0.280804196907227d0,
     & 0.311787893371681d0,0.343652018581199d0,0.376249607069379d0,
     & 0.409430068746358d0,0.443039908692630d0,0.476923458937795d0,
     & 0.510923619141852d0,0.544882602966637d0,0.578642686833388d0,
     & 0.612046957702879d0,0.644940056482826d0,0.677168913660367d0,
     & 0.708583473773287d0,0.739037405370911d0,0.768388793173291d0,
     & 0.796500809215335d0,0.823242359861509d0,0.848488705699473d0,
     & 0.872122051474283d0,0.894032103425607d0,0.914116591681553d0,
     & 0.932281755859211d0,0.948442793061117d0,0.962524270148234d0,
     & 0.974460510755284d0,0.984196001885642d0,0.991686051091295d0,
     & 0.996899609726306d0/

       data (x50i(i),i=1,50)/0.000369570433119,0.002347134219197,
     & 0.006115571770995,0.011687370845523,0.019056082830034,
     & 0.028203824821908,0.039103519166679,0.051719928318451,
     & 0.066010270305174,0.081924660297733,0.099406475155824,
     & 0.118392685022626,0.138814173979560,0.160596061423103,
     & 0.183658030519096,0.207914667166511,0.233275811191679,
     & 0.259646920438599,0.286929447739133,0.315021230287611,
     & 0.343816890623593,0.373208248195249,0.403084740304811,
     & 0.433333851108568,0.463841547245645,0.494492718594754,
     & 0.525171622601351,0.555762330575936,0.586149174335211,
     & 0.616217191539941,0.645852568075473,0.674943075822135,
     & 0.703378504172559,0.731051083670879,0.757855900174449,
     &  0.783691297972073,0.808459270333679,0.832065836015480,
     & 0.854421400302878,0.875441099243377,0.895045125809644,
     &  0.913159036852631,0.929714039891706,0.944647259136955,
     & 0.957901980923149,0.969427880858720,0.979181241908920,
     & 0.987125200529378,0.993230209396203,0.997476276077283/

       data(a50i(i),i=1,50)/0.068504136804909,0.105151577875052,
     & 0.121321128315539,0.128125986014206,0.129586021833051,
     & 0.127695678926300,0.123617666227836,0.118092991974382,
     & 0.111630546226554,0.104588588861948,0.097234493302991,
     & 0.089763620635523,0.082328865185079,0.075041385762480,
     & 0.067990500989898,0.061237978454549,0.054833878003905,
     & 0.048808520509080,0.043186012720081,0.037975680674068,
     & 0.033183830323399,0.028805524868256,0.024834759426066,
     & 0.021257033108662,0.018058012776293,0.015217112796852,
     & 0.012714715943295,0.010526827220442,0.008630943425451,
     &  0.007001753076647,0.005615774001052,0.004448011377657,
     & 0.003475504785264,0.002674820983222,0.002024664852130,
     & 0.001504066542102,0.001094214437320,0.000777193141826,
     & 0.000537192101604,0.000359657965941,0.000232021425044,
     & 0.000143141823139,0.000083678850721,0.000045726421112,
     & 0.000022936769507,0.000010262444070,0.000003919661657,
     & 0.000001180004688,0.000000238999866,0.000000021149283/


      call callstack_push('OCT')

      pi = acos(-1.d0)
      eps=3.d0*1.d-14
      xnradp1=dble(int_numradpts+1)
   
      alpha=0.6D0
      FbyT=5.0D0/2.0D0
      SbyT=7.0D0/2.0D0


c Use Slater radii for radial integration scaling
       if(int_radscal.eq.int_radscal_slater)then
           rad=atmrad(atomchrg(icntr))
       elseif (int_radscal .eq.3)then
           rad=xbsl(atomchrg(icntr))
       elseif (int_radscal .eq.4)then
            rad=TA(atomchrg(icntr))
       elseif (int_radscal .eq.5)then
           rad=multiEX(atomchrg(icntr))
c all molecular centers use the same integration radius
       else
           rad=1.d0
       endif

c If using rigid partitioning
      if(int_parttyp.eq.int_parttyp_rigid)then

        write(*,9000)
 9000   format(t3,'@OCT - rigid partitioning not yet supported.')
        call errex
        return

cSSSc Determine the radial limit of integration
cSSS
cSSS        rintmin=1.d10
cSSS        do 40 iangpt=1,maxangpts
cSSS          do 20 iatom=1,natoms
cSSSc           if(iatom.eq.icntr)goto 20
cSSS            rtmp(iatom)=0.d0
cSSS            do 10 i=1,3
cSSS              rtmp(iatom)=rtmp(iatom)+(atmvc(icntr,iatom,i)+
cSSS     &            rstep*gridxyz(i,iangpt,maxanggrd))**2
cSSS   10       continue
cSSS            rtmp(iatom)=sqrt(rtmp(iatom))
cSSS   20     continue
cSSS
cSSS          do 30 iatom=1,natoms
cSSS            if(iatom.eq.icntr)goto 30
cSSS            cosangl=(rtmp(icntr)**2 + rij(icntr,iatom)**2 -
cSSS     &          rtmp(iatom)**2)/
cSSS     &          (2.d0*rtmp(icntr)*rij(icntr,iatom))
cSSS            if(cosangl.gt.1.d-10)then
cSSS              rint(iatom,iangpt)=
cSSS     &            relsiz(icntr,iatom)*rij(icntr,iatom)/cosangl
cSSS            else
cSSS              rint(iatom,iangpt)=1.d10
cSSS            endif
cSSS            if(rint(iatom,iangpt).lt.rintmin)
cSSS     &          rintmin=rint(iatom,iangpt)
cSSS   30     continue
cSSS   40   continue
cSSS
cSSSc Determine radial point and weight
cSSS
cSSS        iset=0
cSSS        radpt(1)=rad / ((xnradp1-1.d0)**2)
cSSS        bndry=(2.d0*rad*xnradp1) / ((xnradp1-1.d0)**3)
cSSS        rwt(1,1)=bndry*radpt(1)**2
cSSS        do 140 iradpt=2,int_numradpts
cSSS          xiradpt=dble(iradpt)
cSSS          radpt(iradpt)=(rad*xiradpt**2) / ((xnradp1-xiradpt)**2)
cSSS          dis=(2.d0*rad*xnradp1*xiradpt) / ((xnradp1-xiradpt)**3)
cSSS          bndry2=bndry+dis
cSSS
cSSSc The radial boundary is less than the shortest radial line,
cSSSc determine the radial weight in the normal fashion
cSSS          if(bndry2.lt.rintmin)then
cSSS            rwt(iradpt,1)=dis*radpt(iradpt)**2
cSSS            bndry=bndry2
cSSS
cSSS          else
cSSS
cSSSc The radial boundary is greater than the shortest radial line, set
cSSSc the angular grid for a smooth transition to the top angular grid
cSSSc which is used throughout the interstitial region
cSSS
cSSS            if(iset.eq.1)goto 110
cSSS            radgrid(iradpt-1)=itopgrd
cSSS            do 50 iangpt=2,maxangpts
cSSS              rwt(iradpt-1,iangpt)=rwt(iradpt-1,1)
cSSS   50       continue
cSSS            icnt=2
cSSS            do 80 ii=1,itopgrd-1
cSSS              if(radgrid(iradpt-icnt).lt.itopgrd-ii)
cSSS     &            radgrid(iradpt-icnt)=itopgrd-ii
cSSS              do 60 iangpt=2,grdangpts(radgrid(iradpt-icnt))
cSSS                rwt(iradpt-icnt,iangpt)=rwt(iradpt-icnt,1)
cSSS   60         continue
cSSS              icnt=icnt+1
cSSS              if(radgrid(iradpt-icnt).lt.itopgrd-ii)
cSSS     &            radgrid(iradpt-icnt)=itopgrd-ii
cSSS              do 70 iangpt=2,grdangpts(radgrid(iradpt-icnt))
cSSS                rwt(iradpt-icnt,iangpt)=rwt(iradpt-icnt,1)
cSSS   70         continue
cSSS              icnt=icnt+1
cSSS   80       continue
cSSS
cSSSc Fill out the angular weights for the rest of the radial points
cSSS            do 100 ii=1,iradpt-icnt
cSSS              do 90 iangpt=2,grdangpts(radgrid(ii))
cSSS                rwt(ii,iangpt)=rwt(ii,1)
cSSS   90         continue
cSSS  100       continue
cSSS
cSSS            iset=1
cSSS
cSSS  110       continue
cSSS
cSSS            radgrid(iradpt)=itopgrd
cSSS
cSSSc Radial lines will now end at different lengths for different angles
cSSS            do 130 iangpt=1,grdangpts(itopgrd)
cSSS
cSSSc Determine if the radial boundary is less than the current radial line
cSSS
cSSS              rinttmp=1.d10
cSSS              iset2=0
cSSS              do 120 iatom=1,natoms
cSSS                if(iatom.eq.icntr)goto 120
cSSS                if(rint(iatom,iangpt).lt.bndry2)then
cSSS                  if(rint(iatom,iangpt).lt.rinttmp)then
cSSS                    rinttmp=rint(iatom,iangpt)
cSSS                    iset2=1
cSSS                  endif
cSSS                endif
cSSS  120         continue
cSSS
cSSSc The radial boundary is less than the current radial line,
cSSSc determine the radial weight in the normal fashion
cSSS
cSSS              if(iset2.eq.0)then
cSSS                rwt(iradpt,iangpt)=dis*radpt(iradpt)**2
cSSS              else
cSSS
cSSSc The current radial line is greater than the radial boundary, adjust
cSSSc the radial weight accordingly
cSSS
cSSS                dis2=rinttmp-bndry
cSSS                if(dis2.lt.0.d0)then
cSSS                  rwt(iradpt,iangpt)=0.d0
cSSS                else
cSSS                  rwt(iradpt,iangpt)=dis2*radpt(iradpt)**2
cSSS                endif
cSSS              endif
cSSS  130       continue
cSSS
cSSS            bndry=bndry2
cSSS
cSSS          endif
cSSS
cSSS  140   continue

      else

c For fuzzy partitioning (or atoms)

c For Handy radial integration

        if (int_radtyp.eq.int_radtyp_handy) then
          do 160 iradpt=1,int_numradpts
            xiradpt=dble(iradpt)
            radpt(iradpt)=(rad*xiradpt**2) / ((xnradp1-xiradpt)**2)
            rwt(iradpt,1)=(2.d0*rad**3*xnradp1*xiradpt**5) /
     &          ((xnradp1-xiradpt)**7)
            do 150 iangpt=2,grdangpts(radgrid(iradpt))
              rwt(iradpt,iangpt)=rwt(iradpt,1)
  150       continue
  160     continue

        

       elseif (int_radtyp.eq.3) then
           do 260 iradpt=1,int_numradpts
             diradpt=dble(iradpt)
             xiradpt=cos(diradpt*pi/xnradp1)
             radpt(iradpt)=(rad*(one+xiradpt))/(one-xiradpt)
             w1=(two*pi*rad**three)/(xnradp1)
             w2=(one+xiradpt)**FbyT
             w3=    (one-xiradpt)**SbyT
           rwt(iradpt,1)=w1*w2/w3
           do 250 iangpt=2,grdangpts(radgrid(iradpt))
            rwt(iradpt,iangpt)=rwt(iradpt,1)
  250      continue
  260      continue

 


       elseif (int_radtyp.eq.4) then
          do 660 iradpt=1,int_numradpts
             diradpt=dble(iradpt)
             xiradpt=cos(diradpt*pi/xnradp1)
             w1=-rad/DLOG(two)
             w2=(one+xiradpt)**alpha
             w3=DLOG((one-xiradpt)/two)
             radpt(iradpt)=w1*w2*w3
             w4=pi*(rad**three)/xnradp1
         w5=(one+xiradpt)**(three*alpha)/(DLOG(two))**three
            w6=DSQRT((one+xiradpt)/(one-xiradpt))
            w7=(DLOG( (one-xiradpt)/two))**two
            w8=DSQRT((one-xiradpt)/(one+xiradpt))
c            w9=(DLOG(two/(one-xiradpt)))**three
c            rwt(iradpt,1)=w4*w5*W9*(w6*w7+alpha*w8)

           w9=(DLOG((one-xiradpt)/two))**three
                 rwt(iradpt,1)=w4*w5*(w6*w7-alpha*w8*w9)
         
        do 650 iangpt=2,grdangpts(radgrid(iradpt))
            rwt(iradpt,iangpt)=rwt(iradpt,1)
  650      continue
  660      continue


         elseif (int_radtyp .eq.5)then
           if (int_numradpts .eq.23) then
         do i=1,23
             radpt(i)=-Rad*DLOG(x23i(i))
             rwt(i,1)=(a23i(i)/x23i(i))*Rad**three
            do 750 iangpt=2,grdangpts(radgrid(i))
               rwt(i,iangpt)=rwt(i,1)
  750       continue
         end do
           else if (int_numradpts .eq. 50) then
         do i=1,50
            radpt(i)=-Rad*DLOG(x50i(i))
            rwt(i,1)=(a50i(i)/x50i(i))*Rad**three
            do 850 iangpt=2,grdangpts(radgrid(i))
               rwt(i,iangpt)=rwt(i,1)
  850       continue     
         end do 
           else if (int_numradpts .eq. 45) then
         do i=1,45
             radpt(i)=-Rad*DLOG(x45i(i))
             rwt(i,1)=(a45i(i)/x45i(i))*Rad**three
             do 950 iangpt=2,grdangpts(radgrid(i))
                rwt(i,iangpt)=rwt(i,1)
  950        continue
         end do
           else if (int_numradpts .eq. 40) then
         do i=1,40
             radpt(i)=-Rad*DLOG(x40i(i))
             rwt(i,1)=(a40i(i)/x40i(i))*Rad**three
             do 1950 iangpt=2,grdangpts(radgrid(i))
                rwt(i,iangpt)=rwt(i,1)
 1950        continue
         end do
           else if (int_numradpts .eq. 35) then
         do i=1,35
             radpt(i)=-Rad*DLOG(x35i(i))
             rwt(i,1)=(a35i(i)/x35i(i))*Rad**three
             do 2950 iangpt=2,grdangpts(radgrid(i))
               rwt(i,iangpt)=rwt(i,1)
 2950        continue
         end do
            else if (int_numradpts .eq. 30) then
         do i=1,30
             radpt(i)=-Rad*DLOG(x30i(i))
             rwt(i,1)=(a30i(i)/x30i(i))*Rad**3
             do 3950 iangpt=2,grdangpts(radgrid(i))
                rwt(i,iangpt)=rwt(i,1)
 3950        continue
         end do
            end if
 
       else

c Determine inner radial integration (Gauss-Legendre)

          if(natoms.eq.1)then
            xinlngth=int_radlimit
            numradpts=int_numradpts
            goto 510
          endif

          numradpts=int_numradpts/2
          xinlngth=1.d10

          do 310 i=1,natoms
            if(icntr.eq.i) goto 310
            if(rij(icntr,i).lt.xinlngth) then
              jcntr=i
              xinlngth=rij(icntr,i)
            endif
  310     continue

          if(xinlngth.gt.3.d0)xinlngth=2.d0

c          xinlngth=0.5d0*xinlngth
          xinlngth=(bslrd(icntr)/(bslrd(icntr)+bslrd(jcntr)))
     &      *xinlngth

  510     continue

          xn=dble(numradpts)
          midpt=(numradpts+1)/2
          xmid=0.5d0*xinlngth
          xinl=0.5d0*xinlngth

          do 320 i=1,midpt
            xi=dble(i)
            zz=cos(pi*(xi-0.25d0)/(xn+0.5d0))

  330       continue

            p1=1.d0
            p2=0.d0

            do 340 j=1,numradpts
              xj=dble(j)
              p3=p2
              p2=p1
              p1=((2.d0*xj-1.d0)*zz*p2-(xj-1.d0)*p3)/xj
  340       continue

            pp=xn*(zz*p1-p2)/(zz*zz-1.d0)
            zz1=zz
            zz=zz1-p1/pp
            if(abs(zz-zz1).gt.eps)go to 330
            radpt(i)=xmid-xinl*zz
            radpt(numradpts+1-i)=xmid+xinl*zz
            rwttmp=2.d0*xinl/((1.d0-zz*zz)*pp*pp)
            rwt(i,1)=(radpt(i)**2)*rwttmp
            rwt(numradpts+1-i,1)=(radpt(numradpts+1-i)**2)*rwttmp
            do 350 iangpt=2,grdangpts(radgrid(i))
              rwt(i,iangpt)=rwt(i,1)
  350       continue
            do 360 iangpt=2,grdangpts(radgrid(numradpts+1-i))
              rwt(numradpts+1-i,iangpt)=rwt(numradpts+1-i,1)
  360       continue
  320     continue

          if(natoms.eq.1) goto 520

c Determine outer radial integration (Gauss-Legendre)

          xoulngth=11.d0

          xmid=0.5d0*(xinlngth+xoulngth)
          xinl=0.5d0*(xoulngth-xinlngth)

          do 420 i=1,midpt
            xi=dble(i)
            zz=cos(pi*(xi-0.25d0)/(xn+0.5d0))

  430       continue

            p1=1.d0
            p2=0.d0

            do 440 j=1,numradpts
              xj=dble(j)
              p3=p2
              p2=p1
              p1=((2.d0*xj-1.d0)*zz*p2-(xj-1.d0)*p3)/xj
  440       continue

            pp=xn*(zz*p1-p2)/(zz*zz-1.d0)
            zz1=zz
            zz=zz1-p1/pp
            if(abs(zz-zz1).gt.eps)go to 430
            radpt(numradpts+i)=xmid-xinl*zz
            radpt(2*numradpts+1-i)=xmid+xinl*zz
            rwttmp=2.d0*xinl/((1.d0-zz*zz)*pp*pp)
            rwt(numradpts+i,1)=(radpt(numradpts+i)**2)*rwttmp
            rwt(2*numradpts+1-i,1)=(radpt(2*numradpts+1-i)**2)*
     &          rwttmp
            do 450 iangpt=2,grdangpts(radgrid(numradpts+i))
              rwt(numradpts+i,iangpt)=rwt(numradpts+i,1)
  450       continue
            do 460 iangpt=2,grdangpts(radgrid(2*numradpts+1-i))
              rwt(2*numradpts+1-i,iangpt)=rwt(2*numradpts+1-i,1)
  460       continue
  420     continue

  520     continue

c Determine outer radial integration (Handy method)
c          xinradp1=dble(numradpts+1)
c          do 370 iradpt=1,numradpts
c            xiradpt=dble(iradpt)
c            radpt(iradpt+numradpts)=(xiradpt**2)/((xinradp1-xiradpt)**2)
c     &          + xinlngth
c            rwt(iradpt+numradpts,1)=((2.d0*xinradp1*xiradpt) /
c     &          ((xinradp1-xiradpt)**3))*radpt(iradpt+numradpts)**2
c            do 380 iangpt=2,grdangpts(radgrid(iradpt+numradpts))
c              rwt(iradpt+numradpts,iangpt)=rwt(iradpt+numradpts,1)
c 380        continue
c 370      continue

        endif

      endif

      call callstack_pop
      return
      end
