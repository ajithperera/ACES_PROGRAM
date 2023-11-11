










      subroutine integint(nshellatom,nshellprim,cdnt,alpha,
     & rsqrd,valprim,
     &   pcoeff,valmo,totwt,gradprim,gradmo,xnat,
     &   xocc,ksvxc,factor,ifctr,valao,valgradao,fourind,scrnn,xcoeff,
     &   coord,compmemb,angfct,order,init,iradpt,iangpt,
     &   times,ksdint,dograd)

c This module determines the MO, density, energies and
c matrix elements involving the KS potential at a given
c grid point. It features nonstandard way of calculating <i|Vxc|j>.
c If F[n]=int f (n,|gradn|) dr, then <i|Vxc|j>=<i|df/dn|j> +
c int (df/d|gradn) (1/|grad n|) (grad n).grad(i*j) dr.
c 12/25/98 SI

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






c ***NOTE*** This is a genuine (though not serious) limit on what Aces3 can do.
c     12 => s,p,d,f,g,h,i,j,k,l,m,n
      integer maxangshell
      parameter (maxangshell=12)




      integer           iuhf
      common /iuhf_com/ iuhf
      save   /iuhf_com/




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








c This common block keeps track of timing information. It is only required
c in the crapsi and crapso routines for the time being as only the starting
c and initial times are determined. In the future, it may be nice to find
c the time required to do different pieces of a calculation, and the data
c here is sufficiently flexible to do so.

c timein   the time of the first call to timer
c timenow  the time of the previous/current call to timer (timenow is set
c          to the time returned, so up to the actual call, timenow actually
c          contains the value from the previous call)
c timetot  total time elapsed since the first call to timer
c timenew  total time elapsed since the last call to timer

      double precision            timein,timenow,timetot,timenew
      common /timeinfo/ timein,timenow,timetot,timenew
      save /timeinfo/





      integer
     &    nshellatom(natoms),nshellprim(maxshell,natoms),
     &    compmemb(natoms),ifctr,kol,ii,jj,
     &    angfct(numangfct,3),order(naobasfn),iradpt,iangpt

      double precision
     &    cdnt(natoms,3),alpha(totprim),rsqrd(natoms),
     &    valprim(totprim),pcoeff(totprim,nbastot,2),
     &    valmo(nbastot,2),
     &    totwt,gradprim(totprim,3),
     &    gradmo(nbastot,3,2),xnat(nbastot),
     &    xocc(nbastot,2),ksvxc(naobasfn,naobasfn,iuhf+1),factor,
     &    valao(naobasfn),valgradao(naobasfn,3),
     &    fourind((nbastot*(nbastot+1)/2)*
     &            (nbastot*(nbastot+1)/2+1)/2,iuhf+1),
     &    scrnn((nbastot*(nbastot+1))/2),
     &    xcoeff(naobasfn*totprim),
     &    coord(3,natoms)
      double precision
     &   ksdint(naobasfn,naobasfn,3,iuhf+1)
      double precision pi

      logical init
      logical dograd

      double precision
     &    times(20),time

      integer
     &    iatom,ifct,i,iang,imomfct,ispin,xyz,l,m,n,naonao,j

      double precision
     &    xyzn(-2:maxangshell+2,3),denscomp,texp,roab(2),
     &    gradcomp(3,2),grad(3,2),twoa,xao(2),
     &    xxderao(2),xcderao(2),xcxc(2),
     &    comp(3),foura,func

      double precision
     &   xtmp(fun_num_exch,2),ctmp(fun_num_corr,2),
     &   xdertmp(fun_num_exch,2),cdertmp(fun_num_corr,2),
     &   vtmp(2),vxdertmp(2),vcdertmp(2),
     &   abs_grcomp(2),xxdercomp(2),xcdercomp(2),
     &   lyptmp(2)

      double precision
     &    tfkef,wkef,ddot


      call callstack_push('INTEG')

      call timer(1)
      time=timenow

      pi = acos(-1.d0)

c Initialize integral and functional values

      if ((.not.int_ks).and.init) then
         do i=1,fun_num_exch
            fun_exch(i)=0.d0
         end do
         do i=1,fun_num_corr
            fun_corr(i)=0.d0
         end do
         do i=1,fun_num_hyb
            fun_hyb(i)=0.d0
         end do
         totele = 0.d0
         etottf = 0.d0
         etotw  = 0.d0
         init   = .false.
      end if

      ifct=1
      do iatom=1,natoms
c        Angular momentum
c                         xyzn(i,x) = x**i   if i>=0
c                                   = 0      if i<0
         do xyz=1,3
            xyzn(-2,xyz)=0.d0
            xyzn(-1,xyz)=0.d0
            xyzn( 0,xyz)=1.d0
         end do
         do iang=1,maxangmom+2
            do xyz=1,3
               xyzn(iang,xyz)=xyzn(iang-1,xyz)*cdnt(iatom,xyz)
            end do
         end do

c Loop over angular momentum for each atom and find the density and
c value of certain factors for each function.
c
c Each function is of the form:
c   x**l y**m z**n exp(-ar**2)
c
c The gradient is of the form:
c   [ lx**(l-1) - 2ax**(l+1) ] y**m z**n exp(-ar**2)
         do iang=1,nshellatom(iatom)
            l=angfct(iang,1)
            m=angfct(iang,2)
            n=angfct(iang,3)
            comp(1) = xyzn(m+0,2)*xyzn(n+0,3)
            comp(2) = xyzn(l+0,1)*xyzn(n+0,3)
            comp(3) = xyzn(l+0,1)*xyzn(m+0,2)

c           Density 0 component
            denscomp= xyzn(l+0,1)*comp(1)
c           Gradient -1 component
            gradcomp(1,1) = xyzn(l-1,1)*comp(1)*l
            gradcomp(2,1) = xyzn(m-1,2)*comp(2)*m
            gradcomp(3,1) = xyzn(n-1,3)*comp(3)*n

c           Gradient +1 component
            gradcomp(1,2) = xyzn(l+1,1)*comp(1)
            gradcomp(2,2) = xyzn(m+1,2)*comp(2)
            gradcomp(3,2) = xyzn(n+1,3)*comp(3)

            do imomfct=1,nshellprim(iang,iatom)
               texp=dexp(-alpha(ifct)*rsqrd(iatom))
               twoa=2.d0*alpha(ifct)
               foura=twoa*twoa

               valprim(ifct)=denscomp*texp
 
               gradprim(ifct,1)=(gradcomp(1,1)-twoa*gradcomp(1,2))*texp
               gradprim(ifct,2)=(gradcomp(2,1)-twoa*gradcomp(2,2))*texp
               gradprim(ifct,3)=(gradcomp(3,1)-twoa*gradcomp(3,2))*texp

               ifct=ifct+1
            end do

c        end do iang=1,nshellatom(iatom)
         end do
c     end do iatom=1,natoms
      end do
      call timer(1)
      times(1)=times(1)+(timenow-time)
      time=timenow
        
      do ispin=1,iuhf+1
c      o contract functions with alpha/beta MO coefficients
         call xgemm('t','n',nbastot,1,totprim,
     &              1.d0, pcoeff(1,1,ispin),totprim,
     &                   valprim,          totprim,
     &              0.d0,valmo(1,ispin),   nbastot)

c      o add to get alpha/beta MO values
         call vecprd(valmo(1,ispin),valmo(1,ispin),xnat,nbastot)
c      o form the alpha/beta density at this point
         roab(ispin)=ddot(nbastot,xnat,1,xocc(1,ispin),1)
         if (roab(ispin).ge.thresh) then
            do xyz=1,3
c            o contract x/y/z gradients with alpha/beta coefficients
               call xgemm('t','n',nbastot,1,totprim,
     &                    1.d0, pcoeff(1,1,ispin),  totprim,
     &                         gradprim(1,xyz),    totprim,
     &                    0.d0,gradmo(1,xyz,ispin),nbastot)
c            o add to get alpha/beta MO x/y/z gradients at this point
               call vecprd(valmo(1,ispin),gradmo(1,xyz,ispin),
     &                     xnat,nbastot)
c            o form alpha/beta x/y/z density gradient at this point
               grad(xyz,ispin)=2.d0*
     &         ddot(nbastot,xnat,1,xocc(1,ispin),1)
            end do
         else
            roab(ispin)=0.d0
            do xyz=1,3
               grad(xyz,ispin)=roab(ispin)
            end do
         end if
c     end do ispin=1,iuhf+1
      end do

      if (iuhf.eq.0) then
         roab(2)=roab(1)
         do i=1,3
            grad(i,2)=grad(i,1)
         end do
      end if

      call timer(1)
      times(2)=times(2)+(timenow-time)
      time=timenow

c The common block
         
      ro     = roab(1)+roab(2)
      if (ro.lt.thresh) then
         call timer(1)
         times(3)=times(3)+(timenow-time)
         time=timenow
         call callstack_pop
         return
      end if

      roinv  = 1.d0/ro
      rom    = roab(1)-roab(2)
      roa    = roab(1)
      rob    = roab(2)
      rs     = (4.d0*pi*ro/3.d0)**(-1.d0/3.d0)
      zeta   = (roa-rob)/ro
      gradx  = grad(1,1) + grad(1,2)
      grady  = grad(2,1) + grad(2,2)
      gradz  = grad(3,1) + grad(3,2)
      gradxm = grad(1,1) - grad(1,2)
      gradym = grad(2,1) - grad(2,2)
      gradzm = grad(3,1) - grad(3,2)
      gro2   = gradx**2 + grady**2 + gradz**2
      gro    = sqrt(gro2)
      gro2a  = grad(1,1)**2 + grad(2,1)**2 + grad(3,1)**2
      gro2b  = grad(1,2)**2 + grad(2,2)**2 + grad(3,2)**2

      grdaa=grad(1,1)*grad(1,1)+grad(2,1)*grad(2,1)+grad(3,1)*grad(3,1)
      grdbb=grad(1,2)*grad(1,2)+grad(2,2)*grad(2,2)+grad(3,2)*grad(3,2)
      grdab=grad(1,1)*grad(1,2)+grad(2,1)*grad(2,2)+grad(3,1)*grad(3,2)

c Kohn-Sham potential
      if (int_ks) then

c Contract functions with coefficients and add to get AO values at this point.

       call xgemm('t','n',naobasfn,1,totprim,1.d0,xcoeff,totprim,
     &    valprim,totprim,0.d0,valao,naobasfn)

c Contract functions with coeff. and add to get the values of gradients of AO
c needed for computing matrix elements with the KS potential when integration
c by parts is used.

       do xyz=1,3
          call xgemm('t','n',naobasfn,1,totprim,1.d0,xcoeff,totprim,
     &    gradprim(1,xyz),totprim,0.d0,valgradao(1,xyz),naobasfn) 
       end do

       xtmp(fun_exch_hf,1)=0.d0
       xtmp(fun_exch_hf,2)=0.d0
       xdertmp(fun_exch_hf,1)=0.d0
       xdertmp(fun_exch_hf,2)=0.d0
       xdertmp(fun_exch_lda,1)=0.d0
       xdertmp(fun_exch_lda,2)=0.d0
       cdertmp(fun_corr_vwn,1)=0.d0
       cdertmp(fun_corr_vwn,2)=0.d0

       if (int_kspot.eq.fun_special) then

c Check what is needed and call it

          if (coef_pot_exch(fun_exch_lda).ne.0.d0) then
              call pot_exch_lda(xtmp(fun_exch_lda,1),
     &                          xtmp(fun_exch_lda,2))
          end if

          if (coef_pot_exch(fun_exch_becke).ne.0.d0) then
              call pot_exch_becke(xtmp(fun_exch_becke,1),
     &                            xdertmp(fun_exch_becke,1),
     &                            xtmp(fun_exch_becke,2),
     &                            xdertmp(fun_exch_becke,2))
          end if

          if (coef_pot_exch(fun_exch_pbe).ne.0.d0) then
              call pot_exch_pbe(xtmp(fun_exch_pbe,1),
     &                          xdertmp(fun_exch_pbe,1),
     &                          xtmp(fun_exch_pbe,2),
     &                          xdertmp(fun_exch_pbe,2))
          end if

          if (coef_pot_exch(fun_exch_pw91).ne.0.d0) then
              call pot_exch_pw91(xtmp(fun_exch_pw91,1),
     &                           xdertmp(fun_exch_pw91,1),
     &                           xtmp(fun_exch_pw91,2),
     &                           xdertmp(fun_exch_pw91,2))
          end if

          if (coef_pot_corr(fun_corr_vwn).ne.0.d0) then
               call pot_corr_vwn(ctmp(fun_corr_vwn,1),
     &                           ctmp(fun_corr_vwn,2))
          end if

          if (coef_pot_corr(fun_corr_lyp).ne.0.d0) then
                 call pot_corr_lyp(ctmp(fun_corr_lyp,1),
     &                             cdertmp(fun_corr_lyp,1),
     &                             lyptmp(1),
     &                             ctmp(fun_corr_lyp,2),
     &                             cdertmp(fun_corr_lyp,2),
     &                             lyptmp(2))
          end if 

          if (coef_pot_corr(fun_corr_pbe).ne.0.d0) then
                 call pot_corr_pbe(ctmp(fun_corr_pbe,1),
     &                             cdertmp(fun_corr_pbe,1),
     &                             ctmp(fun_corr_pbe,2),
     &                             cdertmp(fun_corr_pbe,2)) 
          end if                                 


          if (coef_pot_corr(fun_corr_pw91).ne.0.d0) then
               call pot_corr_pw91(ctmp(fun_corr_pw91,1),
     &                            cdertmp(fun_corr_pw91,1),
     &                            ctmp(fun_corr_pw91,2),
     &                            cdertmp(fun_corr_pw91,2))    
          end if

       else 

          if (int_kspot.eq.fun_hyb_b3lyp) then

c We know exactly what we need and call the components

          call pot_corr_vwn(ctmp(fun_corr_vwn,1),
     &                      ctmp(fun_corr_vwn,2))
          call pot_corr_lyp(ctmp(fun_corr_lyp,1),
     &                      cdertmp(fun_corr_lyp,1),
     &                      lyptmp(1), 
     &                      ctmp(fun_corr_lyp,2),
     &                      cdertmp(fun_corr_lyp,2),
     &                      lyptmp(2))

          call pot_exch_lda(xtmp(fun_exch_lda,1),
     &                      xtmp(fun_exch_lda,2))
           call pot_exch_becke(xtmp(fun_exch_becke,1),
     &                         xdertmp(fun_exch_becke,1),
     &                         xtmp(fun_exch_becke,2),
     &                         xdertmp(fun_exch_becke,2))
          end if
       end if

        abs_grcomp(1)=dsqrt(gro2a)
        abs_grcomp(2)=dsqrt(gro2b)

       do 500 ispin=1,iuhf+1
         if (roab(ispin).lt.thresh) go to 500

         vtmp(ispin)=
     &        ddot(fun_num_exch,xtmp(1,ispin),1,coef_pot_exch,1)+
     &        ddot(fun_num_corr,ctmp(1,ispin),1,coef_pot_corr,1)

         vxdertmp(ispin)=
     &        ddot(fun_num_exch,xdertmp(1,ispin),1,coef_pot_exch,1)
     &        + lyptmp(ispin)*coef_pot_corr(fun_corr_lyp)

         vcdertmp(ispin)=
     &        ddot(fun_num_corr,cdertmp(1,ispin),1,coef_pot_corr,1)

         xao(ispin)=totwt*factor*vtmp(ispin)
         xxderao(ispin)=totwt*factor*vxdertmp(ispin)
         xcderao(ispin)=totwt*factor*vcdertmp(ispin)
         call dsyr('l',naobasfn,xao(ispin),valao,1,
     &                                 ksvxc(1,1,ispin),naobasfn)

         do xyz=1,3
            xxdercomp(ispin)=
     &             xxderao(ispin)*grad(xyz,ispin)/abs_grcomp(ispin)

            xcdercomp(ispin)=
     &             xcderao(ispin)*(grad(xyz,1)+grad(xyz,2))/gro

            xcxc(ispin)=xxdercomp(ispin)+xcdercomp(ispin)

            call dsyr2('l',naobasfn,xcxc(ispin),valao,1,
     &             valgradao(1,xyz),1,ksvxc(1,1,ispin),naobasfn)

         end do 
         if (dograd) then
            do xyz=1,3
               call dgemm('n','t',naobasfn,naobasfn,1,
     &                     -xao(ispin),valao,                naobasfn,
     &                                 valgradao(1,xyz),     naobasfn,
     &                     1.d0,        ksdint(1,1,xyz,ispin),naobasfn)
            end do 
         end if

 500  continue
 
      else

       if (int_tdks) then
          call xgemm('t','n',naobasfn,1,totprim,1.d0,xcoeff,totprim,
     &         valprim,totprim,0.d0,valao,naobasfn)

c Reorder and transform to SA AO at each point of the grid. Slow but
c memory more efficient. Note that we have evaluated the MOS at each
c point of the grid in order to calculate the density. We could directly
c get the integrals in MO basis set and eliminate the transformations in
c vtran. It has to be checked and investigated further!!!
 
          do i=1,naobasfn
              scrnn(i)=valao(order(i))
          end do

          call xgemm('t','n',nbastot,1,naobasfn,1.d0,ksvxc,naobasfn,
     &         scrnn,naobasfn,0.d0,valao,nbastot)

c Transformations done. Calculate only unique combinations.

          naonao=(nbastot*(nbastot+1))/2

          m=1
          do i=1,nbastot
             do j=i,nbastot
               scrnn(m)=valao(i)*valao(j)
               m=m+1
             end do
          end do

          do ispin=1,iuhf+1
               call ker_exch_lda(vxdertmp(1),vxdertmp(2))
               call ker_corr_vwn(ro,vcdertmp(1),vcdertmp(2))
               xao(ispin)=totwt*factor*(vxdertmp(ispin)+vcdertmp(ispin))
c                xao(ispin)=totwt*factor

               m=1
               do i=1,naonao
                  do j=i,naonao
                     fourind(m,ispin)=xao(ispin)*scrnn(i)*scrnn(j)
                     m=m+1
                  end do
               end do
 
           end do
        end if

       totele= totele+totwt*ro

       call func_corr_vwn(func)
       fun_corr(fun_corr_vwn) = fun_corr(fun_corr_vwn) + totwt*func

       call func_corr_lyp(func)
       fun_corr(fun_corr_lyp) = fun_corr(fun_corr_lyp) + totwt*func

       call func_corr_pbe(func)
       fun_corr(fun_corr_pbe) = fun_corr(fun_corr_pbe) + totwt*func

       call func_corr_pw91(func)
       fun_corr(fun_corr_pw91) = fun_corr(fun_corr_pw91) + totwt*func

       call func_corr_wl(func)
       fun_corr(fun_corr_wl) = fun_corr(fun_corr_wl) + totwt*func

       call func_corr_wi(func)
       fun_corr(fun_corr_wi) = fun_corr(fun_corr_wi) + totwt*func

       call func_exch_lda(func)
       fun_exch(fun_exch_lda) = fun_exch(fun_exch_lda) + totwt*func

       call func_exch_becke(func)
       fun_exch(fun_exch_becke) = fun_exch(fun_exch_becke) + totwt*func

       call func_exch_pbe(func)
       fun_exch(fun_exch_pbe) = fun_exch(fun_exch_pbe) + totwt*func

       call func_exch_pw91(func)
       fun_exch(fun_exch_pw91) = fun_exch(fun_exch_pw91) + totwt*func

       etottf= etottf+totwt*tfkef(roab(1),roab(2),grdaa,grdab,grdbb)
       etotw = etotw +totwt*wkef(roab(1),roab(2),grdaa,grdab,grdbb)

      end if
      call timer(1)
      times(3)=times(3)+(timenow-time)
      time=timenow
      call callstack_pop
      return

      end
