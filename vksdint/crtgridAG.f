










      subroutine crtgridAG(icntr,atomchrg,atmvc,rij,aij,cdnt,rsqrd,
     &    rrtmp,wtintr,atmwt,xx,gridwt,rwt,iangpt,iradpt,grid,
     &    dw,dp,dz,pp,dmmuji)

c This routine determines the cartesian coordinates for arbitrary grid
c point

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




c This contains the global string for identifying the current subroutine
c or function (provided the programmer set it).  cf. tools/callstack.F
c BE GOOD AND RESET CURR ON EXIT!

      character*64                callstack_curr,callstack_prev
      common /callstack_curr_com/ callstack_curr,callstack_prev
      save   /callstack_curr_com/


      integer icntr,atomchrg(natoms)
      double precision
     &    cdnt(natoms,3),atmvc(natoms,natoms,3),
     &    rsqrd(natoms),rij(natoms,natoms),
     &    aij(natoms,natoms),wtintr(natoms),
     &    rrtmp(natoms),atmwt,xx(3)

      integer id, i,iatom,jatom
      double precision
     &  a,xmmuji, wttot,zmuij,xmuij,fuzzy,cutij


      double precision
     &  th,th1

      double precision pi,dw(3,natoms),W,DP(3,natoms,natoms),
     &     DZ(3,natoms),
     &       Z, PP(natoms),swt,rwt(int_numradpts,maxangpts),
     &       gridwt(maxangpts,numgrid),g4,ddcutij,ddcutji,
     &  mmuji,one,
     &       dmmuji(3), dcutij,dcutji,dcoeff,B,F2,F4   

       integer iangpt,iradpt,ia,ja,ijatom,jjatom,iiatom,k,grid
       one=1.0d0
       th=1.0d-12
       th1=1.0d-08
       pi=acos(-1.d0)
       callstack_curr='CRTGRID'

       call dzero(dw,natoms*3)
       call  dzero(dz,natoms*3)
       call dzero(dp,natoms*natoms*3)
       do 180 iatom=1,natoms
        rsqrd(iatom)=0.d0
        do 170 i=1,3
          cdnt(iatom,i)=atmvc(icntr,iatom,i)+xx(i)
          rsqrd(iatom)=rsqrd(iatom)+cdnt(iatom,i)**2
  170   continue
        rrtmp(iatom)=dsqrt(rsqrd(iatom))
  180  continue

       if(int_parttyp.eq.rigid)then
        atmwt=1.d0
       else
        wttot=0.d0
        do 200 iatom=1,natoms
          if (atomchrg(iatom).eq.110) then
            wtintr(iatom)=0.d0
            goto 200
          else
            wtintr(iatom)=1.d0
          endif
          do 190 jatom=1,natoms
            if (atomchrg(jatom).eq.110) goto 190
            if (iatom.ne.jatom) then 
            zmuij=(rrtmp(iatom)-rrtmp(jatom))/
     &          rij(iatom,jatom)
            xmuij=zmuij +aij(iatom,jatom)*(1.d0-zmuij**2)
            f4=fuzzy(int_fuzzyiter,xmuij)
            cutij=0.5d0*(1.d0-f4)
            wtintr(iatom)=wtintr(iatom)*cutij
            end if  
  190  continue
          wttot=wttot+wtintr(iatom)
  200  continue
        atmwt=wtintr(icntr)/wttot
       endif

c       if(icntr.eq.2)then
c        write(*,*) 'irad=',iradpt,'ia=',iangpT,'icntr=',icntr
       do 611 ja=1,natoms
         if(abs(aij(ja,1)+1.0d0) .lt. 1.0d-05) then
              pp(ja)=0.0d0
              goto 611
          else
              pp(ja)=1.0d0
         end if
       do 612 ia=1,natoms
c         write(*,*) 'ia=',ia
          if(ia .eq. ja) goto 612
           IF(ABS(AIJ(JA,IA)-1.0D0).LT.1.0D-05) GOTO 612
          mmuji=(rrtmp(ja)-rrtmp(ia))/rij(ja,ia)
          xmmuji=mmuji+aij(ja,ia)*(1.0d0-mmuji*mmuji) 
          F4=xmmuji
          G4=one
            do k=1,int_fuzzyiter
               g4=g4*(one-f4*f4)
               f4=f4*(1.5d0-0.5D0*f4*f4)
            end do
            f2=0.5d00*f4
            dcutij=0.5d00+f2
            dcutji=0.5d00-f2
           pp(ja)=pp(ja)*dcutji
            if(ja .eq. icntr) goto 612
            id=ja
            b=mmuji/(rij(ja,ia))**2
            a=1.0d0-2.0d0*aij(ja,ia)*mmuji
      dmmuji(1)=a*(-cdnt(id,1)/(rrtmp(id)*rij(ja,ia))-B*atmvc(ja,ia,1)) 
      dmmuji(2)=a*(-cdnt(id,2)/(rrtmp(id)*rij(ja,ia))-B*atmvc(ja,ia,2))
      dmmuji(3)=a*(-cdnt(id,3)/(rrtmp(id)*rij(ja,ia))-B*atmvc(ja,ia,3))

c      write(*,*)'d',dmmuji(3),B,atmvc(ja,ia,3),cdnt(ja,3)/rrtmp(ja),ja,ia

               dcoeff=-81.0d0/32.0d0*g4
              if(abs(dcutij).gt. th) then
                  ddcutij=dcoeff/dcutij
                  dp(1,ia,ja)=-ddcutij*dmmuji(1)
                  dp(2,ia,ja)=-ddcutij*dmmuji(2)
                  dp(3,ia,ja)=-ddcutij*dmmuji(3)
c              write(*,*) 'dp=',dp(3,ia,ja),ia,ja
c        write(*,*) 'ij=',dcoeff,ddcutij,'a=',a,'b=',b,ia,ja
              else
                 if(abs(g4) .gt. th1) then
                    call errex
                 endif
              endif
              
              if(abs(dcutji).gt. th) then
                  ddcutji=dcoeff/dcutji
                  dp(1,ja,ja)=dp(1,ja,ja)+ddcutji*dmmuji(1)
                  dp(2,ja,ja)=dp(2,ja,ja)+ddcutji*dmmuji(2) 
                  dp(3,ja,ja)=dp(3,ja,ja)+ddcutji*dmmuji(3)
c             write(*,*) 'dpp=',dp(3,ja,ja)

              else
                  if(abs(g4).gt. th1) then
                    call errex
                  end if
              end if
  612  continue
  611  continue
             Z=0.d0
             do iiatom=1,natoms
                z=z+pp(iiatom)
             end do 
       swt=rwt(iradpt,iangpt)*gridwt(iangpt,grid)
c       write(*,*) 'swt=',swt

            W=Pp(icntr)/Z*swt
c       write(*,*)'w=',W,Z

c      write(*,*) '1=',dp(3,1,1),'2=',dp(3,1,2),'3=',dp(3,1,3)
c      write(*,*) '4=',dp(3,2,1),'5=',dp(3,2,2),'6=',dp(3,2,3)
c      write(*,*) '7=',dp(3,3,1),'8=',dp(3,3,2),'9=',dp(3,3,3)

          
      do 613 jjatom=1,natoms
        if(jjatom .eq. icntr) goto 613
        do  ijatom=1,natoms
        dz(1,jjatom)=dz(1,jjatom)+DP(1,ijatom,jjatom)*Pp(ijatom)
        dz(2,jjatom)=dz(2,jjatom)+DP(2,ijatom,jjatom)*Pp(ijatom)
        dz(3,jjatom)=dz(3,jjatom)+DP(3,ijatom,jjatom)*Pp(ijatom)
        end do
        dw(1,jjatom)=w*(dp(1,icntr,jjatom)-dz(1,jjatom)/z)
        dw(2,jjatom)=w*(dp(2,icntr,jjatom)-dz(2,jjatom)/z)
        dw(3,jjatom)=w*(dp(3,icntr,jjatom)-dz(3,jjatom)/z)

        dw(1,icntr)=dw(1,icntr)-dw(1,jjatom)
        dw(2,icntr)=dw(2,icntr)-dw(2,jjatom)
        dw(3,icntr)=dw(3,icntr)-dw(3,jjatom)
c        write(*,*) 'dw=',dw(1,icntr),dw(2,icntr),dw(3,icntr),icntr
c      write(*,*) 'dww=',dw(1,jjatom),dw(2,jjatom),dw(3,jjatom),jjatom
  613  continue
c       endif
       return
       end
