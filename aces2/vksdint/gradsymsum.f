










       subroutine gradsymsum(ncount,gradE,gradF,cord,pop,equ,
     &                       pmapatom,gradS)

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




        double precision
     &  CORD(3,natoms),zero,gradE(natoms,3,natoms),xtemp,
     &  gradF(natoms,3,natoms),gradS(natoms,3)


         integer igrid,xyz,ncount,iatom,equ(natoms,natoms)
         integer r,ri,i,ref,pop(ncount),iorbit,pmapatom(natoms)

          zero=0.0d0
       if ( natoms .ne. ncount) then
        do iorbit=1,ncount
         do i=0,pop(iorbit)-1
            ref=kscore(polist+iorbit-1)
            r=equ(ref,ref)
            ri=equ(ref,ref+i)
             do xyz=1,3
               gradF(r,xyz,r)=gradF(r,xyz,r)+gradE(r,xyz,ri)
             end do
         end do
        end do

c          write(*,*) ' grad='
c             do i=1,natoms
c           write(*,100) gradF(i,1,i),gradF(i,2,i),gradF(i,3,i)
c         end do
 100    format (3f16.10)

       do iorbit=1,ncount
         do i=1,pop(iorbit)-1
            ref=kscore(polist+iorbit-1)
            r=equ(ref,ref)
            ri=equ(ref,ref+i)
            do xyz=1,3
        if(cord(xyz,ri) .eq. zero .and. cord(xyz,r) .eq. zero) then
                   xtemp=1.d0  
            else if( cord(xyz,ri) .eq. zero) then
                    xtemp=zero
               else
                    xtemp=cord(xyz,r)/cord(xyz,ri)
               end if
            gradF(ri,xyz,ri)=gradF(r,xyz,r)*xtemp            
           end do
           end do
        end do

c       write(*,*) ' full grad'
c       do i=1,natoms
c           write(*,100) gradF(i,1,i),gradF(i,2,i),gradF(i,3,i)
c         end do

         do iatom=1,natoms
           do iorbit=1,ncount
              do i=1,pop(iorbit)-1
                ref=kscore(polist+iorbit-1)
                r=equ(ref,ref)
                ri=equ(ref,ref+i)
                if(iatom .eq. r  .OR.  iatom .eq. ri) then
                   pmapatom(iatom)=iorbit 
                end if
              end do
           end do
         end do

        do iatom=1,natoms
           do igrid=1,natoms
c              write(*,*) pmapatom(igrid),igrid,pmapatom(iatom),iatom
              if(pmapatom(igrid) .ne. pmapatom(iatom)) then
c              write(*,*) 'ig=',igrid,'ia=',iatom
                do xyz=1,3
                gradS(iatom,xyz)=gradS(iatom,xyz)+gradE(iatom,xyz,igrid) 
                end do
              end if
           end do
         end do

c          write(*,*)'chek'
c           do i=1,natoms
c           write(*,100) gradS(i,1),grads(i,2),gradS(i,3)
c         end do

        do xyz=1,3
        do iatom=1,natoms
           gradS(iatom,xyz)=gradS(iatom,xyz)+gradF(iatom,xyz,iatom)
        end do
        end do

        else
           do iatom =1 ,natoms
             do igrid=1,natoms
               do xyz=1,3
                gradS(iatom,xyz)=gradS(iatom,xyz)+gradE(iatom,xyz,igrid)
               end do
             end do
           end do
        end if 
 
        return
        end 






