











      subroutine dintproc(ksdint,density,v,gradient,
     &                    nao,nso,
     &                    dscr,iscr,
     &                    order)

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





      integer           iuhf
      common /iuhf_com/ iuhf
      save   /iuhf_com/




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



      integer
     &   nao, nso,
     &   order(nao), iscr(nao)

      double precision
     &   ksdint(nao,nao,3,iuhf+1),
     &   density(nao,nao,iuhf+1),
     &   v(nao,nao),
     &   dscr(nao,nao),
     &   gradient(natoms*3)

      integer
     &   n, xyz, i, j, offset, spin

      double precision  ddot
      real sdot

      call callstack_push('DINTPROC')

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c 1) Get the gradient.
c
c 2) Since we need the trace of the product PC'VC, transform the density
c    from SO/comp ordering to AO/zmat ordering once instead of
c    transforming every potential (of which there are 3*natoms).
c
c 3) Reorder ksdint from AO/angular to AO/zmat with the map KSSRTORD
c    which was written to JOBARC earlier.
c 
c 4) Separate the centers by looping over all atoms and extract from
c    ksdint only the elements that contribute to the particular
c    atom (loop index). Observe:
c
c       ksdint(mu,nu,q,s) = < mu | Vxc,s | d(nu)/dq >
c
c       v(mu,nu) =   < d(mu)/dq | Vxc,s | nu >
c                  + < mu | Vxc,s | d(nu)/dq >
c 
c 4) Add the dot product to the appropriate gradient element.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c   o get the gradient
      call getrec(1,'JOBARC','GRADIENT',natoms*3*iintfp,gradient)


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c   o get the density matrix
      call mkdens(dscr,density)

c   o if this is a UHF calculation and if # SOs is less than # AOs then
c     we need to stretch out the matrix.
      if ((iuhf.eq.1).and.(nso.lt.nao)) call stretch(density,nso,nao)

c   o get the SO/comp to AO/zmat map
      call getrec(1,'JOBARC','CMP2ZMAT',nao*nso*iintfp,dscr)

c   o create CPC' using v (currently unused) as tmp
      do spin=1,iuhf+1
         call dsymm('R','L',nao,nso,
     &               1.d0,density(1,1,spin),nso,
     &                    dscr,             nao,
     &               0.d0,v,                nao)
         call dgemm('N','T',nao,nao,nso,
     &               1.d0,v,                nao,
     &                    dscr,             nao,
     &               0.d0,density(1,1,spin),nao)
      end do

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c   o get the AO/zmat to AO/ang map
      call getrec(1,'JOBARC','KSSRTORD',nao,order)

c   o reverse from AO/ang to AO/zmat ordering
      do spin=1,iuhf+1
         do xyz=1,3
            call mat_reorder
     &         (2,1,nao,nao,ksdint(1,1,xyz,spin),dscr,iscr,order)
         end do
      end do

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c   o get the AO center assignments in zmat ordering
      call getrec(1,'JOBARC','CNTERBF0',nao,order)

c   o get the computational to ZMAT ordering map
      call getrec(1,'JOBARC','MAP2ZMAT',natoms,iscr)


      do n=1,natoms
         do spin=1,iuhf+1
            do xyz=1,3

c            o build v
               do j=1,nao
                  if (order(j).eq.iscr(n)) then
                     do i=1,nao
                        v(i,j)=ksdint(i,j,xyz,spin)
                     end do
                  else
                     do i=1,nao
                        v(i,j)=0.d0
                     end do
                  end if
               end do

c            o add the dot product to the gradient
               offset=3*n-3+xyz
               if (iuhf.eq.0) then
                  gradient(offset)=
     &            gradient(offset)
     &            +4*ddot(nao*nao,density(1,1,spin),1,v,1)
               else
                  gradient(offset)=
     &            gradient(offset)
     &            +2*ddot(nao*nao,density(1,1,spin),1,v,1)
               end if
            end do
         end do
      end do

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c   o What are we waiting for?
      call putrec(1,'JOBARC','GRADIENT',natoms*3*iintfp,gradient)

      call callstack_pop
      return
      end

