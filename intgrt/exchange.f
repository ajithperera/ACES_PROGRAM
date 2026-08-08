










      subroutine exchange(buf,ibuf,dens,dens1,
     &    cdens,imap,ehfkin,ehfnatr,ehfx,
     &    ehfcoul,ehar,escf,iflg)

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




c ilnbuf : the length of a buffer to read in chunks of integral files
c (There is no practical reason why this would be anything but 600.)

      integer           ilnbuf
      common /vmol_com/ ilnbuf
      save   /vmol_com/

      external sb_bd_vmol


c TODO: these need proper gfname lookup calls
      integer iiiiio, ijijio, iijjio, ijklio, vpoutio
      character*(*) iiiifil, ijijfil, iijjfil, ijklfil, vpoutfil
      parameter (iiiiio=10)
      parameter (ijijio=21)
      parameter (iijjio=22)
      parameter (ijklio=23)
      parameter (vpoutio=30)
      parameter (iiiifil='IIII')
      parameter (ijijfil='IJIJ')
      parameter (iijjfil='IIJJ')
      parameter (ijklfil='IJKL')
      parameter (vpoutfil='VPOUT')


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



c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end


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




      integer maxirrep,num2comb,max2comb
      parameter (maxirrep=8)
      parameter (num2comb=22)
      parameter (max2comb=25)

c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end
      integer        nirrep, numbasir(8),
     &               irpsz1(36),irpsz2(28),irpds1(36),irpds2(56),
     &               old_irpoff(9), irrorboff(9), dirprd(8,8),
     &               old_iwoff1(37), old_iwoff2(29),
     &               inewvc(maxbasfn), idxvec(maxbasfn),
     &               irrtrilen(9), irrtrioff(8),
     &               irrsqrlen(9), irrsqroff(8)
      common /symm2/ nirrep, numbasir,
     &               irpsz1,    irpsz2,    irpds1,    irpds2,
     &               old_irpoff,    irrorboff,    dirprd,
     &               old_iwoff1,     old_iwoff2,
     &               inewvc,           idxvec,
     &               irrtrilen,    irrtrioff,
     &               irrsqrlen,    irrsqroff
      save   /symm2/

      integer             occup(8,2),totocc(2),totocca,totoccb,
     &                    maxirrtri,maxirrsqr,irrtritot,irrsqrtot
      common /sym_ks_com/ occup,     totocc,   totocca,totoccb,
     &                    maxirrtri,maxirrsqr,irrtritot,irrsqrtot
      save   /sym_ks_com/




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



      integer ibuf(ilnbuf),imap(nbastot),iflg
      double precision
     &    buf(ilnbuf),dens(nbastot,nbastot,iuhf+1),
     &    cdens(nbastot,nbastot,iuhf+1),
     &    dens1(nbastot*(nbastot+1)/2,iuhf+1),
     &    ehfkin,ehfnatr,ehfx,ehfcoul,ehar,escf,ierr

      integer
     &    iupki,iupkj,iupkk,iupkl,i,icnt,roff,isym,j,int,
     &    nut,k,l,i23,i24,i25,ind,index,ispin,alpha,beta,unit(3),
     &    iunit,nunit
      double precision
     &    ecoulab,ecoulaa,ecoulbb,eexchaa,eexchbb,x,tmp1a,tmp2a,tmp2b,
     &    tmp3a,tmp3b,fnuc,field,dipnuc,ehf1e,factor

      integer indxt,irrep,ibot,iposp,length,nbasir,indx1,indx2,iposf

      character*8 label
      character*8 dum1,dum2,dum3,title
      character*32 junk
      character*20 files(3)
      logical exist

      indxt(i,j)=i+(j*(j-1))/2

c Statement functions for packing and unpacking indices.
      iupki(int)=iand(int,ialone)
      iupkj(int)=iand(ishft(int,-ibitwd),ialone)
      iupkk(int)=iand(ishft(int,-2*ibitwd),ialone)
      iupkl(int)=iand(ishft(int,-3*ibitwd),ialone)

      call callstack_push('EXCHANGE')
      factor=1.d0
      alpha=1
      beta=2
      if (iuhf.eq.0) then
        factor=2.d0
        beta=1
      endif

c Make sure that all integral files are present
      nunit=3
      if (c1symmet.eq.1) nunit=1
      unit(1)=iiiiio
      unit(2)=ijijio
      unit(3)=iijjio
      files(1)=iiiifil
      files(2)=ijijfil
      files(3)=iijjfil
      do 10 iunit=1,nunit
        inquire(file=files(iunit),exist=exist)
        if (.not.exist) then
          write (*,9200) files(iunit)
          call errex
        endif
   10 continue
 9200 format('@EXCHANGE-F, file ',a12,' not found')

c fill alpha density array for contraction with one electron integrals
      icnt=0
      roff=0
      do 710 isym=1,nirrep
        do 712 i=1,numbasir(isym)
          do 715 j=1,i
            icnt=icnt+1
            do 711 ispin=1,iuhf+1
              dens1(icnt,ispin)=dens(roff+j,roff+i,ispin)
              if (i.ne.j) dens1(icnt,ispin)=2.d0*dens1(icnt,ispin)
  711       continue
  715     continue
  712   continue
        roff=roff+numbasir(isym)
  710 continue

c multiply off diagonal elements of the alpha density matrix by two
c for coulomb contraction with two electron integrals

      do 55 ispin=1,iuhf+1
        do 50 i=1,nbastot
          do 40 j=1,nbastot
            cdens(j,i,ispin)=2.d0*dens(j,i,ispin)
   40     continue
          cdens(i,i,ispin)=dens(i,i,ispin)
   50   continue
   55 continue

c read integrals and contract them with the appropriate density
c matrix for one-electron, coulomb and exchange energies

      open(iiiiio,file=iiiifil,status='old',form='unformatted',
     &    access='sequential')

c one electron energy

      ehf1e=0.d0
      call locate(iiiiio,'ONEHAMIL')
   60 read(iiiiio)buf,ibuf,nut
      do 80 ispin=1,iuhf+1
        do 70 int=1,nut
          ehf1e=ehf1e+dens1(ibuf(int),ispin)*buf(int)
   70   continue
   80 continue
      if (nut.ne.-1) goto 60
      ehf1e=ehf1e*factor

c kinetic energy

      ehfkin=0.d0
      call locate(iiiiio,'KINETINT')
   90 read(iiiiio)buf,ibuf,nut
      do 110 ispin=1,iuhf+1
        do 100 int=1,nut
          ehfkin=ehfkin+dens1(ibuf(int),ispin)*buf(int)
  100   continue
  110 continue
      if (nut.ne.-1) goto 90
      ehfkin=ehfkin*factor

      close(iiiiio,status='keep')

c two electron energy

      ecoulab=0.d0
      ecoulaa=0.d0
      ecoulbb=0.d0
      eexchaa=0.d0
      eexchbb=0.d0

      do 400 iunit=1,nunit

        open(unit(iunit),file=files(iunit),status='old',
     &      form='unformatted',access='sequential')
        call locate(unit(iunit),'TWOELSUP')
    1   read(unit(iunit))buf,ibuf,nut
        do 7 int=1,nut
          i=iupki(ibuf(int))
          j=iupkj(ibuf(int))
          k=iupkk(ibuf(int))
          l=iupkl(ibuf(int))
          x=buf(int)

          tmp1a=0.d0
          tmp2a=0.d0
          tmp2b=0.d0
          tmp3a=0.d0
          tmp3b=0.d0

c coulomb energy
          if (iunit.eq.1 .or. iunit.eq.3) then
c           alpha beta
            tmp1a=x*cdens(j,i,alpha)*cdens(l,k,beta)
            if (i.ne.k .or. j.ne.l)
     &          tmp1a=tmp1a+x*cdens(l,k,alpha)*cdens(j,i,beta)
c           alpha alpha
            tmp2a=x*cdens(j,i,alpha)*cdens(l,k,alpha)
            tmp2b=x*cdens(j,i,beta)*cdens(l,k,beta)

c           off diagonal elements
            if (i.ne.k .or. j.ne.l) then
              tmp2a=2.d0*tmp2a
              tmp2b=2.d0*tmp2b
            endif
          endif

c exchange energy
          if (iunit.eq.1 .or. iunit.eq.2) then
c           alpha alpha
            tmp3a=x*dens(l,i,alpha)*dens(k,j,alpha)
            tmp3b=x*dens(l,i,beta)*dens(k,j,beta)

            if (i.ne.j .or. k.ne.l) then
              tmp3a=tmp3a+x*dens(k,i,alpha)*dens(l,j,alpha)
              tmp3b=tmp3b+x*dens(k,i,beta)*dens(l,j,beta)
            endif

            if (i.ne.j .and. k.ne.l) then
              tmp3a=2.d0*tmp3a
              tmp3b=2.d0*tmp3b
            endif

c           off diagonal elements
            if (i.ne.k .or. j.ne.l) then
              tmp3a=2.d0*tmp3a
              tmp3b=2.d0*tmp3b
            endif
          endif

          ecoulab=ecoulab+tmp1a
          ecoulaa=ecoulaa+tmp2a
          ecoulbb=ecoulbb+tmp2b
          eexchaa=eexchaa+tmp3a
          eexchbb=eexchbb+tmp3b

    7   continue
        if (nut.ne.-1) goto 1
        close(unit(iunit),status='keep')
  400 continue

      ehfcoul=ecoulab+0.5d0*(ecoulaa+ecoulbb)
      ehfx=-0.5d0*(eexchaa+eexchbb)

c if finite field then nuclear dipole correction
      fnuc=0.d0
      if(iflags(23).ne.0.or.iflags(24).ne.0.or.iflags(25).ne.0) then
        open(vpoutio,file=vpoutfil,form='unformatted',status='old')
        rewind(vpoutio)
c ----------------------------------------------------------------------
c makmap()
        length=(nbastot*(nbastot+1))/2
        call izero(imap,length)
        ibot=1
        iposp=0
        do irrep=1,nirrep
           nbasir=numbasir(irrep)
           do indx1=ibot,ibot+nbasir-1
              do indx2=ibot,indx1
                 iposp=iposp+1
                 iposf=indxt(indx2,indx1)
                 imap(iposf)=iposp
              end do
           end do
           ibot=ibot+nbasir
        end do
c ----------------------------------------------------------------------

  312   i23=iflags(23)
        i24=iflags(24)
        i25=iflags(25)
        if(iflags(23).ne.0)then
          label='     X  '
          ind=iflags(23)
          iflags(23)=0
        elseif(iflags(24).ne.0)then
          label='     Y  '
          ind=iflags(24)
          iflags(24)=0
        elseif(iflags(25).ne.0)then
          label='     Z  '
          ind=iflags(25)
          iflags(25)=0
        else
          goto 23
        endif
        field=dble(ind)*1.d-6
c ----------------------------------------------------------------------
c seeklb()
        rewind(vpoutio)
        title=' '
        do while (title.ne.label)
           read(vpoutio,end=120,err=121)dum1,dum2,dum3,title
           goto 122
  120      write(*,*) '@EXCHANGE: end-of-file on unit ',vpoutio
           call errex
  121      write(*,*) '@EXCHANGE: error reading from unit ',vpoutio
           call errex
  122      continue
        end do
c ----------------------------------------------------------------------
        backspace(vpoutio)
        read(vpoutio)junk,dipnuc
        fnuc=fnuc+field*dipnuc

  123   read(vpoutio)buf,ibuf,nut
        if (nut.eq.-1) goto 312
        do 313 index=1,nut
          roff=imap(ibuf(index))
          if (roff.ne.0)
     &        ehf1e=ehf1e+(dens1(roff,alpha)+
     &        dens1(roff,beta))*buf(index)*field
  313   continue
        goto 123

   23   close(vpoutio,status='keep')
        iflags(23)=i23
        iflags(24)=i24
        iflags(25)=i25
      endif

      ehfnatr=ehf1e-ehfkin
      escf=ehf1e+ehfcoul+ehfx+nucrep+fnuc
      if (int_ks) then

c     escf is the expectation value of H w.r.t. tke KS
c          single determinant.

           call putrec(1,'JOBARC','KSSCFENG',iintfp,escf)
      end if

      ehar=escf-ehfx

      if (int_printscf.gt.int_printlev) then
        if(iflg.eq.0)then
          write(*,9000) ehfnatr
          write(*,9010) ehfkin
          write(*,9030) ehfcoul
          write(*,9060) ehfx
          write(*,9020) ehf1e
          write(*,9080) escf
        else
          write(*,9100) ehfnatr
          write(*,9110) ehfkin
          write(*,9120) ehf1e
          write(*,9130) ehfcoul
          write(*,9160) ehfx
          write(*,9180) escf
        endif
      endif

 9000 format(/'The SCF nuclear-electron attraction energy =',f20.12)
 9010 format('The SCF kinetic energy                     =',f20.12)
 9020 format('The SCF one electron energy                =',f20.12)
 9030 format('The SCF coulomb energy                     =',f20.12)
 9060 format('The SCF exchange energy                    =',f20.12)
 9080 format('The SCF total energy                       =',f20.12)

 9100 format(/'The correlated nuclear-electron attraction energy =',
     &    f20.12)
 9110 format('The correlated kinetic energy                     =',
     &    f20.12)
 9120 format('The correlated one electron energy                =',
     &    f20.12)
 9130 format('The correlated coulomb energy                     =',
     &    f20.12)
 9160 format('The correlated exchange energy                    =',
     &    f20.12)
 9180 format('The correlated density total SCF energy           =',
     &    f20.12)

      call callstack_pop
      return
      end
