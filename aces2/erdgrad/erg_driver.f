





























































































































































































































      Subroutine Erg_driver(Icore,Maxcor,Iuhf)
C
C xerdgrad driver. The basis, symmetry and SO basis are built exactly
C as xerdints builds them (Erd_rdmol, Erd_symgrp, Erd_sobld,
C Erd_cmpord, Erd_soord), so the SOs are those the SCF vectors refer
C to. The gradient is then
C
C    nuclear repulsion + one-electron (Erg_one) + two-electron (Erg_two)
C
C over the atoms of MOL_ERDOED, symmetrized (Erg_symz) and written to
C the GRADIENT record in the computational (VMOL) center order, as
C xvdint writes it.
C
      Implicit Double Precision (A-H,O-Z)
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
c erdints.com : begin
c
c Basis set, geometry and symmetry data for xerdints (integrals from
c Flocke's OED/ERD, symmetry adapted to the ACES II SO basis without
c VMOL or VMOL2JA). Filled by Erd_rdmol and Erd_symgrp.
c
c Atoms are those of MOL_ERDOED (all real centers, ZMAT order).
c   NATM          : number of atoms
c   ZATM(i)       : nuclear charge
c   XYZATM(3,i)   : coordinates (bohr, computational orientation)
c
c A shell is one angular momentum block on one atom, possibly with
c several (general) contractions, as ERD wants it.
c   NSHL          : number of shells
c   ISHATM(s)     : atom of shell s
c   LSHL(s)       : angular momentum
c   NPRSH(s)      : number of primitives
c   NCTSH(s)      : number of contractions
c   IOFALP(s)     : ALPHA(IOFALP(s)+1 ...) are the exponents
c   IOFCC(s)      : CCOEF(IOFCC(s)+1 ...) are the coefficients,
c                   NPRSH x NCTSH, column major
c   IOFCB(s)      : ICCBEG/ICCEND(IOFCB(s)+1 ...) are the first/last
c                   nonzero primitive of each contraction
c   NFNSH(s)      : number of functions (components x contractions)
c   IOFFN(s)      : offset of the first function of s in ERD order;
c                   within a shell the components run fastest
c                   (contraction major)
c   IOFCRT(s)     : offset of the first Cartesian function of s in
c                   the ZMAT AO order (= ERD Cartesian order)
c   NFNERD        : total number of functions in ERD order
c   NCART         : total number of Cartesian functions
c   SPHER         : spherical harmonics (ERD in spherical mode)
c
c Symmetry (D2h and subgroups; every operation is a sign triple).
c   NGEN          : number of generators (0-3)
c   NOPS          : number of operations, 2**NGEN
c   IOPSGN(k,r)   : sign of coordinate k under operation r (r=0 is E);
c                   bit g of r set means generator g+1 is included
c   IATMAP(i,r)   : image of atom i under operation r
c   ISHMAP(s,r)   : image of shell s under operation r
c
      integer mxshl, mxprmt, mxcct, mxctrt
      parameter (mxshl=5000, mxprmt=50000, mxcct=500000,
     &           mxctrt=20000)

      integer natm, nshl, ishatm, lshl, nprsh, nctsh, iofalp, iofcc,
     &        iofcb, nfnsh, ioffn, iofcrt, nfnerd, ncart, iccbeg,
     &        iccend
      double precision zatm, xyzatm, alpha, ccoef
      logical spher
      common /erdbas/ zatm(mxatms), xyzatm(3,mxatms),
     &                alpha(mxprmt), ccoef(mxcct),
     &                natm, nshl, ishatm(mxshl), lshl(mxshl),
     &                nprsh(mxshl), nctsh(mxshl), iofalp(mxshl),
     &                iofcc(mxshl), iofcb(mxshl), nfnsh(mxshl),
     &                ioffn(mxshl), iofcrt(mxshl), nfnerd, ncart,
     &                iccbeg(mxctrt), iccend(mxctrt), spher

c   CCRAW(i)      : the coefficients as read (MOL_ERDOED, i.e. the MOL
c                   file's), same layout as CCOEF; VMOL writes these
c                   for vprops (Erd_wrvp)
      double precision ccraw
      common /erdraw/ ccraw(mxcct)

      integer ngen, nops, iopsgn, iatmap, ishmap
      common /erdsym/ ngen, nops, iopsgn(3,0:7), iatmap(mxatms,0:7),
     &                ishmap(mxshl,0:7)

c SO basis (Erd_sobld). A shell orbit is a shell and its images; each
c SO is a +-1 combination of one function (same component and
c contraction) on every member shell of an orbit, i.e. VMOL's
c unnormalized SOs. SOs are numbered irrep by irrep.
c   NIRR          : number of irreps (= NOPS)
c   NSO           : number of SOs (= NFNERD)
c   NBFIRR(g)     : SOs in irrep g (1 based, ACES irrep g)
c   NORB          : number of shell orbits
c   NORBM(k)      : members of orbit k (orbit size, 1 2 4 or 8)
c   IORBSH(m,k)   : member shells, IORBSH(1,k) is the representative
c   IORBOP(m,k)   : an operation taking IORBSH(1,k) to IORBSH(m,k)
c   KORBSH(s)     : orbit of shell s
c   IORBSO(k)     : LORBSO(IORBSO(k)+1 ...) are the SOs of orbit k,
c                   NORBM(k)*NFNSH(rep) of them
c   ISOIRR(p)     : irrep of SO p (1 based)
c   ISOORB(p)     : orbit of SO p
c   ISOFN(p)      : function (component x contraction, ERD order)
c                   within the orbit's shells
c   ISOSGN(m,p)   : sign of member m in SO p
c   FNFAC(i)      : ERD function i times FNFAC(i) is the ACES (VMOL)
c                   AO; for Cartesian components x^a y^b z^c this is
c                   sqrt((2a-1)!!(2b-1)!!(2c-1)!!), 1 for spherical ones
c   WVSPH(m,j,L)  : spherical basis, L >= 2: VMOL harmonic j =
c                   sum_m WVSPH(m,j,L) ERD spherical function m
c                   (Erd_vsph). The basis functions SOs are made of are
c                   VMOL's: Erd_fntr turns ERD integrals into them.
      integer mxfn, mxorb
      parameter (mxfn=20000, mxorb=mxshl)
      integer nirr, nso, nbfirr, norb, norbm, iorbsh, iorbop, korbsh,
     &        iorbso, lorbso, isoirr, isoorb, isofn, isosgn
      common /erdso/ nirr, nso, nbfirr(8), norb, norbm(mxorb),
     &               iorbsh(8,mxorb), iorbop(8,mxorb), korbsh(mxshl),
     &               iorbso(mxorb), lorbso(mxfn), isoirr(mxfn),
     &               isoorb(mxfn), isofn(mxfn), isosgn(8,mxfn)
      double precision fnfac, wvsph
      common /erdsof/ fnfac(mxfn), wvsph(13,13,2:6)

c erdints.com : end


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

      Integer Icore(Maxcor), Npr(8), Ioffg(8), Nbg(8)
      Character*3 Kgen(3)
      Logical Lecp, Lcorr, Debug

      Lecp = (Iflags(71) .Ne. 0)
      If (Iflags(11) .Gt. 1) Then
         Write(6,"(a)") " @Erg_driver: only RHF and UHF references"
     &                  //" so far"
         Call Errex
      Endif
      If (Lecp) Then
         Write(6,"(a)") " @Erg_driver: no ECP gradients yet"
         Call Errex
      Endif

      Call Erd_rdmol(Kgen,Tla)
      Call Erd_symgrp(Kgen)
      Call Erd_sobld
      Write(6,"(a,i4,a,i5,a,i6,a,l2)") " @Erg_driver: atoms ", Natm,
     &     " shells ", Nshl, " functions ", Nfnerd, " spherical", Spher
      Write(6,"(a,i2,3(1x,a3))") " @Erg_driver: generators ", Ngen,
     &     (Kgen(I),I=1,Ngen)

C Computational order of the AOs and SOs (as in Erd_driver). The
C ZMAT center -> atom map (Iz2a) and computational center -> ZMAT
C center map (Imap) are kept at K000 and K010.

      Call Getrec(20,'JOBARC','NATOMS  ',1,Natoms)
      Ncrt = Ncart
      K000 = 1
      K010 = K000 + Natoms
      J000 = K010 + Natoms
      J010 = J000 + Ncrt
      J020 = J010 + Natoms
      J030 = J020 + 3*Natoms*Iintfp
      J040 = J030 + 3*8*Iintfp
      J050 = J040 + Natoms
      J060 = J050 + Natoms
      J070 = J060 + Natoms
      J080 = J070 + Ncrt
      J090 = J080 + Ncrt
      J100 = J090 + Ncrt
      J110 = J100 + Ncrt
      J120 = J110 + Natoms
      J130 = J120 + Natoms
      If (J130 .Gt. Maxcor) Call Insmem("Erg_driver",J130,Maxcor)
      Call Erd_cmpord(Icore(J000),Ncrt,Icore(J010),Natoms,
     &                Icore(J020),Icore(J030),Icore(J040),Icore(J050),
     &                Icore(J060),Icore(J070),Icore(J080),Icore(J090),
     &                Icore(J100),Icore(J110),Icore(J120))
      Do I = 0, Natoms-1
         Icore(K000+I) = Icore(J010+I)
         Icore(K010+I) = Icore(J040+I)
      Enddo
      If (J010+13*Nso .Gt. Maxcor) Call Insmem("Erg_driver",
     &                                         J010+13*Nso,Maxcor)
      Call Erd_soord(Icore(J000),Ncrt,Icore(J010),Icore(J010+Nso),
     &               Icore(J010+2*Nso),11*Nso)

C Densities over ERD functions (Erg_dens), kept with B to the end.

      Lcorr = (Iflags(2) .Ne. 0)
      N   = Nfnerd
      N2  = N*N*Iintfp
      I0b = J000
      I0ha = I0b  + N2
      I0hb = I0ha + N2
      I0xa = I0hb + N2
      I0xb = I0xa + N2
      I0t  = I0xb + N2
      I0w  = I0t  + N2
      I0c  = I0w  + N2
      I0o  = I0c  + N2
      I0e  = I0o  + N*Iintfp
      I0s  = I0e  + N*Iintfp
      Iend = I0s  + 2*N2
      If (Iend .Gt. Maxcor) Call Insmem("Erg_driver",Iend,Maxcor)
      Call Erg_bso(Icore(I0b),N)
      Call Erg_dens(Iuhf,Lcorr,N,Icore(I0b),Icore(I0ha),Icore(I0hb),
     &              Icore(I0xa),Icore(I0xb),Icore(I0t),Icore(I0w),
     &              Icore(I0c),Icore(I0o),Icore(I0e),Icore(I0s))

C Nonseparable two-particle density (correlated gradients only). The
C MOINTS lists are over the SOs of the SCF, or for spherical functions
C over the Cartesian SOs (T = CMP2CART takes the SCF SOs to them):
C their expansion in ERD functions is Bx = B, or B T+ with the
C pseudoinverse T+ = (T(T) T)**-1 T(T) (Gam(Cart) = T Gam T(T)).
C Then the SOs each shell enters (Erg_bcol), the SO Gamma (Erg_gam)
C and the buffers that transform one shell quartet of it (Erg_gnsp).

      Ng   = N
      I0bx = I0b
      Lgso = 0
      Lw   = 1
      I0rg = I0c
      I0nc = I0rg
      I0ic = I0rg
      I0lc = I0rg
      I0px = I0rg
      I0g2 = I0rg
      If (Lcorr) Then
         Do I = 1, Nirr
            Nbg(I) = Nbfirr(I)
         Enddo
         If (Spher) Then
            Call Getrec(20,'JOBARC','NAOBASFN',1,Ng)
            Call Getrec(20,'JOBARC','FAOBASIR',Nirr,Nbg)
         Endif
         I0rg = I0c
         I0bx = I0rg + Ng
         If (Spher) Then
            I0ct = I0bx + N*Ng*Iintfp
            I0m  = I0ct + Ng*N*Iintfp
            I0p  = I0m  + N*N*Iintfp
            I0iv = I0p  + N*Ng*Iintfp
            Iend = I0iv + N
            If (Iend .Gt. Maxcor) Call Insmem("Erg_driver",Iend,Maxcor)
            Call Getrec(20,'JOBARC','CMP2CART',Ng*N*Iintfp,Icore(I0ct))
            Call Erg_pinv(N,Ng,Icore(I0ct),Icore(I0m),Icore(I0p),
     &                    Icore(I0iv))
            Call Dgemm('N','N',N,Ng,N,1.0D0,Icore(I0b),N,Icore(I0p),N,
     &                 0.0D0,Icore(I0bx),N)
            I0nc = I0bx + N*Ng*Iintfp
         Else
            I0bx = I0b
            I0nc = I0rg + Ng
         Endif
         Ip = 0
         Do I = 1, Nirr
            Do J = 1, Nbg(I)
               Icore(I0rg+Ip) = I
               Ip = Ip + 1
            Enddo
         Enddo
         I0ic = I0nc + Nshl
         I0lc = I0ic + Nshl
         Call Erg_bcol(1,Icore(I0bx),N,Ng,Icore(I0nc),Icore(I0ic),
     &                 Icore(I0lc),Ltot)
         Call Erg_bcol(2,Icore(I0bx),N,Ng,Icore(I0nc),Icore(I0ic),
     &                 Icore(I0lc),Ltot)
         Mxcol = 0
         Do Is = 1, Nshl
            Mxcol = Max(Mxcol,Icore(I0nc+Is-1))
         Enddo
         Lw   = 2*Mxcol**4
         I0px = I0lc + Ltot
         I0g2 = I0px + Ng*Ng
         Call Erg_gam(1,Icore(I0px),Ng,Nirr,Nbg,Icore(I0rg),Npr,Ioffg,
     &                Lgso,Dum,Dum)
         Iend = I0g2 + (Lgso+Ng*Ng)*Iintfp
         If (Iend .Gt. Maxcor) Call Insmem("Erg_driver",Iend,Maxcor)
         Call Erg_gam(2,Icore(I0px),Ng,Nirr,Nbg,Icore(I0rg),Npr,Ioffg,
     &                Lgso,Icore(I0g2),Icore(I0g2+Lgso*Iintfp))
         Write(6,"(a,i12)") " @Erg_driver: SO Gamma elements stored ",
     &                      Lgso
      Endif
      I0w1 = I0g2 + Lgso*Iintfp
      I0w2 = I0w1 + Lw*Iintfp

C Gradient. Scratch for OED/ERD after the gradient arrays and the
C Gamma buffer of the largest shell quartet.

      Mxf  = 0
      Do Is = 1, Nshl
         Mxf = Max(Mxf,Nfnsh(Is))
      Enddo
      Lgam = Mxf**4
      I0g  = I0w2 + Lw*Iintfp
      I0gs = I0g  + 3*Natm*Iintfp
      I0g1 = I0gs + 3*Natm*Iintfp
      I0q  = I0g1 + 3*Natm*Iintfp
      I0z  = I0q  + Lgam*Iintfp
      Left = Maxcor - I0z + 1
      If (Left .Le. 0) Call Insmem("Erg_driver",I0z,Maxcor)
      Lzscr = (Left/2)/Iintfp
      I0i   = I0z + Lzscr*Iintfp
      Liscr = Maxcor - I0i + 1

C One-electron part and nuclear repulsion.

      Call Zero(Icore(I0g1),3*Natm)
      Call Erg_one(Icore(I0g1),Icore(I0t),Icore(I0w),N,
     &             Icore(I0z),Lzscr,Icore(I0i),Liscr)
      Call Erg_symz(Icore(I0g1),Icore(I0g))
      Call Erg_nucrep(Icore(I0g))
      Call Dcopy(3*Natm,Icore(I0g),1,Icore(I0g1),1)

C Two-electron part. With the file ERDGRAD.DEBUG the separable and
C nonseparable parts are also computed and printed separately.

      Cns = 0.0D0
      If (Lcorr) Cns = 4.0D0
      Inquire(File='ERDGRAD.DEBUG',Exist=Debug)
      If (Debug .And. Lcorr) Then
         Do Ipart = 1, 2
            Call Zero(Icore(I0g),3*Natm)
            Call Erg_two(Icore(I0g),Icore(I0ha),Icore(I0hb),
     &           Icore(I0xa),Icore(I0xb),N,Icore(I0z),Lzscr,
     &           Icore(I0i),Liscr,Icore(I0q),Lgam,Ipart.Eq.1,
     &           Dble(Ipart-1),Icore(I0bx),Ng,Icore(I0rg),
     &           Icore(I0nc),Icore(I0ic),Icore(I0lc),Icore(I0px),Npr,
     &           Ioffg,Icore(I0g2),Icore(I0w1),Icore(I0w2),Lw)
            Call Erg_symz(Icore(I0g),Icore(I0gs))
            If (Ipart .Eq. 1) Then
               Write(6,"(a)") " @Erg_driver: one-electron + nuclear"
               Call Erg_prgrd(Icore(I0g1),Natm)
               Write(6,"(a)") " @Erg_driver: separable two-electron"
            Else
               Write(6,"(a)") " @Erg_driver: nonseparable two-electron"
     &                        //" (factor 1)"
            Endif
            Call Erg_prgrd(Icore(I0gs),Natm)
         Enddo
      Endif

      Call Zero(Icore(I0g),3*Natm)
      Call Erg_two(Icore(I0g),Icore(I0ha),Icore(I0hb),Icore(I0xa),
     &             Icore(I0xb),N,Icore(I0z),Lzscr,Icore(I0i),Liscr,
     &             Icore(I0q),Lgam,.True.,Cns,Icore(I0bx),Ng,
     &             Icore(I0rg),Icore(I0nc),Icore(I0ic),Icore(I0lc),
     &             Icore(I0px),Npr,Ioffg,Icore(I0g2),Icore(I0w1),
     &             Icore(I0w2),Lw)
      Call Erg_symz(Icore(I0g),Icore(I0gs))
      Call Daxpy(3*Natm,1.0D0,Icore(I0g1),1,Icore(I0gs),1)

      Call Erg_wrgrd(Icore(I0gs),Natm,Natoms,Icore(K000),Icore(K010),
     &               Icore(I0z))

      Return
      End
