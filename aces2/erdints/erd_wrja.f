










      Subroutine Erd_wrja(Icmpao,Ncrt,U,Ui,Csy,Scr,Ipiv,Lwork)
C
C JOBARC records of the SO basis that vmol2ja (driver.F, v2ja.F)
C writes from VMOL's SYMTRANS and CSYMTRAN: NIRREP, NUMBASIR, NBASTOT,
C NAOBASFN, FAOBASIR, AO2SO (SO -> computational AO, rows are the
C Cartesian AOs; the name is the reverse of what it does), FULLSOAO,
C FULLAOSO, AO2SOINV, the same pair with the AOs in ZMAT order
C (CMP2ZMAT, ZMAT2CMP), and for spherical bases CART3CMP, CART2CMP,
C CMP3CART, CMP2CART.
C
C U = FULLSOAO is square: the Nso SOs, then (spherical basis only) the
C symmetry adapted contaminants r^2j Y_(L-2j) that complete the
C Cartesian space, all made of VMOL's functions (Erd_vharm). Csy (= VMOL's CSYMTRAN) is the purely Cartesian SO
C basis. As in vmol2ja the inverses are taken with LAPACK.
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




      Parameter (Lmax = 7, Mxyz = (Lmax+1)*(Lmax+2)/2)
      Integer Icmpao(Ncrt), Ipiv(Ncrt), Nfao(8), Irc(8), Isg(8,8),
     &        Nconir(8),
     &        Iparv(3,Mxyz), Ixyz(3,Mxyz), Ipar(3)
      Dimension U(Ncrt,Ncrt), Ui(Ncrt,Ncrt), Csy(Ncrt,Ncrt),
     &          Scr(Lwork), Cv(Mxyz,Mxyz), W(13,13)

      Call Putrec(20,'JOBARC','NIRREP  ',1,Nirr)
      Call Putrec(20,'JOBARC','NUMBASIR',Nirr,Nbfirr)
      Call Putrec(20,'JOBARC','NBASTOT ',1,Nso)
      Call Putrec(20,'JOBARC','NAOBASFN',1,Ncrt)

C SOs -> AOs.

      Call Dzero(U,Ncrt*Ncrt)
      Do Ip = 1, Nso
         Io  = Isoorb(Ip)
         Is0 = Iorbsh(1,Io)
         L   = Lshl(Is0)
         Nk  = (L+1)*(L+2)/2
         Nc  = Nfnsh(Is0)/Nctsh(Is0)
         It  = (Isofn(Ip)-1)/Nc + 1
         Mc  = Mod(Isofn(Ip)-1,Nc) + 1
         If (Spher .And. L .Ge. 2)
     &      Call Erd_vsph(L,Cv,Mxyz,W,13,Iparv)
         Do M = 1, Norbm(Io)
            I0 = Iofcrt(Iorbsh(M,Io)) + (It-1)*Nk
            If (Spher .And. L .Ge. 2) Then
               Do K = 1, Nk
                  U(Icmpao(I0+K),Ip) = Dble(Isosgn(M,Ip))*Cv(K,Mc)
               Enddo
            Else
               U(Icmpao(I0+Mc),Ip) = Dble(Isosgn(M,Ip))
            Endif
         Enddo
      Enddo

C Contaminants (irrep by irrep) and the Cartesian SO basis.

      Do Ig = 1, Nirr
         Nfao(Ig) = 0
      Enddo
      If (Spher) Then
         Ncol = Nso
         Do Ig = 1, Nirr
            Nconir(Ig) = 0
         Enddo
         Do Ig = 1, Nirr
            Do Io = 1, Norb
               Is0 = Iorbsh(1,Io)
               L   = Lshl(Is0)
               Nk  = (L+1)*(L+2)/2
               If (L .Ge. 2) Then
                  Call Erd_vsph(L,Cv,Mxyz,W,13,Iparv)
                  Ns = 2*L + 1
                  Do It = 1, Nctsh(Is0)
                     Do Ic = Ns+1, Nk
                        Call Erd_irrfn(Io,Iparv(1,Ic),Irc,Isg)
                        Do J = 1, Norbm(Io)
                           If (Irc(J) .Eq. Ig) Then
                              Ncol = Ncol + 1
                              Nconir(Ig) = Nconir(Ig) + 1
C VMOL's phase: the member that comes first in the computational
C order gets the + sign (as in Erd_soord).
                              Kmn = Ncrt + 1
                              Do M = 1, Norbm(Io)
                                 I0 = Iofcrt(Iorbsh(M,Io))
     &                              + (It-1)*Nk
                                 If (Icmpao(I0+1) .Lt. Kmn) Then
                                    Kmn = Icmpao(I0+1)
                                    Isg0 = Isg(M,J)
                                 Endif
                              Enddo
                              Do M = 1, Norbm(Io)
                                 I0 = Iofcrt(Iorbsh(M,Io))
     &                              + (It-1)*Nk
                                 Do K = 1, Nk
                                    U(Icmpao(I0+K),Ncol) =
     &                              Dble(Isg(M,J)*Isg0)*Cv(K,Ic)
                                 Enddo
                              Enddo
                           Endif
                        Enddo
                     Enddo
                  Enddo
               Endif
            Enddo
         Enddo
         If (Ncol .Ne. Ncrt) Then
            Write(6,"(a,2i6)") " @Erd_wrja: contaminant count error ",
     &                         Ncol, Ncrt
            Call Errex
         Endif

C VMOL lists the contaminants of an irrep by their first AO.

         Call Erd_consrt(U,Ncrt,Nso,Nconir,Nirr,Ipiv,Scr)

         Call Dzero(Csy,Ncrt*Ncrt)
         Ncol = 0
         Do Ig = 1, Nirr
            Do Io = 1, Norb
               Is0 = Iorbsh(1,Io)
               L   = Lshl(Is0)
               Call Erd_monos(L,Ixyz,Nk)
               Do It = 1, Nctsh(Is0)
                  Do K = 1, Nk
                     Do Kk = 1, 3
                        Ipar(Kk) = Mod(Ixyz(Kk,K),2)
                     Enddo
                     Call Erd_irrfn(Io,Ipar,Irc,Isg)
                     Do J = 1, Norbm(Io)
                        If (Irc(J) .Eq. Ig) Then
                           Ncol = Ncol + 1
                           Nfao(Ig) = Nfao(Ig) + 1
                           Do M = 1, Norbm(Io)
                              I0 = Iofcrt(Iorbsh(M,Io)) + (It-1)*Nk
                              Csy(Icmpao(I0+K),Ncol) =
     &                           Dble(Isg(M,J))
                           Enddo
                        Endif
                     Enddo
                  Enddo
               Enddo
            Enddo
         Enddo
      Endif
      If (Spher) Call Erd_csyord(Csy,Ncrt,Nfao,Nirr,Ipiv,Scr)
      Call Putrec(20,'JOBARC','FAOBASIR',Nirr,Nfao)

      Call Putrec(20,'JOBARC','AO2SO   ',Ncrt*Nso*Iintfp,U)
      Call Putrec(20,'JOBARC','FULLSOAO',Ncrt*Ncrt*Iintfp,U)
      Call Erd_minv(U,Ui,Ncrt,Ipiv,Scr,Lwork)
      Call Putrec(20,'JOBARC','FULLAOSO',Ncrt*Ncrt*Iintfp,Ui)
      Call Erd_blkcp(Ui,Ncrt,Nso,Ncrt,Scr)
      Call Putrec(20,'JOBARC','AO2SOINV',Nso*Ncrt*Iintfp,Scr)

C Cartesian SOs <-> computational SOs (spherical basis only).

      If (Spher) Then
         Call Xgemm('N','N',Ncrt,Ncrt,Ncrt,1.0D0,Ui,Ncrt,Csy,Ncrt,
     &              0.0D0,Scr,Ncrt)
         Call Putrec(20,'JOBARC','CART3CMP',Ncrt*Ncrt*Iintfp,Scr)
         Call Dcopy(Ncrt*Ncrt,Scr,1,Csy,1)
         Call Erd_blkcp(Csy,Ncrt,Nso,Ncrt,Scr)
         Call Putrec(20,'JOBARC','CART2CMP',Nso*Ncrt*Iintfp,Scr)
         Call Dcopy(Ncrt*Ncrt,Csy,1,Scr,1)
         Call Erd_minv(Scr,Csy,Ncrt,Ipiv,Ui,Ncrt*Ncrt)
         Call Putrec(20,'JOBARC','CMP3CART',Ncrt*Ncrt*Iintfp,Csy)
         Call Putrec(20,'JOBARC','CMP2CART',Ncrt*Nso*Iintfp,Csy)
         Call Erd_minv(U,Ui,Ncrt,Ipiv,Scr,Lwork)
      Endif

C Same transformations with the AOs in ZMAT order.

      Do Ip = 1, Nso
         Do Iao = 1, Ncrt
            Scr(Iao + (Ip-1)*Ncrt) = U(Icmpao(Iao),Ip)
         Enddo
      Enddo
      Call Putrec(20,'JOBARC','CMP2ZMAT',Ncrt*Nso*Iintfp,Scr)
      Do Iao = 1, Ncrt
         Do Ip = 1, Nso
            Scr(Ip + (Iao-1)*Nso) = Ui(Ip,Icmpao(Iao))
         Enddo
      Enddo
      Call Putrec(20,'JOBARC','ZMAT2CMP',Nso*Ncrt*Iintfp,Scr)

      Return
      End
