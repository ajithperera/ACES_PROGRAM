










      Subroutine Erd_so2e(Lu,Tla,Lfull,W1,W2,Lw,Zscr,Lzscr,Iscr,Liscr,
     &                    Asym,Nwr)
C
C Two-electron SO integrals from ERD, written in VMOL's format: unit
C Lu(1) gets the (ii|ii) class (IIII), Lu(2) (ii|jj) (IIJJ), Lu(3)
C (ij|ij) (IJIJ) and Lu(4) (ij|kl) (IJKL), all four irreps different.
C Each unique integral is written once, with p <= q, r <= s in its
C label (p,q|r,s); integrals smaller than Tla are dropped.
C
C The loop is over quartets of shell orbits (K1 K2|K3 K4), K1 >= K2,
C K3 >= K4, (K1 K2) >= (K3 K4). For each, ERD gives the AO batches of
C the member shells, and four quarter transformations (Erd_qtr) turn
C every member index into an SO index. Only totally symmetric SO
C quartets are kept.
C
C Lfull = .True. sums over all members of the first orbit, and
C returns in Asym the largest SO integral that symmetry forbids (it
C must vanish: a check of the SO basis). Otherwise only the
C representative of the first orbit is used and the sum is scaled by
C the orbit size, which is exact for the totally symmetric integrals.
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




      Parameter (Lnbuf = 600, Mxpk = 5000)
      Integer Lu(4), Nwr(4), Iscr(Liscr), Nut(4), Ibuf(Lnbuf,4)
      Integer List(4), Icb(Mxpk), Ice(Mxpk), Kq(4), Ns(4), Nf(4),
     &        Nm(4), Isr(4)
      Dimension W1(Lw), W2(Lw), Zscr(Lzscr), Buf(Lnbuf,4),
     &          Alp(Mxpk), Cc(Mxpk)
      Logical Lfull, Screen
      Data Screen /.True./

      Ipack(I,J,K,L) = Ior(Ior(Ior(I,Ishft(J,Ibitwd)),
     &                 Ishft(K,2*Ibitwd)),Ishft(L,3*Ibitwd))

      Asym = 0.0D0
      Do I = 1, 4
         Nut(I) = 0
         Nwr(I) = 0
      Enddo

      Do K1 = 1, Norb
      Do K2 = 1, K1
      Do K3 = 1, K1
         K4top = K3
         If (K3 .Eq. K1) K4top = K2
      Do K4 = 1, K4top
         Kq(1) = K1
         Kq(2) = K2
         Kq(3) = K3
         Kq(4) = K4
         Do I = 1, 4
            Isr(I) = Iorbsh(1,Kq(I))
            Nf(I)  = Nfnsh(Isr(I))
            Nm(I)  = Norbm(Kq(I))
         Enddo
         If (Lfull) Then
            M1hi = Nm(1)
         Else
            M1hi = 1
         Endif
         Ntot = Nf(1)*M1hi*Nf(2)*Nm(2)*Nf(3)*Nm(3)*Nf(4)*Nm(4)
         If (Ntot .Gt. Lw) Then
            Write(6,"(a,2i12)") " @Erd_so2e: work space too small ",
     &                          Ntot, Lw
            Call Errex
         Endif
         Nd1 = Nf(1)*M1hi
         Nd2 = Nf(2)*Nm(2)
         Nd3 = Nf(3)*Nm(3)
         Nd4 = Nf(4)*Nm(4)

C AO integrals over the member shells, stored as
C W1(f1,m1,f2,m2,f3,m3,f4,m4).

         Do M1 = 1, M1hi
         Do M2 = 1, Nm(2)
         Do M3 = 1, Nm(3)
         Do M4 = 1, Nm(4)
            List(1) = Iorbsh(M1,K1)
            List(2) = Iorbsh(M2,K2)
            List(3) = Iorbsh(M3,K3)
            List(4) = Iorbsh(M4,K4)
            Call Erd_eri(List,Zscr,Lzscr,Iscr,Liscr,Alp,Cc,Icb,Ice,
     &                   Mxpk,Screen,Nbatch,Nfirst)
            Do I4 = 1, Nf(4)
            Do I3 = 1, Nf(3)
            Do I2 = 1, Nf(2)
            Do I1 = 1, Nf(1)
               Iw = I1 + Nf(1)*(M1-1)
     &            + Nd1*(I2-1 + Nf(2)*(M2-1)
     &            + Nd2*(I3-1 + Nf(3)*(M3-1)
     &            + Nd3*(I4-1 + Nf(4)*(M4-1))))
               If (Nbatch .Gt. 0) Then
                  W1(Iw) = Zscr(Nfirst + I1-1 + Nf(1)*(I2-1 +
     &                     Nf(2)*(I3-1 + Nf(3)*(I4-1))))
               Else
                  W1(Iw) = 0.0D0
               Endif
            Enddo
            Enddo
            Enddo
            Enddo
         Enddo
         Enddo
         Enddo
         Enddo

C Members -> SOs for indices 2, 3 and 4. The first index is done
C while the integrals are picked out below: it only needs the signs
C of the members of orbit K1 that were computed (all of them, or the
C representative scaled by the orbit size).

         Call Erd_qtr(W1,W2,Nd1,Nf(2),Nm(2),Nd3*Nd4,K2,1,1.0D0)
         Call Erd_qtr(W2,W1,Nd1*Nd2,Nf(3),Nm(3),Nd4,K3,1,1.0D0)
         Call Erd_qtr(W1,W2,Nd1*Nd2*Nd3,Nf(4),Nm(4),1,K4,1,1.0D0)
         Wt = Dble(Nm(1))/Dble(M1hi)

C Pick out the unique, totally symmetric SO integrals. W2 is now
C (f1,m1,f2,j2,f3,j3,f4,j4), f running fastest.

         Do J4 = 1, Nm(4)
         Do I4 = 1, Nf(4)
            Is = Lorbso(Iorbso(K4) + (I4-1)*Nm(4) + J4)
         Do J3 = 1, Nm(3)
         Do I3 = 1, Nf(3)
            Ir = Lorbso(Iorbso(K3) + (I3-1)*Nm(3) + J3)
         Do J2 = 1, Nm(2)
         Do I2 = 1, Nf(2)
            Iq = Lorbso(Iorbso(K2) + (I2-1)*Nm(2) + J2)
            Ib = Nd1*(I2-1 + Nf(2)*(J2-1) + Nd2*(I3-1 + Nf(3)*(J3-1)
     &         + Nd3*(I4-1 + Nf(4)*(J4-1))))
         Do J1 = 1, Nm(1)
         Do I1 = 1, Nf(1)
            Ip = Lorbso(Iorbso(K1) + (I1-1)*Nm(1) + J1)
            X = 0.0D0
            Do M1 = 1, M1hi
               X = X + Dble(Isosgn(M1,Ip))*W2(Ib + I1 + Nf(1)*(M1-1))
            Enddo
            X = X*Wt
            Ixor = Ieor(Ieor(Isoirr(Ip)-1,Isoirr(Iq)-1),
     &                  Ieor(Isoirr(Ir)-1,Isoirr(Is)-1))
            If (Ixor .Ne. 0) Then
               Asym = Max(Asym,Dabs(X))
               Goto 10
            Endif
            If (K1 .Eq. K2 .And. Ip .Lt. Iq) Goto 10
            If (K3 .Eq. K4 .And. Ir .Lt. Is) Goto 10
            Ipq = Max(Ip,Iq)*(Max(Ip,Iq)-1)/2 + Min(Ip,Iq)
            Irs = Max(Ir,Is)*(Max(Ir,Is)-1)/2 + Min(Ir,Is)
            If (K1 .Eq. K3 .And. K2 .Eq. K4 .And. Ipq .Lt. Irs) Goto 10
            If (Dabs(X) .Lt. Tla) Goto 10

            I = Min(Ip,Iq)
            J = Max(Ip,Iq)
            K = Min(Ir,Is)
            L = Max(Ir,Is)
            Ii = Isoirr(I)
            Ij = Isoirr(J)
            Ik = Isoirr(K)
            Il = Isoirr(L)
            If (Ii .Eq. Ij .And. Ik .Eq. Il) Then
               If (Ii .Eq. Ik) Then
                  Ifl = 1
               Else
                  Ifl = 2
               Endif
            Else If (Ii .Eq. Ik .And. Ij .Eq. Il) Then
               Ifl = 3
            Else
               Ifl = 4
            Endif
            Nut(Ifl) = Nut(Ifl) + 1
            Buf(Nut(Ifl),Ifl)  = X
            Ibuf(Nut(Ifl),Ifl) = Ipack(I,J,K,L)
            Nwr(Ifl) = Nwr(Ifl) + 1
            If (Nut(Ifl) .Eq. Lnbuf) Then
               Write(Lu(Ifl)) (Buf(Ii,Ifl),Ii=1,Lnbuf),
     &                        (Ibuf(Ii,Ifl),Ii=1,Lnbuf), Nut(Ifl)
               Nut(Ifl) = 0
            Endif
 10         Continue
         Enddo
         Enddo
         Enddo
         Enddo
         Enddo
         Enddo
         Enddo
         Enddo
      Enddo
      Enddo
      Enddo
      Enddo

      Do Ifl = 1, 4
         If (Lu(Ifl) .Gt. 0) Then
            If (Nut(Ifl) .Ne. 0) Write(Lu(Ifl))
     &         (Buf(Ii,Ifl),Ii=1,Lnbuf), (Ibuf(Ii,Ifl),Ii=1,Lnbuf),
     &         Nut(Ifl)
            Nut(Ifl) = -1
            Write(Lu(Ifl)) (Buf(Ii,Ifl),Ii=1,Lnbuf),
     &                     (Ibuf(Ii,Ifl),Ii=1,Lnbuf), Nut(Ifl)
         Else If (Nwr(Ifl) .Ne. 0) Then
            Write(6,"(a,i2)") " @Erd_so2e: integrals for closed file",
     &                        Ifl
            Call Errex
         Endif
      Enddo

      Return
      End
