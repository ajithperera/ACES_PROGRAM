










      Subroutine Erg_gnsp(List,Fac,Gam,Bx,N,Ng,Irg,Ncol,Icol0,Lcol,
     &                    Ipx,Npr,Ioffg,Gso,W1,W2,Lw)
C
C Add Fac times the nonseparable two-particle density of one shell
C quartet, over ERD functions, to Gam (N1 x N2 x N3 x N4, first index
C fastest):
C
C    Gam(i,j,k,l) += Fac sum Bx(i,p) Bx(j,q) Bx(k,r) Bx(l,s) Gam(pq,rs)
C
C Gam(pq,rs) from Erg_gam over the Ng SOs of the lists, Bx their
C expansion in ERD functions. The functions of shell s only enter the
C Ncol(s) SOs Lcol(Icol0(s)+1 ...) (Erg_bcol). The SO block is
C gathered into W1 and transformed one index at a time.
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

      Integer List(4), Irg(Ng), Ncol(*), Icol0(*), Lcol(*),
     &        Ipx(Ng,Ng), Npr(8), Ioffg(8), Nf(4), Nb(4), Io(4), I0(4)
      Dimension Gam(*), Bx(N,Ng), Gso(*), W1(Lw), W2(Lw)

      Do I = 1, 4
         Is    = List(I)
         Nf(I) = Nfnsh(Is)
         Nb(I) = Ncol(Is)
         Io(I) = Icol0(Is)
         I0(I) = Ioffn(Is)
      Enddo
      Nso4 = Nb(1)*Nb(2)*Nb(3)*Nb(4)
      If (2*Nso4 .Gt. Lw) Then
         Write(6,"(a)") " @Erg_gnsp: SO quartet buffer too small"
         Call Errex
      Endif

C Gather the SO block, W1(p,q,r,s).

      Ind = 0
      Do L = 1, Nb(4)
         Is4 = Lcol(Io(4)+L)
         Do K = 1, Nb(3)
            Is3 = Lcol(Io(3)+K)
            Igrs = Ieor(Irg(Is3)-1,Irg(Is4)-1) + 1
            Irs  = Ipx(Is3,Is4)
            Do J = 1, Nb(2)
               Is2 = Lcol(Io(2)+J)
               Do I = 1, Nb(1)
                  Is1 = Lcol(Io(1)+I)
                  Ind = Ind + 1
                  Igpq = Ieor(Irg(Is1)-1,Irg(Is2)-1) + 1
                  If (Igpq .Eq. Igrs) Then
                     W1(Ind) = Gso(Ioffg(Igpq) +
     &                         (Ipx(Is1,Is2)-1)*Npr(Igpq) + Irs)
                  Else
                     W1(Ind) = 0.0D0
                  Endif
               Enddo
            Enddo
         Enddo
      Enddo

C B blocks of the four shells go to the end of W2 (rows I0+1..,
C columns the SOs they enter); the transformed blocks alternate between
C W1 and W2.
C  1: W2(i,qrs)   = B1(i,p) W1(p,qrs)

      Lb  = Lw
      Ib1 = Lb - Nf(1)*Nb(1) + 1
      Call Erg_bblk(Bx,N,Ng,Lcol,I0(1),Nf(1),Io(1),Nb(1),W2(Ib1))
      M = Nb(2)*Nb(3)*Nb(4)
      Call Dgemm('N','N',Nf(1),M,Nb(1),1.0D0,W2(Ib1),Nf(1),W1,Nb(1),
     &           0.0D0,W2,Nf(1))

C  2: W1(i,j,rs)  = W2(i,q,rs) B2(j,q)

      Ib2 = Lw - Nf(2)*Nb(2) + 1
      Call Erg_bblk(Bx,N,Ng,Lcol,I0(2),Nf(2),Io(2),Nb(2),W1(Ib2))
      M  = Nb(3)*Nb(4)
      Do Irs = 1, M
         Call Dgemm('N','T',Nf(1),Nf(2),Nb(2),1.0D0,
     &              W2(1+(Irs-1)*Nf(1)*Nb(2)),Nf(1),W1(Ib2),Nf(2),
     &              0.0D0,W2(1+Lw/2+(Irs-1)*Nf(1)*Nf(2)),Nf(1))
      Enddo
      Nij = Nf(1)*Nf(2)
      Call Dcopy(Nij*M,W2(1+Lw/2),1,W1,1)

C  3: W2(ij,k,s)  = W1(ij,r,s) B3(k,r)

      Ib3 = Lw - Nf(3)*Nb(3) + 1
      Call Erg_bblk(Bx,N,Ng,Lcol,I0(3),Nf(3),Io(3),Nb(3),W2(Ib3))
      Do Is = 1, Nb(4)
         Call Dgemm('N','T',Nij,Nf(3),Nb(3),1.0D0,
     &              W1(1+(Is-1)*Nij*Nb(3)),Nij,W2(Ib3),Nf(3),
     &              0.0D0,W2(1+(Is-1)*Nij*Nf(3)),Nij)
      Enddo

C  4: Gam(ijk,l) += Fac W2(ijk,s) B4(l,s)

      Ib4 = Lw - Nf(4)*Nb(4) + 1
      Call Erg_bblk(Bx,N,Ng,Lcol,I0(4),Nf(4),Io(4),Nb(4),W1(Ib4))
      Nijk = Nij*Nf(3)
      Call Dgemm('N','T',Nijk,Nf(4),Nb(4),Fac,W2,Nijk,W1(Ib4),Nf(4),
     &           1.0D0,Gam,Nijk)

      Return
      End
