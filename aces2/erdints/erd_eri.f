










      Subroutine Erd_eri(List,Zscr,Lzscr,Iscr,Liscr,Alp,Cc,Icb,Ice,
     &                   Mxpk,Screen,Nbatch,Nfirst)
C
C One ERD batch (List(1) List(2)|List(3) List(4)) of contracted
C electron repulsion integrals, Zscr(Nfirst ...) with the functions of
C List(1) running fastest, over the basis functions the SOs are made
C of (Erd_fntr). The scratch needed is checked against
C ERD__MEMORY_ERI_BATCH every time.
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

      Integer List(4), Iscr(Liscr), Icb(Mxpk), Ice(Mxpk)
      Dimension Zscr(Lzscr), Alp(Mxpk), Cc(Mxpk), X(4), Y(4), Z(4)
      Logical Screen

      Call Erd_pack(4,List,Alp,Nalp,Cc,Ncc,Icb,Ice,Ncsum,Mxpk)
      Do I = 1, 4
         Ia = Ishatm(List(I))
         X(I) = Xyzatm(1,Ia)
         Y(I) = Xyzatm(2,Ia)
         Z(I) = Xyzatm(3,Ia)
      Enddo
      I1 = List(1)
      I2 = List(2)
      I3 = List(3)
      I4 = List(4)

      Call Erd__memory_eri_batch(Nalp,Ncc,
     &     Nctsh(I1),Nctsh(I2),Nctsh(I3),Nctsh(I4),
     &     Nprsh(I1),Nprsh(I2),Nprsh(I3),Nprsh(I4),
     &     Lshl(I1),Lshl(I2),Lshl(I3),Lshl(I4),
     &     X(1),Y(1),Z(1),X(2),Y(2),Z(2),X(3),Y(3),Z(3),X(4),Y(4),Z(4),
     &     Alp,Cc,Spher,Imin,Iopt,Jmin,Jopt)
      If (Imin .Gt. Liscr .Or. Jmin .Gt. Lzscr) Then
         Write(6,"(a,4i12)") " @Erd_eri: ERD scratch too small ",
     &                       Imin, Liscr, Jmin, Lzscr
         Call Errex
      Endif
      Imax = Max(Imin,Min(Iopt,Liscr))
      Jmax = Max(Jmin,Min(Jopt,Lzscr))

      Call Erd__gener_eri_batch(Imax,Jmax,Nalp,Ncc,Ncsum,
     &     Nctsh(I1),Nctsh(I2),Nctsh(I3),Nctsh(I4),
     &     Nprsh(I1),Nprsh(I2),Nprsh(I3),Nprsh(I4),
     &     Lshl(I1),Lshl(I2),Lshl(I3),Lshl(I4),
     &     X(1),Y(1),Z(1),X(2),Y(2),Z(2),X(3),Y(3),Z(3),X(4),Y(4),Z(4),
     &     Alp,Cc,Icb,Ice,Spher,Screen,Iscr,Nbatch,Nfirst,Zscr)

C ERD functions -> the basis functions of the SOs (VMOL's).

      If (Nbatch .Gt. 0) Then
         N1 = Nfnsh(I1)
         N2 = Nfnsh(I2)
         N3 = Nfnsh(I3)
         N4 = Nfnsh(I4)
         Call Erd_fntr(Zscr(Nfirst),1,N1,N2*N3*N4,I1)
         Call Erd_fntr(Zscr(Nfirst),N1,N2,N3*N4,I2)
         Call Erd_fntr(Zscr(Nfirst),N1*N2,N3,N4,I3)
         Call Erd_fntr(Zscr(Nfirst),N1*N2*N3,N4,1,I4)
      Endif

      Return
      End
