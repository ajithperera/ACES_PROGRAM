










      Subroutine Erd_sphchk(Zscr,Lzscr)
C
C Check of the spherical functions against OED itself: for every
C shell with L >= 2 the spherical-mode overlap of the shell with
C itself, turned into VMOL's harmonics (W^T S W, Erd_vsph), must equal
C Cv^T S' Cv, S' the Cartesian-mode overlap scaled to ACES AOs (Fac)
C and Cv VMOL's harmonics over those AOs; contraction block by block.
C Stops if it does not.
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

      Parameter (Lmax = 7, Mxyz = (Lmax+1)*(Lmax+2)/2, Mxpk = 5000,
     &           Mxf = 30*Mxyz)
      Dimension Zscr(Lzscr), Alp(Mxpk), Cc(Mxpk), Cv(Mxyz,Mxyz),
     &          W(2*Lmax+1,2*Lmax+1), Fac(Mxyz), Sc(Mxf*Mxf)
      Integer List(2), Icb(Mxpk), Ice(Mxpk), Iparv(3,Mxyz)
      Logical Lsph, Lcrt

      Lcrt = .False.
      Lsph = .True.
      Errmx = 0.0D0
      Liscr = Lzscr/2
      Lz    = Lzscr - Liscr
      Do Is = 1, Nshl
         L = Lshl(Is)
         If (L .Lt. 2) Goto 100
         List(1) = Is
         List(2) = Is
         Call Erd_pack(2,List,Alp,Nalp,Cc,Ncc,Icb,Ice,Ncsum,Mxpk)
         Ia = Ishatm(Is)
         X = Xyzatm(1,Ia)
         Y = Xyzatm(2,Ia)
         Z = Xyzatm(3,Ia)
         Nc = Nctsh(Is)
         Np = Nprsh(Is)
         Nk = (L+1)*(L+2)/2
         Ns = 2*L + 1
         Call Erd_cmpfac(L,Lcrt,Fac,N)
         Call Erd_vsph(L,Cv,Mxyz,W,2*Lmax+1,Iparv)

         Call Oed__memory_ovl_batch(Nalp,Ncc,Nc,Nc,Np,Np,L,L,X,Y,Z,
     &        X,Y,Z,Alp,Cc,Lcrt,Imin,Iopt,Jmin,Jopt)
         If (Imin .Gt. Liscr .Or. Jmin .Gt. Lz) Call Errex
         Call Oed__gener_ovl_batch(Liscr,Lz,Nalp,Ncc,Ncsum,Nc,Nc,Np,Np,
     &        L,L,X,Y,Z,X,Y,Z,Alp,Cc,Icb,Ice,Lcrt,.False.,
     &        Zscr(Lz+1),Nbc,Nfc,Zscr)
         N1 = Nk*Nc
         Do I = 1, N1*N1
            Sc(I) = Zscr(Nfc+I-1)
         Enddo
         Call Oed__memory_ovl_batch(Nalp,Ncc,Nc,Nc,Np,Np,L,L,X,Y,Z,
     &        X,Y,Z,Alp,Cc,Lsph,Imin,Iopt,Jmin,Jopt)
         If (Imin .Gt. Liscr .Or. Jmin .Gt. Lz) Call Errex
         Call Oed__gener_ovl_batch(Liscr,Lz,Nalp,Ncc,Ncsum,Nc,Nc,Np,Np,
     &        L,L,X,Y,Z,X,Y,Z,Alp,Cc,Icb,Ice,Lsph,.False.,
     &        Zscr(Lz+1),Nbs,Nfs,Zscr)
         N2 = Ns*Nc
         Do It = 1, Nc
         Do Iu = 1, Nc
            Do M = 1, Ns
            Do M2 = 1, Ns
               Sum = 0.0D0
               Do K = 1, Nk
               Do K2 = 1, Nk
                  Sum = Sum + Cv(K,M)*Cv(K2,M2)*Fac(K)*Fac(K2)*
     &                  Sc(K + (It-1)*Nk + N1*(K2-1 + (Iu-1)*Nk))
               Enddo
               Enddo
               Ssph = 0.0D0
               Do J = 1, Ns
               Do J2 = 1, Ns
                  Ssph = Ssph + W(J,M)*W(J2,M2)*Zscr(Nfs + J-1 +
     &                   (It-1)*Ns + N2*(J2-1 + (Iu-1)*Ns))
               Enddo
               Enddo
               Errmx = Max(Errmx,Dabs(Sum-Ssph))
            Enddo
            Enddo
         Enddo
         Enddo
 100     Continue
      Enddo

      Write(6,"(2a,1p,e10.2)") " @Erd_sphchk: spherical functions vs",
     &     " OED, largest overlap error", Errmx
      If (Errmx .Gt. 1.0D-10) Then
         Write(6,"(a)") " @Erd_sphchk: Cartesian -> spherical map wrong"
         Call Errex
      Endif

      Return
      End
