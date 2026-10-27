










      Subroutine Erg_two(G,Dt,Da,Db,N,Zscr,Lzscr,Iscr,Liscr,Gam,Lgam)
C
C Two-electron part of the gradient for a single determinant,
C
C    dE/dX += 1/2 sum Gam(mu,nu,la,si) d(mu nu|la si)/dX,
C
C    Gam = D(mu,nu) D(la,si)
C        - 1/2 sum_spin [Ds(mu,la) Ds(nu,si) + Ds(mu,si) Ds(nu,la)],
C
C (D = Da + Db), over the symmetry unique canonical shell quartets
C only. Each quartet carries its permutational degeneracy and its
C orbit size; the result is symmetrized later (Erg_symz). The
C derivative is taken with respect to every atom that carries a
C center of the quartet but the last one: all centers on that atom
C are differentiated together, and the last atom gets minus the sum
C (translational invariance). One-center quartets are skipped.
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

      Parameter (Mxpk = 5000)
      Dimension G(3,Natm), Dt(N,N), Da(N,N), Db(N,N), Zscr(Lzscr),
     &          Gam(Lgam)
      Integer Iscr(Liscr)
      Integer List(4), Icb(Mxpk), Ice(Mxpk), Iat(4), Iatd(4),
     &        Ider(3,4)
      Dimension Alp(Mxpk), Cc(Mxpk), X(4), Y(4), Z(4), Gq(3)
      Logical Screen, Lrep
      Data Screen /.False./

      Do Is1 = 1, Nshl
      Do Is2 = 1, Is1
      Ij = Is1*(Is1-1)/2 + Is2
      Do Is3 = 1, Is1
      Do Is4 = 1, Is3
      Kl = Is3*(Is3-1)/2 + Is4
      If (Kl .Gt. Ij) Goto 100
         List(1) = Is1
         List(2) = Is2
         List(3) = Is3
         List(4) = Is4

C Distinct atoms of the quartet.

         Nd = 0
         Do I = 1, 4
            Iat(I) = Ishatm(List(I))
            Lnew = 1
            Do J = 1, Nd
               If (Iatd(J) .Eq. Iat(I)) Lnew = 0
            Enddo
            If (Lnew .Eq. 1) Then
               Nd = Nd + 1
               Iatd(Nd) = Iat(I)
            Endif
         Enddo
         If (Nd .Eq. 1) Goto 100

         Call Erg_orbit(List,4,Lrep,Nimg)
         If (.Not. Lrep) Goto 100
         Fac = 0.5D0*Dble(Nimg)
         If (Is1 .Ne. Is2) Fac = 2.0D0*Fac
         If (Is3 .Ne. Is4) Fac = 2.0D0*Fac
         If (Ij  .Ne. Kl ) Fac = 2.0D0*Fac

         N1 = Nfnsh(Is1)
         N2 = Nfnsh(Is2)
         N3 = Nfnsh(Is3)
         N4 = Nfnsh(Is4)
         If (N1*N2*N3*N4 .Gt. Lgam) Then
            Write(6,"(a)") " @Erg_two: Gamma buffer too small"
            Call Errex
         Endif
         I0 = Ioffn(Is1)
         J0 = Ioffn(Is2)
         K0 = Ioffn(Is3)
         L0 = Ioffn(Is4)
         Ind = 0
         Do L = L0+1, L0+N4
         Do K = K0+1, K0+N3
         Do J = J0+1, J0+N2
         Do I = I0+1, I0+N1
            Ind = Ind + 1
            Gam(Ind) = Dt(I,J)*Dt(K,L) - 0.5D0*(
     &                 Da(I,K)*Da(J,L) + Da(I,L)*Da(J,K) +
     &                 Db(I,K)*Db(J,L) + Db(I,L)*Db(J,K))
         Enddo
         Enddo
         Enddo
         Enddo

         Call Erd_pack(4,List,Alp,Nalp,Cc,Ncc,Icb,Ice,Ncsum,Mxpk)
         Do I = 1, 4
            X(I) = Xyzatm(1,Iat(I))
            Y(I) = Xyzatm(2,Iat(I))
            Z(I) = Xyzatm(3,Iat(I))
         Enddo

         Do Ia = 1, Nd-1
            Do Kc = 1, 3
               Do I = 1, 4
                  Do K = 1, 3
                     Ider(K,I) = 0
                  Enddo
                  If (Iat(I) .Eq. Iatd(Ia)) Ider(Kc,I) = 1
               Enddo
               Call Erd__memory_eri_derv_batch(Nalp,Ncc,
     &              Nctsh(Is1),Nctsh(Is2),Nctsh(Is3),Nctsh(Is4),
     &              Nprsh(Is1),Nprsh(Is2),Nprsh(Is3),Nprsh(Is4),
     &              Lshl(Is1),Lshl(Is2),Lshl(Is3),Lshl(Is4),
     &              X(1),Y(1),Z(1),X(2),Y(2),Z(2),
     &              X(3),Y(3),Z(3),X(4),Y(4),Z(4),
     &              Ider(1,1),Ider(2,1),Ider(3,1),
     &              Ider(1,2),Ider(2,2),Ider(3,2),
     &              Ider(1,3),Ider(2,3),Ider(3,3),
     &              Ider(1,4),Ider(2,4),Ider(3,4),
     &              Alp,Cc,Spher,Imin,Iopt,Jmin,Jopt)
               Call Erg_scrchk("Erg_two",Imin,Liscr,Jmin,Lzscr)
               Imax = Max(Imin,Min(Iopt,Liscr))
               Jmax = Max(Jmin,Min(Jopt,Lzscr))
               Call Erd__gener_eri_derv_batch(Imax,Jmax,Nalp,Ncc,Ncsum,
     &              Nctsh(Is1),Nctsh(Is2),Nctsh(Is3),Nctsh(Is4),
     &              Nprsh(Is1),Nprsh(Is2),Nprsh(Is3),Nprsh(Is4),
     &              Lshl(Is1),Lshl(Is2),Lshl(Is3),Lshl(Is4),
     &              X(1),Y(1),Z(1),X(2),Y(2),Z(2),
     &              X(3),Y(3),Z(3),X(4),Y(4),Z(4),
     &              Ider(1,1),Ider(2,1),Ider(3,1),
     &              Ider(1,2),Ider(2,2),Ider(3,2),
     &              Ider(1,3),Ider(2,3),Ider(3,3),
     &              Ider(1,4),Ider(2,4),Ider(3,4),
     &              Alp,Cc,Icb,Ice,Spher,Screen,Iscr,
     &              Nbatch,Nfirst,Zscr)
               Gq(Kc) = 0.0D0
               If (Nbatch .Gt. 0) Gq(Kc) = Ddot(N1*N2*N3*N4,Gam,1,
     &                                         Zscr(Nfirst),1)
            Enddo
            Do Kc = 1, 3
               G(Kc,Iatd(Ia)) = G(Kc,Iatd(Ia)) + Fac*Gq(Kc)
               G(Kc,Iatd(Nd)) = G(Kc,Iatd(Nd)) - Fac*Gq(Kc)
            Enddo
         Enddo
  100 Continue
      Enddo
      Enddo
      Enddo
      Enddo

      Return
      End
