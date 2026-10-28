










      Subroutine Erg_one(G,Dt,W,N,Zscr,Lzscr,Iscr,Liscr)
C
C One-electron part of the gradient,
C
C    dE/dX += sum D(mu,nu) d[T+V](mu,nu)/dX - sum W(mu,nu) dS(mu,nu)/dX,
C
C over the symmetry unique shell pairs only: each pair is weighted
C with the size of its orbit and the result is symmetrized later
C (Erg_symz). Dt is the total density and W the energy weighted
C density, both over ERD functions. S and T depend on the two basis
C centers only (dX on the second = -dX on the first); V also depends
C on the nuclei, so its derivative is taken for every atom, with the
C basis centers on that atom and its nucleus differentiated together.
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
      Dimension G(3,Natm), Dt(N,N), W(N,N), Zscr(Lzscr)
      Integer Iscr(Liscr)
      Integer List(2), Icb(Mxpk), Ice(Mxpk), Id(3,3)
      Dimension Alp(Mxpk), Cc(Mxpk), Xn(Mxatms), Yn(Mxatms),
     &          Zn(Mxatms)
      Logical Screen, Lrep
      Data Screen /.False./
      Data Id /1,0,0, 0,1,0, 0,0,1/

      Do I = 1, Natm
         Xn(I) = Xyzatm(1,I)
         Yn(I) = Xyzatm(2,I)
         Zn(I) = Xyzatm(3,I)
      Enddo

      Do Is1 = 1, Nshl
      Do Is2 = 1, Is1
         List(1) = Is1
         List(2) = Is2
         Call Erg_orbit(List,2,Lrep,Nimg)
         If (.Not. Lrep) Goto 100
         Fac = Dble(Nimg)
         If (Is1 .Ne. Is2) Fac = 2.0D0*Fac
         Call Erd_pack(2,List,Alp,Nalp,Cc,Ncc,Icb,Ice,Ncsum,Mxpk)
         Ia1 = Ishatm(Is1)
         Ia2 = Ishatm(Is2)
         X1 = Xyzatm(1,Ia1)
         Y1 = Xyzatm(2,Ia1)
         Z1 = Xyzatm(3,Ia1)
         X2 = Xyzatm(1,Ia2)
         Y2 = Xyzatm(2,Ia2)
         Z2 = Xyzatm(3,Ia2)
         Nc1 = Nctsh(Is1)
         Nc2 = Nctsh(Is2)
         Np1 = Nprsh(Is1)
         Np2 = Nprsh(Is2)
         L1  = Lshl(Is1)
         L2  = Lshl(Is2)
         Nf1 = Nfnsh(Is1)
         Nf2 = Nfnsh(Is2)
         I0  = Ioffn(Is1)
         J0  = Ioffn(Is2)

C Overlap and kinetic energy: differentiate the first center.

         If (Ia1 .Ne. Ia2) Then
         Do K = 1, 3
         Do Iop = 1, 2
            If (Iop .Eq. 1) Then
               Call Oed__memory_ovl_derv_batch(Nalp,Ncc,Nc1,Nc2,
     &              Np1,Np2,L1,L2,X1,Y1,Z1,X2,Y2,Z2,
     &              Id(1,K),Id(2,K),Id(3,K),0,0,0,Alp,Cc,Spher,
     &              Imin,Iopt,Jmin,Jopt)
            Else
               Call Oed__memory_kin_derv_batch(Nalp,Ncc,Nc1,Nc2,
     &              Np1,Np2,L1,L2,X1,Y1,Z1,X2,Y2,Z2,
     &              Id(1,K),Id(2,K),Id(3,K),0,0,0,Alp,Cc,Spher,
     &              Imin,Iopt,Jmin,Jopt)
            Endif
            Call Erg_scrchk("Erg_one",Imin,Liscr,Jmin,Lzscr)
            Imax = Max(Imin,Min(Iopt,Liscr))
            Jmax = Max(Jmin,Min(Jopt,Lzscr))
            If (Iop .Eq. 1) Then
               Call Oed__gener_ovl_derv_batch(Imax,Jmax,Nalp,Ncc,
     &              Ncsum,Nc1,Nc2,Np1,Np2,L1,L2,X1,Y1,Z1,X2,Y2,Z2,
     &              Id(1,K),Id(2,K),Id(3,K),0,0,0,
     &              Alp,Cc,Icb,Ice,Spher,Screen,Iscr,
     &              Nbatch,Nfirst,Zscr)
            Else
               Call Oed__gener_kin_derv_batch(Imax,Jmax,Nalp,Ncc,
     &              Ncsum,Nc1,Nc2,Np1,Np2,L1,L2,X1,Y1,Z1,X2,Y2,Z2,
     &              Id(1,K),Id(2,K),Id(3,K),0,0,0,
     &              Alp,Cc,Icb,Ice,Spher,Screen,Iscr,
     &              Nbatch,Nfirst,Zscr)
            Endif
            If (Nbatch .Gt. 0) Then
               Sum = 0.0D0
               Do J = 1, Nf2
               Do I = 1, Nf1
                  X = Zscr(Nfirst + (J-1)*Nf1 + I - 1)
                  If (Iop .Eq. 1) Then
                     Sum = Sum - W(I0+I,J0+J)*X
                  Else
                     Sum = Sum + Dt(I0+I,J0+J)*X
                  Endif
               Enddo
               Enddo
               G(K,Ia1) = G(K,Ia1) + Fac*Sum
               G(K,Ia2) = G(K,Ia2) - Fac*Sum
            Endif
         Enddo
         Enddo
         Endif

C Nuclear attraction: every atom.

         Do Ic = 1, Natm
         Do K = 1, 3
            Kd1 = 0
            If (Ia1 .Eq. Ic) Kd1 = 1
            Kd2 = 0
            If (Ia2 .Eq. Ic) Kd2 = 1
            Call Oed__memory_nai_derv_batch(Nalp,Ncc,Nc1,Nc2,
     &           Np1,Np2,L1,L2,X1,Y1,Z1,X2,Y2,Z2,
     &           Natm,Xn,Yn,Zn,Ic,
     &           Kd1*Id(1,K),Kd1*Id(2,K),Kd1*Id(3,K),
     &           Kd2*Id(1,K),Kd2*Id(2,K),Kd2*Id(3,K),
     &           Id(1,K),Id(2,K),Id(3,K),Alp,Cc,Spher,
     &           Imin,Iopt,Jmin,Jopt)
            Call Erg_scrchk("Erg_one",Imin,Liscr,Jmin,Lzscr)
            Imax = Max(Imin,Min(Iopt,Liscr))
            Jmax = Max(Jmin,Min(Jopt,Lzscr))
            Call Oed__gener_nai_derv_batch(Imax,Jmax,Nalp,Ncc,
     &           Ncsum,Nc1,Nc2,Np1,Np2,L1,L2,X1,Y1,Z1,X2,Y2,Z2,
     &           Natm,Xn,Yn,Zn,Zatm,Ic,
     &           Kd1*Id(1,K),Kd1*Id(2,K),Kd1*Id(3,K),
     &           Kd2*Id(1,K),Kd2*Id(2,K),Kd2*Id(3,K),
     &           Id(1,K),Id(2,K),Id(3,K),
     &           Alp,Cc,Icb,Ice,Spher,Screen,Iscr,
     &           Nbatch,Nfirst,Zscr)
            If (Nbatch .Gt. 0) Then
               Sum = 0.0D0
               Do J = 1, Nf2
               Do I = 1, Nf1
                  Sum = Sum + Dt(I0+I,J0+J)*
     &                  Zscr(Nfirst + (J-1)*Nf1 + I - 1)
               Enddo
               Enddo
               G(K,Ic) = G(K,Ic) + Fac*Sum
            Endif
         Enddo
         Enddo
  100    Continue
      Enddo
      Enddo

      Return
      End
