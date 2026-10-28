










      Subroutine Erd_symgrp(Kgen)
C
C Set up the computational point group from the VMOL generators in
C MOL_ERDOED. ACES II works in D2h or one of its subgroups, so every
C operation is diag(+-1,+-1,+-1). A generator string lists the
C coordinates it flips ('X  ' is the yz mirror, 'XY ' is C2(z), 'XYZ'
C is the inversion). Operation r (0 <= r < 2**Ngen) is the product of
C the generators whose bit is set in r, the same numbering VMOL uses
C (readin.F, KA array), so irrep Gamma (0 based) has character
C (-1)**popcount(iand(Gamma,r)) and is ACES irrep Gamma+1.
C
C Also maps every atom and shell onto its image under each operation.
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

      Character*3 Kgen(3)
      Character*1 Kxyz(3)
      Integer Igsgn(3,3), Ifirst(Mxatms+1)
      Data Kxyz /'X','Y','Z'/
      Data Tol /1.0D-6/

      Do Ig = 1, 3
         Do K = 1, 3
            Igsgn(K,Ig) = 1
         Enddo
      Enddo
      Do Ig = 1, Ngen
         Do Ic = 1, 3
            Do K = 1, 3
               If (Kgen(Ig)(Ic:Ic) .Eq. Kxyz(K)) Igsgn(K,Ig) = -1
            Enddo
         Enddo
      Enddo

      Nops = 2**Ngen
      Do Ir = 0, Nops-1
         Do K = 1, 3
            Iopsgn(K,Ir) = 1
            Do Ig = 1, Ngen
               If (Btest(Ir,Ig-1)) Iopsgn(K,Ir) = Iopsgn(K,Ir)*
     &                                            Igsgn(K,Ig)
            Enddo
         Enddo
      Enddo

C Atom images.

      Do Ir = 0, Nops-1
         Do Ia = 1, Natm
            Jmatch = 0
            Do Ja = 1, Natm
               Dist = 0.0D0
               Do K = 1, 3
                  Dist = Dist + Dabs(Xyzatm(K,Ja) -
     &                          Dble(Iopsgn(K,Ir))*Xyzatm(K,Ia))
               Enddo
               If (Dist .Lt. Tol .And. Zatm(Ja) .Eq. Zatm(Ia))
     &            Jmatch = Ja
            Enddo
            If (Jmatch .Eq. 0) Then
               Write(6,"(a,i4,a,i2)") " @Erd_symgrp: atom ", Ia,
     &               " has no image under operation ", Ir
               Call Errex
            Endif
            Iatmap(Ia,Ir) = Jmatch
         Enddo
      Enddo

C Shell images. The shells of an atom are contiguous and equivalent
C atoms carry the same basis, so the image of the k-th shell of atom
C a is the k-th shell of the image of a. Check that it really is.

      Do Ia = 1, Natm+1
         Ifirst(Ia) = 0
      Enddo
      Do Is = Nshl, 1, -1
         Ifirst(Ishatm(Is)) = Is
      Enddo
      Do Ir = 0, Nops-1
         Do Is = 1, Nshl
            Ia = Ishatm(Is)
            Ja = Iatmap(Ia,Ir)
            Js = Ifirst(Ja) + (Is - Ifirst(Ia))
            Ibad = 0
            If (Js .Gt. Nshl) Then
               Ibad = 1
            Else If (Ishatm(Js) .Ne. Ja .Or. Lshl(Js) .Ne. Lshl(Is)
     &          .Or. Nprsh(Js) .Ne. Nprsh(Is) .Or. Nctsh(Js) .Ne.
     &          Nctsh(Is)) Then
               Ibad = 1
            Else
               Do Ip = 1, Nprsh(Is)
                  If (Alpha(Iofalp(Js)+Ip) .Ne. Alpha(Iofalp(Is)+Ip))
     &               Ibad = 1
               Enddo
            Endif
            If (Ibad .Ne. 0) Then
               Write(6,"(a,i5,a)") " @Erd_symgrp: shell ", Is,
     &               " has no matching image; equivalent atoms"
               Write(6,"(a)") "  must carry the same basis set."
               Call Errex
            Endif
            Ishmap(Is,Ir) = Js
         Enddo
      Enddo

      Return
      End
