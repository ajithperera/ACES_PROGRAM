










      Subroutine Erd_soord(Icmpao,Ncrt,Key,Inew,Iscr,Nsc)
C
C Put the SOs of each irrep in VMOL's order and phase: VMOL sorts them
C by the computational index of their function on the first member
C (orbit, then l, component, contraction, Erd_cmpord) and gives that
C member the + sign. Spherical functions (L >= 2) are ordered the same
C way with the harmonic (VMOL's order, Erd_vsph) in place of the
C Cartesian component. This reproduces VMOL's SYMTRANS (AO2SO).
C
C Rewrites Isoirr, Isoorb, Isofn, Isosgn and Lorbso (erdso common).
C Key, Inew are Nso long, Iscr is Nsc = 11*Nso long scratch.
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

      Integer Icmpao(Ncrt), Key(Nso), Inew(Nso), Iscr(Nsc)
      Logical Lsph

C Sort key of every SO and the member that carries it.

      Do Ip = 1, Nso
         Io  = Isoorb(Ip)
         Is0 = Iorbsh(1,Io)
         L   = Lshl(Is0)
         Nk  = (L+1)*(L+2)/2
         Nc  = Nfnsh(Is0)/Nctsh(Is0)
         It  = (Isofn(Ip)-1)/Nc + 1
         Mc  = Mod(Isofn(Ip)-1,Nc) + 1
         Kmin = Ncrt + 1
         Mmin = 0
         Lsph = Spher .And. L .Ge. 2
         Do M = 1, Norbm(Io)
            If (Lsph) Then
               Iao = Iofcrt(Iorbsh(M,Io)) + (It-1)*Nk + 1
            Else
               Iao = Iofcrt(Iorbsh(M,Io)) + (It-1)*Nk + Mc
            Endif
            If (Icmpao(Iao) .Lt. Kmin) Then
               Kmin = Icmpao(Iao)
               Mmin = M
            Endif
         Enddo
         If (.Not. Lsph) Then
            Key(Ip) = Kmin*10000
         Else

C VMOL orders spherical functions like Cartesian ones: harmonic, then
C all contractions of that L on the atom. Key: start of the atom's L
C block, harmonic, contraction counted over the atom's L shells.

            Is  = Iorbsh(Mmin,Io)
            Ia  = Ishatm(Is)
            Nt  = 0
            Ito = 0
            Ifs = 0
            Do Js = 1, Nshl
               If (Ishatm(Js) .Eq. Ia .And. Lshl(Js) .Eq. L) Then
                  If (Ifs .Eq. 0) Ifs = Js
                  If (Js .Lt. Is) Ito = Ito + Nctsh(Js)
                  Nt = Nt + Nctsh(Js)
               Endif
            Enddo
            Kbl = Ncrt + 1
            Do M = 1, Norbm(Io)
               Ks = Ishmap(Ifs,Iorbop(M,Io))
               Kbl = Min(Kbl,Icmpao(Iofcrt(Ks)+1))
            Enddo
            Key(Ip) = Kbl*10000 + (Mc-1)*Nt + Ito + It
         Endif
         If (Isosgn(Mmin,Ip) .Lt. 0) Then
            Do M = 1, Norbm(Io)
               Isosgn(M,Ip) = -Isosgn(M,Ip)
            Enddo
         Endif
      Enddo

C New position of every SO: rank of its key within its irrep (the
C irreps keep their blocks).

      Iso0 = 0
      Do Ig = 1, Nirr
         Do Ip = Iso0+1, Iso0+Nbfirr(Ig)
            Irank = 0
            Do Iq = Iso0+1, Iso0+Nbfirr(Ig)
               If (Key(Iq) .Lt. Key(Ip)) Irank = Irank + 1
            Enddo
            Inew(Ip) = Iso0 + Irank + 1
         Enddo
         Iso0 = Iso0 + Nbfirr(Ig)
      Enddo

C Permute the SO tables (Iscr holds the old ones).

      Do Ip = 1, Nso
         Iscr(Ip)       = Isoirr(Ip)
         Iscr(Nso+Ip)   = Isoorb(Ip)
         Iscr(2*Nso+Ip) = Isofn(Ip)
         Do M = 1, 8
            Iscr(3*Nso+(Ip-1)*8+M) = Isosgn(M,Ip)
         Enddo
      Enddo
      Do Ip = 1, Nso
         Jp = Inew(Ip)
         Isoirr(Jp) = Iscr(Ip)
         Isoorb(Jp) = Iscr(Nso+Ip)
         Isofn(Jp)  = Iscr(2*Nso+Ip)
         Do M = 1, 8
            Isosgn(M,Jp) = Iscr(3*Nso+(Ip-1)*8+M)
         Enddo
      Enddo
      Do I = 1, Nso
         Lorbso(I) = Inew(Lorbso(I))
      Enddo

      Return
      End
