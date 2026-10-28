










      Subroutine Erd_cmpord(Icmpao,Ncrt,Iz2a,Natoms,Coord,Scr,Imap,
     &                      Ipop,Imem,Iaocen,Iaoang,Iaocn0,Iaoan0,
     &                      Nbsorb,Naoorb)
C
C Computational (VMOL) order of the Cartesian AOs, and the JOBARC
C records describing it. The rest of ACES II assumes this order
C (vscf/setup.F, for instance, rebuilds it from NBASATOM): orbit by
C orbit (joda's COMPPOPV/COMPMEMB), centers of an orbit in REMAP order
C (MAP2ZMAT, as vmol2ja writes it), and within an orbit
C
C    comp = (orbit offset) + (v-1)*npop + member
C
C where v runs over the functions of one atom in VMOL order: angular
C momentum, then Cartesian component, then all contractions of that
C angular momentum. The ZMAT AO order is the ERD Cartesian order
C (atoms in ZMAT order, shells as in MOL_ERDOED, contraction major).
C
C Output: Icmpao(i) is the computational index of ZMAT AO i.
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




      Integer Icmpao(Ncrt), Iz2a(Natoms), Imap(Natoms), Ipop(*),
     &        Imem(Natoms), Iaocen(Ncrt), Iaoang(Ncrt),
     &        Iaocn0(Ncrt), Iaoan0(Ncrt), Nbsorb(*), Naoorb(*)
      Dimension Coord(3,Natoms), Scr(*)
      Integer Ia2z(Mxatms), Ifsh(Mxatms+1)
      Character*4 Cmppgp

C ZMAT centers (dummies included) <-> MOL_ERDOED atoms (real only).

      Call Getrec(20,'JOBARC','ATOMCHRG',Natoms,Imap)
      Nreal = 0
      Do Iz = 1, Natoms
         Iz2a(Iz) = 0
         If (Imap(Iz) .Ne. 0) Then
            Nreal = Nreal + 1
            Iz2a(Iz) = Nreal
            If (Nreal .Le. Mxatms) Ia2z(Nreal) = Iz
         Endif
      Enddo
      If (Nreal .Ne. Natm) Then
         Write(6,"(a,2i5)") " @Erd_cmpord: atom count mismatch ",
     &                      Nreal, Natm
         Call Errex
      Endif
      Call Getrec(20,'JOBARC','COORD   ',3*Natoms*Iintfp,Coord)
      Do Ia = 1, Natm
         Dist = 0.0D0
         Do K = 1, 3
            Dist = Dist + Dabs(Coord(K,Ia2z(Ia)) - Xyzatm(K,Ia))
         Enddo
         If (Dist .Gt. 1.0D-6) Then
            Write(6,"(2a,i4)") " @Erd_cmpord: COORD and MOL_ERDOED ",
     &                        "disagree for atom ", Ia
            Call Errex
         Endif
      Enddo

C Computational centers, exactly as vmol2ja's driver.F does it.

      Call Getcrec(-1,'JOBARC','COMPPTGP',4,Cmppgp)
      Call Getrec(-1,'JOBARC','COMPORDR',1,Iordgp)
      Call Getrec(-1,'JOBARC','COMPNORB',1,Norbit)
      Call Getrec(-1,'JOBARC','COMPPOPV',Norbit,Ipop)
      Call Getrec(-1,'JOBARC','COMPMEMB',Natoms,Imem)
      Call Remap(Cmppgp,Iordgp,Natoms,Norbit,Coord,Scr,Imem,Ipop,
     &           Imap)
      Call Putrec(20,'JOBARC','MAP2ZMAT',Natoms,Imap)

C First shell of every atom.

      Do Ia = 1, Natm+1
         Ifsh(Ia) = Nshl + 1
      Enddo
      Do Is = Nshl, 1, -1
         Ifsh(Ishatm(Is)) = Is
      Enddo
      Ifsh(Natm+1) = Nshl + 1
      Do Ia = Natm, 1, -1
         If (Ifsh(Ia) .Gt. Ifsh(Ia+1)) Ifsh(Ia) = Ifsh(Ia+1)
      Enddo

      Do I = 1, Ncrt
         Icmpao(I) = 0
      Enddo
      Icomp = 0
      Icen  = 0
      Do Io = 1, Norbit
         Npop  = Ipop(Io)
         Norbm0 = Icomp
         Nsoo  = 0
         Do M = 1, Npop
            Icen = Icen + 1
            Iz   = Imap(Icen)
            If (Iz .Le. 0) Then
               Write(6,"(a,i4)") " @Erd_cmpord: REMAP failed, center",
     &                           Icen
               Call Errex
            Endif
            Ia = Iz2a(Iz)
            If (Ia .Eq. 0) Then
               Write(6,"(a)") " @Erd_cmpord: dummy atom in an orbit"
               Call Errex
            Endif

C VMOL order on this atom: l, then component, then contraction.

            Iv = 0
            Lmx = -1
            Do Is = Ifsh(Ia), Ifsh(Ia+1)-1
               Lmx = Max(Lmx,Lshl(Is))
            Enddo
            Do L = 0, Lmx
               Ncmp = (L+1)*(L+2)/2
               Do Ic = 1, Ncmp
                  Do Is = Ifsh(Ia), Ifsh(Ia+1)-1
                     If (Lshl(Is) .Eq. L) Then
                        Do It = 1, Nctsh(Is)
                           Iv = Iv + 1
                           Iao = Iofcrt(Is) + (It-1)*Ncmp + Ic
                           Ic0 = Norbm0 + (Iv-1)*Npop + M
                           Icmpao(Iao) = Ic0
                           Iaocen(Ic0) = Icen
                           Iaoang(Ic0) = L
                           Iaocn0(Iao) = Iz
                           Iaoan0(Iao) = L
                        Enddo
                     Endif
                  Enddo
               Enddo
            Enddo
            Do Is = Ifsh(Ia), Ifsh(Ia+1)-1
               Nsoo = Nsoo + Nfnsh(Is)
            Enddo
            If (M .Eq. 1) Nv1 = Iv
            If (Iv .Ne. Nv1) Then
               Write(6,"(a)") " @Erd_cmpord: orbit members differ"
               Call Errex
            Endif
         Enddo
         Icomp = Icomp + Npop*Nv1
         Naoorb(Io) = Npop*Nv1
         Nbsorb(Io) = Nsoo
      Enddo

      If (Icomp .Ne. Ncrt) Then
         Write(6,"(a,2i6)") " @Erd_cmpord: AO count mismatch ", Icomp,
     &                      Ncrt
         Call Errex
      Endif
      Do I = 1, Ncrt
         If (Icmpao(I) .Eq. 0) Then
            Write(6,"(a,i6)") " @Erd_cmpord: AO not placed ", I
            Call Errex
         Endif
      Enddo

      Call Putrec(20,'JOBARC','CENTERBF',Ncrt,Iaocen)
      Call Putrec(20,'JOBARC','ANGMOMBF',Ncrt,Iaoang)
      Call Putrec(20,'JOBARC','CNTERBF0',Ncrt,Iaocn0)
      Call Putrec(20,'JOBARC','ANMOMBF0',Ncrt,Iaoan0)
      Call Putrec(20,'JOBARC','NBASATOM',Norbit,Nbsorb)
      Call Putrec(20,'JOBARC','NAOBFORB',Norbit,Naoorb)

      Return
      End
