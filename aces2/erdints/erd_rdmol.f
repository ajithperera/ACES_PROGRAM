










      Subroutine Erd_rdmol(Kgen,Tla)
C
C Read the MOL_ERDOED file written by joda (MKERDOED). Unlike the MOL
C file it lists every real center (in ZMAT order), not just the
C symmetry unique ones. The header line also carries the VMOL
C generators (Kgen), which fix the irrep order, and the integral
C neglect threshold (Tla). Fills the erdbas common block.
C
C The contraction coefficients are renormalized the way ERD/OED
C expect them (normalized contracted function over normalized
C primitives), exactly as a3ecp/molread.F READ_BASIS_INFO does.
C
C Ajith Perera, 10/2026.
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

      Character*3 Kgen(3)
      Character*80 Xline
      Character*4 Atmnam
      Integer Ishl(8)
      Double Precision Salpha(200), Scoef(200,30), Sraw(200,30)
      Logical Bexist
      Data Luin /10/

      Inquire(File='MOL_ERDOED',Exist=Bexist)
      If (.Not. Bexist) Then
         Write(6,"(a)") " @Erd_rdmol: MOL_ERDOED is missing"
         Call Errex
      Endif
      Open(Unit=Luin,File='MOL_ERDOED',Form='Formatted',Status='Old')
      Rewind(Luin)

      Read(Luin,"(10x,9i5)") (Ishl(I),I=1,8),Icart
      Spher = (Icart .Eq. 0)
      Read(Luin,'(A)') Xline
      Read(Luin,'(A)') Xline
      Read(Luin,"(2I5,3A3,1X,E10.2)") Ntyp, Ngen, (Kgen(I),I=1,3),
     &                               Tla
      Read(Luin,'(A)') Xline

      Natm   = 0
      Nshl   = 0
      Nalp   = 0
      Ncc    = 0
      Ncb    = 0
      Nfnerd = 0
      Ncart  = 0

C The header count of atom types includes dummy atoms, which have no
C block, so read blocks until FINISH.

 100  Continue
         Read(Luin,'(A)') Xline
         If (Xline(1:6) .Eq. 'FINISH') Goto 200
         Read(Xline,"(F20.1,8I5)") Znuc, Natyp, Nl, (Ishl(I),I=1,Nl)
         If (Natyp .Ne. 1) Then
            Write(6,"(a)") " @Erd_rdmol: one atom per block expected"
            Call Errex
         Endif
         Iatm0 = Natm
         Do Ia = 1, Natyp
            Natm = Natm + 1
            If (Natm .Gt. Mxatms) Then
               Write(6,"(a,i5)") " @Erd_rdmol: too many atoms, max ",
     &                           Mxatms
               Call Errex
            Endif
            Read(Luin,"(A4,3F20.12)") Atmnam, (Xyzatm(K,Natm),K=1,3)
            Zatm(Natm) = Znuc
         Enddo

C The basis set follows the atom lines and is common to all atoms
C of this type; read it once and replicate it on each atom.

         Do L1 = 1, Nl
            Do Igrp = 1, Ishl(L1)
               Read(Luin,"(2I5)") Np, Nc
               If (Np .Gt. 200 .Or. Nc .Gt. 30) Then
                  Write(6,"(a,2i5)") " @Erd_rdmol: shell too large ",
     &                               Np, Nc
                  Call Errex
               Endif
               Do I = 1, Np
                  Read(Luin,*) Salpha(I), (Scoef(I,J),J=1,Nc)
                  Do J = 1, Nc
                     Sraw(I,J) = Scoef(I,J)
                  Enddo
               Enddo
               Call Erd_normcc(Np,Nc,L1-1,Salpha,Scoef,200)

               Do Ia = 1, Natyp
                  Nshl = Nshl + 1
                  If (Nshl .Gt. Mxshl) Then
                     Write(6,"(a,i6)") " @Erd_rdmol: too many shells,",
     &                                 " max ", Mxshl
                     Call Errex
                  Endif
                  Ishatm(Nshl) = Iatm0 + Ia
                  Lshl(Nshl)   = L1 - 1
                  Nprsh(Nshl)  = Np
                  Nctsh(Nshl)  = Nc
                  Iofalp(Nshl) = Nalp
                  Iofcc(Nshl)  = Ncc
                  Iofcb(Nshl)  = Ncb
                  If (Nalp+Np .Gt. Mxprmt .Or. Ncc+Np*Nc .Gt. Mxcct
     &                .Or. Ncb+Nc .Gt. Mxctrt) Then
                     Write(6,"(a)") " @Erd_rdmol: basis too large"
                     Call Errex
                  Endif
                  Do I = 1, Np
                     Alpha(Nalp+I) = Salpha(I)
                  Enddo
                  Do J = 1, Nc
                     Ibeg = 0
                     Iend = 0
                     Do I = 1, Np
                        Ccoef(Ncc+(J-1)*Np+I) = Scoef(I,J)
                        Ccraw(Ncc+(J-1)*Np+I) = Sraw(I,J)
                        If (Scoef(I,J) .Ne. 0.0D0) Then
                           If (Ibeg .Eq. 0) Ibeg = I
                           Iend = I
                        Endif
                     Enddo
                     Iccbeg(Ncb+J) = Ibeg
                     Iccend(Ncb+J) = Iend
                  Enddo
                  Nalp = Nalp + Np
                  Ncc  = Ncc  + Np*Nc
                  Ncb  = Ncb  + Nc
                  L = L1 - 1
                  Ncmpc = (L+1)*(L+2)/2
                  If (Spher) Then
                     Ncmp = 2*L + 1
                  Else
                     Ncmp = Ncmpc
                  Endif
                  Nfnsh(Nshl) = Ncmp*Nc
               Enddo
            Enddo
         Enddo
      Goto 100
 200  Continue
      Close(Luin)

      Ncart = 0
      Do Is = 1, Nshl
         Ioffn(Is) = Nfnerd
         Nfnerd = Nfnerd + Nfnsh(Is)
         Iofcrt(Is) = Ncart
         L = Lshl(Is)
         Ncart = Ncart + Nctsh(Is)*(L+1)*(L+2)/2
      Enddo

      Return
      End
