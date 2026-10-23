










      Subroutine Erd_wrvp(Lu,Icmpao,Ncrt,U,Ipsd,Ipu,Ipkb,Ipc,Ipsh,
     &                    Ipct,Itrm,Ctrm,Iky)
C
C The basis set block VMOL writes into IIII after the header record
C (vmol/readin.F) and vprops reads (props.F, incent.F, input.F):
C
C   POLYAINP : the centers (each symmetry unique atom, then its
C              images), then VMOL's pseudo-segmented Cartesian AOs,
C              one record per primitive with a nonzero coefficient
C   SYMTRANS : every SO as a sparse combination of those AOs, with
C              VMOL's unnormalized coefficients, then (spherical
C              basis) the contaminants
C   LABBASIS : center and component label of every SO
C   MULPOPDA : centers, charges and labels for Mulliken populations
C
C (CSYMTRAN is left out; only vmol2ja reads it.)
C
C Also the JOBARC records readin.F writes for vcceh (symadapt.f,
C drvsym.F): MULCINP (1 + the generators, as bits, that leave each
C unique atom in place), SYMTYPE (1 + the generators that flip x, y,
C z), NOSYMMOP, NOSYATOM, ATMSYMB0-9 (the unique atom names) and
C NBASCART.
C
C The AOs are numbered in VMOL's order: unique atom (MOL file order),
C shell (l, then group), Cartesian component, contraction, image.
C The SOs come from FULLSOAO (Erd_wrja), whose rows are the
C computational AOs; Ipsd maps those to VMOL's. The unique atoms and
C their names are those of the MOL file, matched to the MOL_ERDOED
C atoms by position; an image is labeled by the first operation r
C (VMOL's KB = r+1) that produces it. Charges are the effective ones
C with ECPs, as in VMOL.
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



c mxcbf.par : begin

c MXCBF := the maximum number of Cartesian basis functions (limited by vmol)

c This parameter is the same as MAXBASFN. Do NOT change this without changing
c maxbasfn.par as well.

      INTEGER MXCBF
      PARAMETER (MXCBF=1000)
c mxcbf.par : end

      Parameter (Ldchr = Mxcbf*4/10)
      Character*(Ldchr) Chrtmp
      Character*1 Chrnum
      Integer Mulnuc(100), Isytyp(3)
      Integer Icmpao(Ncrt), Ipsd(Ncrt), Ipu(Ncrt), Ipkb(Ncrt),
     &        Ipc(Ncrt), Ipsh(Ncrt), Ipct(Ncrt), Itrm(Ncrt), Iky(Ncrt)
      Dimension U(Ncrt,Ncrt), Ctrm(Ncrt)
      Character*4 Unam(Mxatms), Kwo(84)
      Character*80 Xline
      Integer Iua(Mxatms), Nimg(Mxatms), Iimr(8,Mxatms),
     &        Iatu(Mxatms), Ishl(8), Iord(36)
      Dimension Xmol(3,Mxatms), Dum(30)
      Logical Bexist
      Data Lumol /11/, Izero /0/, Ione /1/
      Data Kwo /'S   ', 'X   ', 'Y   ', 'Z   ',
     &          'XX  ', 'XY  ', 'XZ  ', 'YY  ',
     &          'YZ  ', 'ZZ  ', 'F300', 'F210',
     &          'F201', 'F120', 'F111', 'F102',
     &          'F030', 'F021', 'F012', 'F003',
     &          'G400', 'G310', 'G301', 'G220',
     &          'G211', 'G202', 'G130', 'G121',
     &          'G112', 'G103', 'G040', 'G031',
     &          'G022', 'G013', 'G004', 'H500',
     &          'H410', 'H401', 'H320', 'H311',
     &          'H302', 'H230', 'H221', 'H212',
     &          'H203', 'H140', 'H131', 'H122',
     &          'H113', 'H104', 'H050', 'H041',
     &          'H032', 'H023', 'H014', 'H005',
     &          'I600', 'I510', 'I501', 'I420', 'I411', 'I402',
     &          'I330', 'I321', 'I312', 'I303', 'I240', 'I231',
     &          'I222', 'I213', 'I204', 'I150', 'I141', 'I132',
     &          'I123', 'I114', 'I105', 'I060', 'I051', 'I042',
     &          'I033', 'I024', 'I015', 'I006'/

C Unique atoms (names and positions) from the MOL file.

      Inquire(File='MOL',Exist=Bexist)
      If (.Not. Bexist) Then
         Write(6,"(a)") " @Erd_wrvp: MOL is missing"
         Call Errex
      Endif
      Open(Unit=Lumol,File='MOL',Form='Formatted',Status='Old')
      Rewind(Lumol)
      Do I = 1, 5
         Read(Lumol,'(A)') Xline
      Enddo
      Nu = 0
 100  Continue
         Read(Lumol,'(A)') Xline
         If (Xline(1:6) .Eq. 'FINISH') Goto 200
         Read(Xline,"(F20.1,8I5)") Znuc, Natyp, Nl, (Ishl(I),I=1,Nl)
         Do Ia = 1, Natyp
            Nu = Nu + 1
            If (Nu .Gt. Mxatms) Then
               Write(6,"(a)") " @Erd_wrvp: too many atoms in MOL"
               Call Errex
            Endif
            Read(Lumol,"(A4,3F20.12)") Unam(Nu), (Xmol(K,Nu),K=1,3)
         Enddo
         Do L1 = 1, Nl
            Do Igrp = 1, Ishl(L1)
               Read(Lumol,"(2I5)") Np, Nc
               Do I = 1, Np
                  Read(Lumol,*) Ex, (Dum(J),J=1,Nc)
               Enddo
            Enddo
         Enddo
      Goto 100
 200  Continue
      Close(Lumol)

C MOL_ERDOED atom of each unique atom, and its images.

      Do Ia = 1, Natm
         Iatu(Ia) = 0
      Enddo
      Nnuc = 0
      Do Iu = 1, Nu
         Iua(Iu) = 0
         Do Ia = 1, Natm
            D = 0.0D0
            Do K = 1, 3
               D = D + Dabs(Xyzatm(K,Ia) - Xmol(K,Iu))
            Enddo
            If (D .Lt. 1.0D-6) Iua(Iu) = Ia
         Enddo
         If (Iua(Iu) .Eq. 0) Then
            Write(6,"(a,i4)") " @Erd_wrvp: no MOL_ERDOED atom at the "
     &                     // "position of MOL atom ", Iu
            Call Errex
         Endif
         Nimg(Iu) = 0
         Do Ir = 0, Nops-1
            Ib = Iatmap(Iua(Iu),Ir)
            If (Iatu(Ib) .Eq. 0) Then
               Iatu(Ib) = Iu
               Nimg(Iu) = Nimg(Iu) + 1
               Iimr(Nimg(Iu),Iu) = Ir
            Endif
         Enddo
         Nnuc = Nnuc + Nimg(Iu)
      Enddo

C vcceh's symmetry records. Generator g is operation 2**(g-1).

      Do K = 1, 3
         Isytyp(K) = 1
         Do Ig = 1, Ngen
            If (Iopsgn(K,2**(Ig-1)) .Lt. 0)
     &         Isytyp(K) = Isytyp(K) + 2**(Ig-1)
         Enddo
      Enddo
      If (Nu .Gt. 100) Then
         Write(6,"(a)") " @Erd_wrvp: MULCINP holds 100 atoms"
         Call Errex
      Endif
      Do Iu = 1, 100
         Mulnuc(Iu) = 0
      Enddo
      Do Iu = 1, Nu
         Mulnuc(Iu) = 1
         Do Ig = 1, Ngen
            Ir = 2**(Ig-1)
            Ifix = 1
            Do K = 1, 3
               If (Iopsgn(K,Ir) .Lt. 0 .And.
     &             Dabs(Xmol(K,Iu)) .Ge. 1.0D-6) Ifix = 0
            Enddo
            Mulnuc(Iu) = Mulnuc(Iu) + Ifix*Ir
         Enddo
      Enddo
      Call Putrec(20,'JOBARC','MULCINP ',100,Mulnuc)
      Call Putrec(20,'JOBARC','SYMTYPE ',3,Isytyp)
      Call Putrec(20,'JOBARC','NOSYMMOP',1,Ngen)
      Call Putrec(20,'JOBARC','NOSYATOM',1,Nu)
      Call Putrec(20,'JOBARC','NBASCART',1,Ncrt)
      Do J = 0, 9
         Do I = 1, Ldchr
            Chrtmp(I:I) = Char(0)
         Enddo
         Nper = Ldchr/4
         Do Iu = Nper*J+1, Min(Nu,Nper*(J+1))
            I = 4*(Iu-Nper*J-1)
            Chrtmp(I+1:I+4) = Unam(Iu)
         Enddo
         Write(Chrnum,'(I1)') J
         Call Putcrec(20,'JOBARC','ATMSYMB'//Chrnum,Ldchr,Chrtmp)
      Enddo

C VMOL's AOs, and the computational AO -> VMOL AO map.

      Do I = 1, Ncrt
         Ipsd(I) = 0
      Enddo
      Np = 0
      Do Iu = 1, Nu
         Ia = Iua(Iu)
         Do Is = 1, Nshl
            If (Ishatm(Is) .Eq. Ia) Then
               L  = Lshl(Is)
               Nk = (L+1)*(L+2)/2
               Do Ic = 1, Nk
                  Do It = 1, Nctsh(Is)
                     Do M = 1, Nimg(Iu)
                        Ir  = Iimr(M,Iu)
                        Np  = Np + 1
                        Isb = Ishmap(Is,Ir)
                        Iao = Iofcrt(Isb) + (It-1)*Nk + Ic
                        Ipsd(Icmpao(Iao)) = Np
                        Ipu(Np)  = Iu
                        Ipkb(Np) = Ir + 1
                        Ipc(Np)  = Ic
                        Ipsh(Np) = Is
                        Ipct(Np) = It
                     Enddo
                  Enddo
               Enddo
            Endif
         Enddo
      Enddo
      If (Np .Ne. Ncrt) Then
         Write(6,"(a,2i6)") " @Erd_wrvp: AO count mismatch ", Np, Ncrt
         Call Errex
      Endif
      Do I = 1, Ncrt
         If (Ipsd(I) .Eq. 0) Then
            Write(6,"(a,i6)") " @Erd_wrvp: AO not placed ", I
            Call Errex
         Endif
      Enddo

C POLYAINP. Primitives with a coefficient below 1.0D-8 are dropped,
C as VMOL does; Itrm holds the primitive count of every AO.

      Nprm = 0
      Do Ip = 1, Ncrt
         Is = Ipsh(Ip)
         Nps = Nprsh(Is)
         Itrm(Ip) = 0
         Do I = 1, Nps
            C = Ccraw(Iofcc(Is) + (Ipct(Ip)-1)*Nps + I)
            If (Dabs(C) .Gt. 1.0D-8) Itrm(Ip) = Itrm(Ip) + 1
         Enddo
         Nprm = Nprm + Itrm(Ip)
      Enddo

      Write(Lu) '********'
      Write(Lu) 'POLYAINP'
      Write(Lu) Nnuc, Izero, Izero, Izero
      Do Iu = 1, Nu
         Ia = Iua(Iu)
         Do M = 1, Nimg(Iu)
            Ir = Iimr(M,Iu)
            Write(Lu) Unam(Iu), Ir+1,
     &                (Xmol(K,Iu)*Dble(Iopsgn(K,Ir)),K=1,3), Zatm(Ia)
         Enddo
      Enddo
      Write(Lu) Nprm, Ncrt, Izero, Izero
      Write(Lu) (Itrm(Ip),Ip=1,Ncrt), Izero, Izero, Izero
      Do Ip = 1, Ncrt
         Is  = Ipsh(Ip)
         Nps = Nprsh(Is)
         L   = Lshl(Is)
         Nak = L*(L+1)*(L+2)/6
         Do I = 1, Nps
            C = Ccraw(Iofcc(Is) + (Ipct(Ip)-1)*Nps + I)
            If (Dabs(C) .Gt. 1.0D-8)
     &         Write(Lu) Unam(Ipu(Ip)), Ipkb(Ip), Kwo(Nak+Ipc(Ip)),
     &                   Alpha(Iofalp(Is)+I), C
         Enddo
      Enddo

C SYMTRANS: the Nso SOs, then the contaminants. VMOL lists the terms
C of a function image by image, components within an image.

      Write(Lu) '********'
      Write(Lu) 'SYMTRANS'
      Write(Lu) Nso, Izero, Izero, Izero
      Do Icol = 1, Ncrt
         Nt = 0
         Do I = 1, Ncrt
            If (Dabs(U(I,Icol)) .Ge. 1.0D-3) Then
               Ip = Ipsd(I)
               Key = Ipkb(Ip)*100 + Ipc(Ip)
               J = Nt
               Do While (J .Gt. 0)
                  If (Iky(J) .Le. Key) Goto 300
                  Iky(J+1)  = Iky(J)
                  Itrm(J+1) = Itrm(J)
                  Ctrm(J+1) = Ctrm(J)
                  J = J - 1
               Enddo
 300           Continue
               Iky(J+1)  = Key
               Itrm(J+1) = Ip
               Ctrm(J+1) = U(I,Icol)
               Nt = Nt + 1
            Endif
         Enddo
         Write(Lu) Nt, (Itrm(K),Ctrm(K),K=1,Nt), Izero
      Enddo

C LABBASIS and MULPOPDA. An SO is labeled with its unique atom and
C Cartesian component; a spherical one with the component that
C indexes its column of VMOL's table (Erd_vhord).

      Do Ip = 1, Nso
         Io  = Isoorb(Ip)
         Is0 = Iorbsh(1,Io)
         L   = Lshl(Is0)
         Nc  = Nfnsh(Is0)/Nctsh(Is0)
         Mc  = Mod(Isofn(Ip)-1,Nc) + 1
         If (Spher .And. L .Ge. 2) Then
            Call Erd_vhord(L,Iord)
            Mc = Iord(Mc)
         Endif
         Itrm(Ip) = Iatu(Ishatm(Is0))
         Iky(Ip)  = L*(L+1)*(L+2)/6 + Mc
      Enddo

      Write(Lu) '********'
      Write(Lu) 'LABBASIS'
      Do Ip = 1, Nso
         Write(Lu) Ip, Unam(Itrm(Ip)), Kwo(Iky(Ip)), Ione
      Enddo

      Write(Lu) '********'
      Write(Lu) 'MULPOPDA'
      Write(Lu) Nu, Nso, Izero, Izero
      Write(Lu) Nu, (Unam(Iu),Zatm(Iua(Iu)),Iu=1,Nu), Izero
      Do Ip = 1, Nso
         Write(Lu) Ip, Itrm(Ip), Iky(Ip), Izero
      Enddo

      Return
      End
