#include "flags.h"
      Subroutine Erd_driver(Icore,Maxcor)
C
C Top level of xerdints. Icore(Maxcor) is the integer work space
C (a double takes Iintfp integers).
C
C   1. MOL_ERDOED -> shells, atoms, VMOL generators     (Erd_rdmol)
C      (ECPs: ECP1INTS and effective charges after 4)   (Erd_ecp)
C   2. operations, atom and shell images                (Erd_symgrp)
C   3. SO basis from the shell orbits                   (Erd_sobld)
C   4. VMOL computational AO order + AO label records   (Erd_cmpord)
C      and VMOL's SO order and phase                     (Erd_soord)
C   5. SO <-> AO transformation records                 (Erd_wrja)
C   6. OED one-electron integrals -> SO -> IIII         (Erd_ao1e,
C                                                        Erd_so1e)
C   7. ERD two-electron integrals -> SO -> IIII, IIJJ, IJIJ, IJKL
C                                                       (Erd_so2e)
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"
#include "machsp.com"
#include "flags.com"

      Integer Icore(Maxcor), Lu(4), Nwr(4)
      Character*3 Kgen(3)
      Character*8 Title(24), Fname8(4)
      Character*80 Fname
      Logical Debug, Lfull, Lecp
      Data Fname8 /'IIII    ','IIJJ    ','IJIJ    ','IJKL    '/

      Lecp = (Iflags(h_IFLAGS_ecp) .Ne. 0)

C Development aid: a file ERDINTS.DEBUG in the work directory asks
C for dumps of intermediate matrices and for the full (unreduced)
C two-electron symmetry check.

      Inquire(File='ERDINTS.DEBUG',Exist=Debug)
      Lfull = Debug

      Call Erd_rdmol(Kgen,Tla)
      Call Erd_symgrp(Kgen)
      Call Erd_sobld
      Write(6,"(a,i4,a,i5,a,i6,a,l2)") " @Erd_driver: atoms ", Natm,
     &     " shells ", Nshl, " functions ", Nfnerd, " spherical", Spher
      Write(6,"(a,i2,3(1x,a3))") " @Erd_driver: generators ", Ngen,
     &     (Kgen(I),I=1,Ngen)
      Write(6,"(a,i5,a,8i5)") " @Erd_driver: shell orbits ", Norb,
     &     "  SOs per irrep ", (Nbfirr(I),I=1,Nirr)

C Computational AO order and the AO <-> SO records. Only Icmpao
C (at J000) is kept beyond this block.

      Call Getrec(20,'JOBARC','NATOMS  ',1,Natoms)
      Ncrt = Ncart
      J000 = 1
      J010 = J000 + Ncrt
      J020 = J010 + Natoms
      J030 = J020 + 3*Natoms*Iintfp
      J040 = J030 + 3*8*Iintfp
      J050 = J040 + Natoms
      J060 = J050 + Natoms
      J070 = J060 + Natoms
      J080 = J070 + Ncrt
      J090 = J080 + Ncrt
      J100 = J090 + Ncrt
      J110 = J100 + Ncrt
      J120 = J110 + Natoms
      J130 = J120 + Natoms
      If (J130 .Gt. Maxcor) Call Insmem("Erd_driver",J130,Maxcor)
      Call Erd_cmpord(Icore(J000),Ncrt,Icore(J010),Natoms,
     &                Icore(J020),Icore(J030),Icore(J040),Icore(J050),
     &                Icore(J060),Icore(J070),Icore(J080),Icore(J090),
     &                Icore(J100),Icore(J110),Icore(J120))

C ECP integrals (ECP1INTS); from here on Zatm holds the effective
C nuclear charges. Iz2a (J010) and Coord (J020) come from Erd_cmpord.

      If (Lecp) Call Erd_ecp(Natoms,Icore(J010),Icore(J020),
     &                       Icore(J040))

C Nuclear repulsion of the MOL_ERDOED atoms (effective charges with
C ECPs). Like vmol2ja, overwrite joda's NUCREP: at finite-difference
C points and in later optimization cycles joda's is not current.

      Repnuc = 0.0D0
      Do I = 2, Natm
         Do J = 1, I-1
            R = Dsqrt((Xyzatm(1,I)-Xyzatm(1,J))**2 +
     &                (Xyzatm(2,I)-Xyzatm(2,J))**2 +
     &                (Xyzatm(3,I)-Xyzatm(3,J))**2)
            Repnuc = Repnuc + Zatm(I)*Zatm(J)/R
         Enddo
      Enddo
      Call Putrec(20,'JOBARC','NUCREP  ',Iintfp,Repnuc)

      If (J010+13*Nso .Gt. Maxcor) Call Insmem("Erd_driver",
     &                                         J010+13*Nso,Maxcor)
      Call Erd_soord(Icore(J000),Ncrt,Icore(J010),Icore(J010+Nso),
     &               Icore(J010+2*Nso),11*Nso)

      I000 = J010
      I010 = I000 + Ncrt*Ncrt*Iintfp
      I020 = I010 + Ncrt*Ncrt*Iintfp
      I030 = I020 + Ncrt*Ncrt*Iintfp
      I040 = I030 + Ncrt*Ncrt*Iintfp
      I050 = I040 + Ncrt
      If (I050 .Gt. Maxcor) Call Insmem("Erd_driver",I050,Maxcor)
      Call Erd_wrja(Icore(J000),Ncrt,Icore(I000),Icore(I010),
     &              Icore(I020),Icore(I030),Icore(I040),Ncrt*Ncrt)
      If (Spher) Call Erd_sphchk(Icore(I000),Maxcor-I000+1)

C Integral files, opened as VMOL does (IIJJ, IJIJ only with
C symmetry, IJKL only with more than two irreps).

      Lu(1) = 10
      Lu(2) = 81
      Lu(3) = 82
      Lu(4) = 83
      Do I = 1, 4
         If (I .Eq. 1 .Or. (I .Le. 3 .And. Ngen .Ge. 1) .Or.
     &       Ngen .Ge. 2) Then
            Call Gfname(Fname8(I),Fname,Ilength)
            Open(Unit=Lu(I),File=Fname(1:Ilength),Form='Unformatted',
     &           Access='Sequential',Status='Replace')
         Else
            Lu(I) = 0
         Endif
      Enddo

      Do I = 1, 24
         Title(I) = ' '
      Enddo
      Title(1) = 'xerdints'
      Title(2) = ' OED/ERD'
      Idosph = 0
      If (Spher) Idosph = 1
      Write(Lu(1)) Title, Nirr, (Nbfirr(I),I=1,Nirr), Repnuc, Idosph

C The basis set block vprops reads (VMOL's, see Erd_wrvp).

      I010 = J010 + Ncrt*Ncrt*Iintfp
      I020 = I010 + Ncrt
      I030 = I020 + Ncrt
      I040 = I030 + Ncrt
      I050 = I040 + Ncrt
      I060 = I050 + Ncrt
      I070 = I060 + Ncrt
      I080 = I070 + Ncrt
      I090 = I080 + Ncrt*Iintfp
      I100 = I090 + Ncrt
      If (I100 .Gt. Maxcor) Call Insmem("Erd_driver",I100,Maxcor)
      Call Getrec(20,'JOBARC','FULLSOAO',Ncrt*Ncrt*Iintfp,Icore(J010))
      Call Erd_wrvp(Lu(1),Icore(J000),Ncrt,Icore(J010),Icore(I010),
     &              Icore(I020),Icore(I030),Icore(I040),Icore(I050),
     &              Icore(I060),Icore(I070),Icore(I080),Icore(I090))

C One-electron integrals.

      N   = Nfnerd
      N2  = N*N*Iintfp
      I0s = J010
      I0t = I0s + N2
      I0v = I0t + N2
      I0b = I0v + N2
      I0z = I0b + Nso*Nso*Iintfp
      Left = Maxcor - I0z + 1
      If (Left .Le. 0) Call Insmem("Erd_driver",I0z,Maxcor)
      Lzscr = (Left/2)/Iintfp
      Liscr = Left/2
      I0i = I0z + Lzscr*Iintfp

      Call Erd_ao1e(Icore(I0s),Icore(I0t),Icore(I0v),N,
     &              Icore(I0z),Lzscr,Icore(I0i),Liscr)
      If (Debug) Call Erd_dump1e(Icore(I0s),Icore(I0t),Icore(I0v),N)

      Call Erd_so1e(Icore(I0s),N,Icore(I0b),Asym1)
      Call Erd_wr1e(Lu(1),'OVERLAP ',Icore(I0b),Icore(I0z),Icore(I0i))
      Call Daxpy(N*N,1.0D0,Icore(I0t),1,Icore(I0v),1)
      Call Erd_so1e(Icore(I0v),N,Icore(I0b),Asym2)
      Call Erd_wr1e(Lu(1),'ONEHAMIL',Icore(I0b),Icore(I0z),Icore(I0i))
      Call Erd_so1e(Icore(I0t),N,Icore(I0b),Asym3)
      Call Erd_wr1e(Lu(1),'KINETINT',Icore(I0b),Icore(I0z),Icore(I0i))
      Asym = Max(Asym1,Asym2,Asym3)
      Write(6,"(2a,1p,e10.2)") " @Erd_driver: largest one-electron SO",
     &     " integral between different irreps", Asym
      If (Asym .Gt. 1.0D-10) Then
         Write(6,"(a)") " @Erd_driver: the SO basis is not symmetry",
     &                  " adapted"
         Call Errex
      Endif

C Two-electron integrals: two transformation buffers and the ERD
C scratch share what is left.

      I0w  = J010
      Left = Maxcor - I0w + 1
      Lw   = (Left/4)/Iintfp
      I0w2 = I0w  + Lw*Iintfp
      I0z  = I0w2 + Lw*Iintfp
      Lzscr = (Left/4)/Iintfp
      I0i  = I0z + Lzscr*Iintfp
      Liscr = Maxcor - I0i + 1

      Do I = 1, 4
         If (Lu(I) .Gt. 0) Then
            Write(Lu(I)) '********'
            Write(Lu(I)) 'TWOELSUP'
         Endif
      Enddo
      Call Erd_so2e(Lu,Tla,Lfull,Icore(I0w),Icore(I0w2),Lw,
     &              Icore(I0z),Lzscr,Icore(I0i),Liscr,Asym,Nwr)
      Write(6,"(a,4i12)") " @Erd_driver: IIII IIJJ IJIJ IJKL ",
     &     (Nwr(I),I=1,4)
      If (Lfull) Then
         Write(6,"(2a,1p,e10.2)") " @Erd_driver: largest symmetry",
     &        " forbidden two-electron SO integral", Asym
      Endif

      Do I = 1, 4
         If (Lu(I) .Gt. 0) Close(Lu(I),Status='Keep')
      Enddo

      Return
      End
