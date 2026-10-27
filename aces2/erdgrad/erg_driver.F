#include "flags.h"
      Subroutine Erg_driver(Icore,Maxcor,Iuhf)
C
C xerdgrad driver. The basis, symmetry and SO basis are built exactly
C as xerdints builds them (Erd_rdmol, Erd_symgrp, Erd_sobld,
C Erd_cmpord, Erd_soord), so the SOs are those the SCF vectors refer
C to. The gradient is then
C
C    nuclear repulsion + one-electron (Erg_one) + two-electron (Erg_two)
C
C over the atoms of MOL_ERDOED, symmetrized (Erg_symz) and written to
C the GRADIENT record in the computational (VMOL) center order, as
C xvdint writes it.
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"
#include "machsp.com"
#include "flags.com"

      Integer Icore(Maxcor)
      Character*3 Kgen(3)
      Logical Lecp

      Lecp = (Iflags(h_IFLAGS_ecp) .Ne. 0)
      If (Iflags(h_IFLAGS_calc) .Ne. 0 .Or.
     &    Iflags(h_IFLAGS_ref) .Gt. 1) Then
         Write(6,"(a)") " @Erg_driver: only RHF and UHF SCF gradients"
     &                  //" so far"
         Call Errex
      Endif
      If (Lecp) Then
         Write(6,"(a)") " @Erg_driver: no ECP gradients yet"
         Call Errex
      Endif

      Call Erd_rdmol(Kgen,Tla)
      Call Erd_symgrp(Kgen)
      Call Erd_sobld
      Write(6,"(a,i4,a,i5,a,i6,a,l2)") " @Erg_driver: atoms ", Natm,
     &     " shells ", Nshl, " functions ", Nfnerd, " spherical", Spher
      Write(6,"(a,i2,3(1x,a3))") " @Erg_driver: generators ", Ngen,
     &     (Kgen(I),I=1,Ngen)

C Computational order of the AOs and SOs (as in Erd_driver). The
C ZMAT center -> atom map (Iz2a) and computational center -> ZMAT
C center map (Imap) are kept at K000 and K010.

      Call Getrec(20,'JOBARC','NATOMS  ',1,Natoms)
      Ncrt = Ncart
      K000 = 1
      K010 = K000 + Natoms
      J000 = K010 + Natoms
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
      If (J130 .Gt. Maxcor) Call Insmem("Erg_driver",J130,Maxcor)
      Call Erd_cmpord(Icore(J000),Ncrt,Icore(J010),Natoms,
     &                Icore(J020),Icore(J030),Icore(J040),Icore(J050),
     &                Icore(J060),Icore(J070),Icore(J080),Icore(J090),
     &                Icore(J100),Icore(J110),Icore(J120))
      Do I = 0, Natoms-1
         Icore(K000+I) = Icore(J010+I)
         Icore(K010+I) = Icore(J040+I)
      Enddo
      If (J010+13*Nso .Gt. Maxcor) Call Insmem("Erg_driver",
     &                                         J010+13*Nso,Maxcor)
      Call Erd_soord(Icore(J000),Ncrt,Icore(J010),Icore(J010+Nso),
     &               Icore(J010+2*Nso),11*Nso)

C Densities over ERD functions.

      N   = Nfnerd
      N2  = N*N*Iintfp
      I0b = J000
      I0a = I0b + N2
      I0p = I0a + N2
      I0w = I0p + N2
      I0t = I0w + N2
      I0c = I0t + N2
      I0o = I0c + N2
      I0e = I0o + N*Iintfp
      I0s = I0e + N*Iintfp
      I0g = I0s + 2*N2
      If (I0g .Gt. Maxcor) Call Insmem("Erg_driver",I0g,Maxcor)
      Call Erg_bso(Icore(I0b),N)
      Call Erg_dens(Iuhf,N,Icore(I0b),Icore(I0a),Icore(I0p),
     &              Icore(I0w),Icore(I0c),Icore(I0o),Icore(I0e),
     &              Icore(I0s))
      Call Dcopy(N*N,Icore(I0a),1,Icore(I0t),1)
      Call Daxpy(N*N,1.0D0,Icore(I0p),1,Icore(I0t),1)

C Gradient. Scratch for OED/ERD after the gradient arrays and the
C Gamma buffer of the largest shell quartet.

      I0gs = I0c
      I0g  = I0gs + 3*Natm*Iintfp
      Mxf  = 0
      Do Is = 1, Nshl
         Mxf = Max(Mxf,Nfnsh(Is))
      Enddo
      Lgam = Mxf**4
      I0q  = I0g + 3*Natm*Iintfp
      I0z  = I0q + Lgam*Iintfp
      Left = Maxcor - I0z + 1
      If (Left .Le. 0) Call Insmem("Erg_driver",I0z,Maxcor)
      Lzscr = (Left/2)/Iintfp
      I0i   = I0z + Lzscr*Iintfp
      Liscr = Maxcor - I0i + 1

      Call Zero(Icore(I0g),3*Natm)
      Call Erg_one(Icore(I0g),Icore(I0t),Icore(I0w),N,
     &             Icore(I0z),Lzscr,Icore(I0i),Liscr)
      Call Erg_two(Icore(I0g),Icore(I0t),Icore(I0a),Icore(I0p),N,
     &             Icore(I0z),Lzscr,Icore(I0i),Liscr,Icore(I0q),Lgam)
      Call Erg_symz(Icore(I0g),Icore(I0gs))
      Call Erg_nucrep(Icore(I0gs))

      Call Erg_wrgrd(Icore(I0gs),Natm,Natoms,Icore(K000),Icore(K010),
     &               Icore(I0z))

      Return
      End
