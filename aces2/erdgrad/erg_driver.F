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

      Integer Icore(Maxcor), Npr(8), Ioffg(8), Nbg(8)
      Character*3 Kgen(3)
      Logical Lecp, Lcorr, Debug

      Lecp = (Iflags(h_IFLAGS_ecp) .Ne. 0)
      Inquire(File='ERDGRAD.DEBUG',Exist=Debug)
      If (Iflags(h_IFLAGS_ref) .Gt. 2 .Or.
     &    Iflags(h_IFLAGS_qrhf_genera) .Ne. 0) Then
         Write(6,"(a)") " @Erg_driver: only RHF, UHF and ROHF"
     &                  //" references so far (no QRHF)"
         Call Errex
      Endif
      If (Iflags(h_IFLAGS_deriv_lev) .Gt. 1) Then
         Write(6,"(a)") " @Erg_driver: first derivatives only (use"
     &                  //" xvdint for DERIV_LEV=2)"
         Call Errex
      Endif
      Call Erd_rdmol(Kgen,Tla)
      Call Erd_symgrp(Kgen)
      Call Erd_sobld
      Write(6,"(a,i4,a,i5,a,i6,a,l2)") " @Erg_driver: atoms ", Natm,
     &     " shells ", Nshl, " functions ", Nfnerd, " spherical", Spher
      Write(6,"(a,i2,3(1x,a3))") " @Erg_driver: generators ", Ngen,
     &     (Kgen(I),I=1,Ngen)

C ECPs: effective nuclear charges and ecplib set up (Erg_ecpini).

      If (Lecp) Call Erg_ecpini(Icore(1))

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

C Densities over ERD functions (Erg_dens), kept with B to the end.

      Lcorr = (Iflags(h_IFLAGS_calc) .Ne. 0)
      N   = Nfnerd
      N2  = N*N*Iintfp
      I0b = J000
      I0ha = I0b  + N2
      I0hb = I0ha + N2
      I0xa = I0hb + N2
      I0xb = I0xa + N2
      I0t  = I0xb + N2
      I0w  = I0t  + N2
      I0d  = I0w  + N2
      I0c  = I0d  + N2
      I0o  = I0c  + N2
      I0e  = I0o  + N*Iintfp
      I0s  = I0e  + N*Iintfp
      Iend = I0s  + 2*N2
      If (Iend .Gt. Maxcor) Call Insmem("Erg_driver",Iend,Maxcor)
      Call Erg_bso(Icore(I0b),N)
      Call Erg_dens(Iuhf,Lcorr,N,Icore(I0b),Icore(I0ha),Icore(I0hb),
     &              Icore(I0xa),Icore(I0xb),Icore(I0t),Icore(I0w),
     &              Icore(I0d),Icore(I0c),Icore(I0o),Icore(I0e),
     &              Icore(I0s))

C Nonseparable two-particle density (correlated gradients only). The
C MOINTS lists are over the SOs of the SCF, or for spherical functions
C over the Cartesian SOs (T = CMP2CART takes the SCF SOs to them):
C their expansion in ERD functions is Bx = B, or B T+ with the
C pseudoinverse T+ = (T(T) T)**-1 T(T) (Gam(Cart) = T Gam T(T)).
C Then the SOs each shell enters (Erg_bcol), the SO Gamma (Erg_gam)
C and the buffers that transform one shell quartet of it (Erg_gnsp).

      Ng   = N
      I0bx = I0b
      Lgso = 0
      Lw   = 1
      I0rg = I0c
      I0nc = I0rg
      I0ic = I0rg
      I0lc = I0rg
      I0px = I0rg
      I0g2 = I0rg
      If (Lcorr) Then
         Do I = 1, Nirr
            Nbg(I) = Nbfirr(I)
         Enddo
         If (Spher) Then
            Call Getrec(20,'JOBARC','NAOBASFN',1,Ng)
            Call Getrec(20,'JOBARC','FAOBASIR',Nirr,Nbg)
         Endif
         I0rg = I0c
         I0bx = I0rg + Ng
         If (Spher) Then
            I0ct = I0bx + N*Ng*Iintfp
            I0m  = I0ct + Ng*N*Iintfp
            I0p  = I0m  + N*N*Iintfp
            I0iv = I0p  + N*Ng*Iintfp
            Iend = I0iv + N
            If (Iend .Gt. Maxcor) Call Insmem("Erg_driver",Iend,Maxcor)
            Call Getrec(20,'JOBARC','CMP2CART',Ng*N*Iintfp,Icore(I0ct))
            Call Erg_pinv(N,Ng,Icore(I0ct),Icore(I0m),Icore(I0p),
     &                    Icore(I0iv))
            Call Dgemm('N','N',N,Ng,N,1.0D0,Icore(I0b),N,Icore(I0p),N,
     &                 0.0D0,Icore(I0bx),N)
            I0nc = I0bx + N*Ng*Iintfp
         Else
            I0bx = I0b
            I0nc = I0rg + Ng
         Endif
         Ip = 0
         Do I = 1, Nirr
            Do J = 1, Nbg(I)
               Icore(I0rg+Ip) = I
               Ip = Ip + 1
            Enddo
         Enddo
         I0ic = I0nc + Nshl
         I0lc = I0ic + Nshl
         Call Erg_bcol(1,Icore(I0bx),N,Ng,Icore(I0nc),Icore(I0ic),
     &                 Icore(I0lc),Ltot)
         Call Erg_bcol(2,Icore(I0bx),N,Ng,Icore(I0nc),Icore(I0ic),
     &                 Icore(I0lc),Ltot)
         Mxcol = 0
         Do Is = 1, Nshl
            Mxcol = Max(Mxcol,Icore(I0nc+Is-1))
         Enddo
         Lw   = 2*Mxcol**4
         I0px = I0lc + Ltot
         I0g2 = I0px + Ng*Ng
         Call Erg_gam(1,Icore(I0px),Ng,Nirr,Nbg,Icore(I0rg),Npr,Ioffg,
     &                Lgso,Dum,Dum)
         Iend = I0g2 + (Lgso+Ng*Ng)*Iintfp
         If (Iend .Gt. Maxcor) Call Insmem("Erg_driver",Iend,Maxcor)
         Call Erg_gam(2,Icore(I0px),Ng,Nirr,Nbg,Icore(I0rg),Npr,Ioffg,
     &                Lgso,Icore(I0g2),Icore(I0g2+Lgso*Iintfp))
         Write(6,"(a,i12)") " @Erg_driver: SO Gamma elements stored ",
     &                      Lgso
      Endif
      I0w1 = I0g2 + Lgso*Iintfp
      I0w2 = I0w1 + Lw*Iintfp

C Gradient. Scratch for OED/ERD after the gradient arrays and the
C Gamma buffer of the largest shell quartet.

      Mxf  = 0
      Do Is = 1, Nshl
         Mxf = Max(Mxf,Nfnsh(Is))
      Enddo
      Lgam = Mxf**4
      I0g  = I0w2 + Lw*Iintfp
      I0gs = I0g  + 3*Natm*Iintfp
      I0g1 = I0gs + 3*Natm*Iintfp
      I0q  = I0g1 + 3*Natm*Iintfp
      I0z  = I0q  + Lgam*Iintfp
      Left = Maxcor - I0z + 1
      If (Left .Le. 0) Call Insmem("Erg_driver",I0z,Maxcor)
      Lzscr = (Left/2)/Iintfp
      I0i   = I0z + Lzscr*Iintfp
      Liscr = Maxcor - I0i + 1

C One-electron part and nuclear repulsion.

      Call Zero(Icore(I0g1),3*Natm)
      Call Erg_one(Icore(I0g1),Icore(I0t),Icore(I0w),N,
     &             Icore(I0z),Lzscr,Icore(I0i),Liscr)
      Call Erg_symz(Icore(I0g1),Icore(I0g))
      Call Erg_nucrep(Icore(I0g))
      If (Lecp) Then
         I0ec = I0q
         I0ep = I0ec + 3*Natoms*Iintfp
         I0ez = I0ep + Ncrt*Ncrt*Iintfp
         I0es = I0ez + Ncrt*N*Iintfp
         I0eg = I0es + Ncrt*N*Iintfp
         Iend = I0eg + 3*Natoms*Iintfp
         If (Iend .Gt. Maxcor) Call Insmem("Erg_driver",Iend,Maxcor)
         If (Debug) Then
            Write(6,"(a)") " @Erg_driver: one-electron + nuclear,"
     &                     //" before the ECP part"
            Call Erg_prgrd(Icore(I0g),Natm)
         Endif
         Call Erg_ecpgrd(Icore(I0g),Icore(I0d),N,Ncrt,Natoms,
     &                   Icore(K000),Icore(I0ec),Icore(I0ez),
     &                   Icore(I0ep),Icore(I0es),Icore(I0eg))
         If (Debug) Then
            Write(6,"(a)") " @Erg_driver: with the ECP part"
            Call Erg_prgrd(Icore(I0g),Natm)
         Endif
      Endif
      Call Dcopy(3*Natm,Icore(I0g),1,Icore(I0g1),1)

C Two-electron part. With the file ERDGRAD.DEBUG the separable and
C nonseparable parts are also computed and printed separately.

      Cns = 0.0D0
      If (Lcorr) Cns = 4.0D0
      If (Debug .And. Lcorr) Then
         Do Ipart = 1, 2
            Call Zero(Icore(I0g),3*Natm)
            Call Erg_two(Icore(I0g),Icore(I0ha),Icore(I0hb),
     &           Icore(I0xa),Icore(I0xb),N,Icore(I0z),Lzscr,
     &           Icore(I0i),Liscr,Icore(I0q),Lgam,Ipart.Eq.1,
     &           Dble(Ipart-1),Icore(I0bx),Ng,Icore(I0rg),
     &           Icore(I0nc),Icore(I0ic),Icore(I0lc),Icore(I0px),Npr,
     &           Ioffg,Icore(I0g2),Icore(I0w1),Icore(I0w2),Lw)
            Call Erg_symz(Icore(I0g),Icore(I0gs))
            If (Ipart .Eq. 1) Then
               Write(6,"(a)") " @Erg_driver: one-electron + nuclear"
               Call Erg_prgrd(Icore(I0g1),Natm)
               Write(6,"(a)") " @Erg_driver: separable two-electron"
            Else
               Write(6,"(a)") " @Erg_driver: nonseparable two-electron"
     &                        //" (factor 1)"
            Endif
            Call Erg_prgrd(Icore(I0gs),Natm)
         Enddo
      Endif

      Call Zero(Icore(I0g),3*Natm)
      Call Erg_two(Icore(I0g),Icore(I0ha),Icore(I0hb),Icore(I0xa),
     &             Icore(I0xb),N,Icore(I0z),Lzscr,Icore(I0i),Liscr,
     &             Icore(I0q),Lgam,.True.,Cns,Icore(I0bx),Ng,
     &             Icore(I0rg),Icore(I0nc),Icore(I0ic),Icore(I0lc),
     &             Icore(I0px),Npr,Ioffg,Icore(I0g2),Icore(I0w1),
     &             Icore(I0w2),Lw)
      Call Erg_symz(Icore(I0g),Icore(I0gs))
      Call Daxpy(3*Natm,1.0D0,Icore(I0g1),1,Icore(I0gs),1)

      Call Erg_wrgrd(Icore(I0gs),Natm,Natoms,Icore(K000),Icore(K010),
     &               Icore(I0z))

      Return
      End
