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
#include "mxatms.par"
#include "erdints.com"

      Character*3 Kgen(3)
      Character*80 Xline
      Character*4 Atmnam
      Integer Ishl(8)
      Double Precision Salpha(200), Scoef(200,30)
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
