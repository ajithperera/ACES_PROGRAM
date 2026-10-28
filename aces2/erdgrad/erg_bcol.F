      Subroutine Erg_bcol(Mode,Bx,N,Ng,Ncol,Icol0,Lcol,Ltot)
C
C For every shell s the SOs (columns of Bx) its functions enter:
C Lcol(Icol0(s)+1 ... Icol0(s)+Ncol(s)). Mode 1 only counts (Ltot),
C mode 2 also fills Lcol.
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"

      Dimension Bx(N,Ng)
      Integer Ncol(Nshl), Icol0(Nshl), Lcol(*)

      Ltot = 0
      Do Is = 1, Nshl
         Icol0(Is) = Ltot
         Nc = 0
         Do Ip = 1, Ng
            Lnz = 0
            Do I = Ioffn(Is)+1, Ioffn(Is)+Nfnsh(Is)
               If (Bx(I,Ip) .Ne. 0.0D0) Lnz = 1
            Enddo
            If (Lnz .Eq. 1) Then
               Nc = Nc + 1
               If (Mode .Eq. 2) Lcol(Ltot+Nc) = Ip
            Endif
         Enddo
         Ncol(Is) = Nc
         Ltot = Ltot + Nc
      Enddo

      Return
      End
