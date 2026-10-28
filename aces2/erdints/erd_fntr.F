      Subroutine Erd_fntr(X,Nl,Ldf,Nt,Is)
C
C Turn one index of X(Nl,Ldf,Nt) from ERD's functions of shell Is
C (the first Nfnsh(Is) positions of the middle dimension) into the
C basis functions the SOs are made of, in place: Cartesian functions
C are scaled by Fnfac, spherical ones (L >= 2) become VMOL's
C harmonics through Wvsph, contraction by contraction.
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"

      Dimension X(Nl,Ldf,Nt), Tmp(13)

      L = Lshl(Is)
      If (.Not. Spher .Or. L .Le. 1) Then
         I0 = Ioffn(Is)
         Do It = 1, Nt
            Do If = 1, Nfnsh(Is)
               F = Fnfac(I0+If)
               If (F .Ne. 1.0D0) Then
                  Do Il = 1, Nl
                     X(Il,If,It) = F*X(Il,If,It)
                  Enddo
               Endif
            Enddo
         Enddo
         Return
      Endif

      Ns = 2*L + 1
      Do It = 1, Nt
         Do Ic = 1, Nctsh(Is)
            I0 = (Ic-1)*Ns
            Do Il = 1, Nl
               Do J = 1, Ns
                  Sum = 0.0D0
                  Do M = 1, Ns
                     Sum = Sum + Wvsph(M,J,L)*X(Il,I0+M,It)
                  Enddo
                  Tmp(J) = Sum
               Enddo
               Do J = 1, Ns
                  X(Il,I0+J,It) = Tmp(J)
               Enddo
            Enddo
         Enddo
      Enddo

      Return
      End
