      Subroutine Erd_qtr(X,Y,Nl,Nf,Nmx,Nt,Iorb,Mlo,Wt)
C
C One quarter of the AO -> SO transformation. X(Nl,Nf,Nmx,Nt) holds a
C shell orbit index as (function f, member m), m = Mlo..Mlo+Nmx-1 of
C orbit Iorb; Y(Nl,Nf,Nn,Nt), Nn = orbit size, holds it as (f, j),
C the j-th SO made from function f:
C
C    Y(.,f,j,.) = Wt * sum_m s(m, SO(f,j)) X(.,f,m,.)
C
C Wt = Nn/Nmx makes the sum over a single member (Nmx = 1) stand for
C the sum over the whole orbit (valid for totally symmetric
C integrals only).
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"

      Dimension X(Nl,Nf,Nmx,Nt), Y(Nl,Nf,Norbm(Iorb),Nt)

      Nn = Norbm(Iorb)
      Do It = 1, Nt
         Do If = 1, Nf
            Do J = 1, Nn
               Ip = Lorbso(Iorbso(Iorb) + (If-1)*Nn + J)
               Do Il = 1, Nl
                  Y(Il,If,J,It) = 0.0D0
               Enddo
               Do K = 1, Nmx
                  W = Wt*Dble(Isosgn(Mlo+K-1,Ip))
                  Do Il = 1, Nl
                     Y(Il,If,J,It) = Y(Il,If,J,It) + W*X(Il,If,K,It)
                  Enddo
               Enddo
            Enddo
         Enddo
      Enddo

      Return
      End
