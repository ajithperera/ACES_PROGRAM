      Subroutine Erg_nucrep(G)
C
C Add the nuclear repulsion gradient, dE/dR(a) = - sum_b Z(a) Z(b)
C (R(a) - R(b)) / |R(a) - R(b)|**3.
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"

      Dimension G(3,Natm)

      Do Ia = 1, Natm
         Do Ib = 1, Natm
            If (Ib .Ne. Ia) Then
               R2 = 0.0D0
               Do K = 1, 3
                  R2 = R2 + (Xyzatm(K,Ia)-Xyzatm(K,Ib))**2
               Enddo
               F = Zatm(Ia)*Zatm(Ib)/(R2*Dsqrt(R2))
               Do K = 1, 3
                  G(K,Ia) = G(K,Ia) - F*(Xyzatm(K,Ia)-Xyzatm(K,Ib))
               Enddo
            Endif
         Enddo
      Enddo

      Return
      End
