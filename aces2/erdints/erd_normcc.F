      Subroutine Erd_normcc(Np,Nc,L,Alp,Coef,Ldc)
C
C Rescale the contraction coefficients of one shell so that each
C contracted function is normalized, assuming normalized primitives
C (ERD/OED supply the primitive normalization). The overlap of two
C normalized primitives of angular momentum L is
C (2 sqrt(a b)/(a+b))**(L+3/2). Same as READ_BASIS_INFO (a3ecp).
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Alp(Np), Coef(Ldc,Nc)

      Pw = Dble(L) + 1.5D0
      Do J = 1, Nc
         Sum = 0.0D0
         Do I = 1, Np
            Do K = 1, Np
               Ovl = (2.0D0*Dsqrt(Alp(I)*Alp(K))/(Alp(I)+Alp(K)))**Pw
               Sum = Sum + Coef(I,J)*Coef(K,J)*Ovl
            Enddo
         Enddo
         Fac = 1.0D0/Dsqrt(Sum)
         Do I = 1, Np
            Coef(I,J) = Coef(I,J)*Fac
         Enddo
      Enddo

      Return
      End
