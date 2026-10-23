










      Subroutine Erd_monos(L,Ixyz,N)
C The N = (L+1)(L+2)/2 monomials x^a y^b z^c of degree L in ERD order
C (a > b > c preference): Ixyz(1..3,i) = (a,b,c).
      Implicit Double Precision (A-H,O-Z)
      Integer Ixyz(3,*)
      N = 0
      Do Lx = L, 0, -1
         Do Ly = L-Lx, 0, -1
            N = N + 1
            Ixyz(1,N) = Lx
            Ixyz(2,N) = Ly
            Ixyz(3,N) = L - Lx - Ly
         Enddo
      Enddo
      Return
      End
