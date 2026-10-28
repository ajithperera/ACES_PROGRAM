      Subroutine Erd_lapl(L,P,Q)
C Q = Laplacian of the polynomial P over the degree L monomials (ERD
C order, Erd_monos); Q is over the degree L-2 monomials. The monomial
C (a,b,c) of degree d sits at s(s+1)/2 + (s-b) + 1, s = d-a.
      Implicit Double Precision (A-H,O-Z)
      Dimension P(*), Q(*)
      Indx(Ia,Ib,Id) = (Id-Ia)*(Id-Ia+1)/2 + (Id-Ia-Ib) + 1
      Do K = 1, (L-1)*L/2
         Q(K) = 0.0D0
      Enddo
      Do Lx = L, 0, -1
         Do Ly = L-Lx, 0, -1
            Lz = L - Lx - Ly
            X = P(Indx(Lx,Ly,L))
            If (Lx .Ge. 2) Then
               K = Indx(Lx-2,Ly,L-2)
               Q(K) = Q(K) + Dble(Lx*(Lx-1))*X
            Endif
            If (Ly .Ge. 2) Then
               K = Indx(Lx,Ly-2,L-2)
               Q(K) = Q(K) + Dble(Ly*(Ly-1))*X
            Endif
            If (Lz .Ge. 2) Then
               K = Indx(Lx,Ly,L-2)
               Q(K) = Q(K) + Dble(Lz*(Lz-1))*X
            Endif
         Enddo
      Enddo
      Return
      End
