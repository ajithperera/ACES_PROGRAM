      Double Precision Function Ovmono(Ia,Ib)
C Overlap of the monomials x^a y^b z^c (Ia) and Ib on one center with
C one exponent, up to a common factor: I(a+a') I(b+b') I(c+c') with
C I(n) = (n-1)!! for even n, 0 for odd n.
      Implicit Double Precision (A-H,O-Z)
      Integer Ia(3), Ib(3)
      Ovmono = 1.0D0
      Do K = 1, 3
         N = Ia(K) + Ib(K)
         If (Mod(N,2) .Eq. 1) Then
            Ovmono = 0.0D0
            Return
         Endif
         Ovmono = Ovmono*Dfact2(N/2)
      Enddo
      Return
      End
