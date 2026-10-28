










      Subroutine Mix_pairs(Scf_mos,Eps,Index1,Index2,Ip1Sym,Ip2sym,
     +                     Nbasis,Type)

      Implicit Double Precision(A-H,O-Z)
      Character*2 Type 

      Dimension Scf_mos(Nbasis,Nbasis)

C p' = eps (p + q), q' = eps (p - q), then normalized, for every pair
C type (OO, VV, OV, VO). OV used to store the results in swapped
C columns. p and q are orthonormal SCF orbitals, so p +/- q are
C orthogonal with norm sqrt(2); normalizing leaves (p +/- q)/sqrt(2)
C times the sign of eps, whatever its magnitude.

      If (Dabs(Eps) .Lt. 1.0D-12) Then
         Write(6,"(a)") " Mix_pairs: eps = 0 gives zero orbitals."
         Call Errex
      Endif
      Scale = Sign(1.0D0/Dsqrt(2.0D0),Eps)

      Do Irow = 1, Nbasis

         P = Scf_mos(Irow,Index1)
         Q = Scf_mos(Irow,Index2)

         Scf_mos(Irow,Index1) = Scale*(P + Q)

         Scf_mos(Irow,Index2) = Scale*(P - Q)

      Enddo

      Return
      End
