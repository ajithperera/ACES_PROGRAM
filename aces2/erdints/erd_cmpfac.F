      Subroutine Erd_cmpfac(L,Spher,Fac,Ncmp)
C
C VMOL/ERD scale factor of each component of an L shell (ERD order):
C ACES AO = Fac * ERD function. ERD normalizes every Cartesian
C component x^a y^b z^c separately; VMOL gives all of them the same
C radial factor (the one that normalizes xy, xyz, ...), which makes
C x^a y^b z^c larger by sqrt((2a-1)!!(2b-1)!!(2c-1)!!). Spherical
C functions are used as ERD makes them (Fac = 1).
C
      Implicit Double Precision (A-H,O-Z)
      Logical Spher
      Dimension Fac(*)

      If (Spher .And. L .Ge. 2) Then
         Ncmp = 2*L + 1
         Do I = 1, Ncmp
            Fac(I) = 1.0D0
         Enddo
         Return
      Endif

      Ncmp = 0
      Do Lx = L, 0, -1
         Do Ly = L-Lx, 0, -1
            Lz = L - Lx - Ly
            Ncmp = Ncmp + 1
            Fac(Ncmp) = Dsqrt(Dfact2(Lx)*Dfact2(Ly)*Dfact2(Lz))
         Enddo
      Enddo

      Return
      End

      Double Precision Function Dfact2(N)
C (2N-1)!!, with (-1)!! = 1.
      Implicit Double Precision (A-H,O-Z)
      Dfact2 = 1.0D0
      Do I = 1, 2*N-1, 2
         Dfact2 = Dfact2*Dble(I)
      Enddo
      Return
      End
