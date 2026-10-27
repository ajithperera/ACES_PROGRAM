










      Subroutine Erg_prgrd(G,Natm)
C
C Print a gradient, atoms of MOL_ERDOED.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension G(3,Natm)

      Do Ia = 1, Natm
         Write(6,"(2x,i4,3f20.12)") Ia, (G(K,Ia),K=1,3)
      Enddo

      Return
      End
