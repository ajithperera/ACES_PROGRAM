      Subroutine Check_Strngs(A,B,D)

      Implicit None
      Integer*8 Len
      Integer*8 I

      Character*5 A
      Character*9 B
      Character*16 D

Cf2py intent (inout) D

      Write(6,"(2a)") "A: ", A
      Write(6,"(2a)") "B: ", B

      D = A//" "//B

      Write(6,"(2a)") "Concatenated: ", D

      Return 
      End

   
