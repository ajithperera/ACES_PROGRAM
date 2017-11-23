      Subroutine reord_vecs(Vector, Scratch, Iorder, Length, Print)

      Implicit Double Precision (A-H, O-Z)
C
      Logical Print
 
      Dimension Vector(Length), Scratch(Length), Iorder(Length)
C
      Call Dcopy(Length, Vector, 1, Scratch, 1)
   

      Do I = 1, Length
      
         Vector(I) = Scratch(Iorder(i))

      Enddo


      Return
      End
