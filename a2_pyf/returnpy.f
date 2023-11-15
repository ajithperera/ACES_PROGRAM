      Subroutine Returnpy(Destn,Ndim)

C This is an experiemntal routine that test in/out python to Fortran
C array.

      Implicit None 

      Integer*8 Ndim
      Integer*8 I
      Double Precision Destn(*)

      Do I = 1, Ndim
         Destn(I) = 1.0D0
      Enddo 

      write(6,"(a)") " The destination array @-returnpy"
      write(6,"(6(1x,F10.5))") (Destn(i),i=1,Ndim)

      Return
      End 
