










      Subroutine Erg_bblk(Bx,N,Ng,Lcol,I0,Nf,Io,Nb,Out)
C
C Out(i,j) = Bx(I0+i, Lcol(Io+j)): the rows of one shell's functions
C and the columns of the SOs they enter.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Bx(N,Ng), Out(Nf,Nb)
      Integer Lcol(*)

      Do J = 1, Nb
         Ip = Lcol(Io+J)
         Do I = 1, Nf
            Out(I,J) = Bx(I0+I,Ip)
         Enddo
      Enddo

      Return
      End
