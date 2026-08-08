










      Subroutine Build_cmatrix(Nxm6, Cmat)

      Implicit Double Precision (A-H, O-Z)

C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
      Integer Constrained 
      Dimension Cmat(Nxm6, Nxm6), Constrained(Maxredunco)
C
      Nxm62 = Nxm6*Nxm6
      Call Dzero(Cmat, Nxm62)
C
C The constrained designates the constrains. It is prepared and
C archived in gen_redinternals. 
C
      Call Getrec( 0, "JOBARC", "CONSTRNS", Ilength, Ijunk)
      Call Getrec(20, "JOBARC", "CONSTRNS", Ilength, Constrained)

      write(6,"(a)") "The constrained array"
      Write(6,"(6(1x,i4))") (Constrained(i),i=1, Ilength)
      Do Idim =1, Nxm6
      
         If (constrained(Idim) .EQ. 1) Then
             Cmat(Idim, Idim) = 1.0D0
         Endif

      Enddo

      write(6,"(a)") "The constrain Matrix"
      call output(Cmat, 1, Nxm6, 1, Nxm6, Nxm6, Nxm6, 1)
      Return
      End
