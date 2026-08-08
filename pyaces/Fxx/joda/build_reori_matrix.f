      Subroutine Build_reori_matrix(Qnew,Qold,AtmMass,W,Natoms)

      Implicit Double Precision(A-H, O-Z)

C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
      Dimension Qnew(3,Natoms) 
      Dimension Qold(3,Natoms) 
      Dimension AtmMass(Natoms) 
      Dimension Scr(6*Mxatms+18)
      Dimension W(3,3)

      Data Done /1.0D0/

      Print*, Natoms 
      Write(6,"(a)") " The new (sym. rotated) geometry"
      Do I = 1, Natoms
         Write(6,"(3(1x,F12.6))") (Qnew(j,i),j=1,3)
      Enddo 
      Write(6,"(a)") " The old (principal axis) geometry"
      Do I = 1, Natoms
         Write(6,"(3(1x,F12.6))") (Qold(j,i),j=1,3)
      Enddo 
C Qucik return for atoms and diatomics. The transfornmation is 
C the unit matrix.

      If (Natoms .Eq. 1 .OR. Natoms .Eq. 2) Then
         Call Zero(W,9)
         W(1,1) = Done 
         W(2,2) = Done 
         W(3,3) = Done  
         Return 
      Endif 

      Call Q2qprime(Qnew,Qold,AtmMass,Scr,W,Natoms)

      Write(6,"(a)") " Q(Sym. ordered)<->Q(Principle axis ordered) "
      Do j = 1,3
      Write(6,"(3(1x,F12.6))") (W(j,i),i=1,3)
      Enddo 

      Return
      End
