










      Subroutine Read_cart_constrains(Constrains, Num_constrains, 
     &                                Not_found)
C
      Implicit double precision (a-h,o-z)

      Character*80 wrk
      Integer Constrains
      Logical Not_found, ZMAT_PRESENT
C
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)

      Dimension Constrains(3*mxatms)
  
      Not_found = .False.

      INQUIRE(FILE='ZMAT',EXIST=ZMAT_PRESENT)
C
      IF (ZMAT_PRESENT)THEN
         OPEN(UNIT=4,FILE='ZMAT',FORM='FORMATTED',STATUS='OLD')
         rewind(4)
      else
         write(6,'(T3,a)')'@Read_constrains, ZMAT file not present.'
         call errex
      ENDIF
c
300   READ(4,'(A)', END=900) WRK
      IF (WRK(1:11) .NE.'*CONSTRAINS') goto 300

      READ(4,'(A)', END=900) WRK
      IF (WRK(1:11) .EQ.'*END') GOTO 999
      BACKSPACE(4)

      READ(4,*,END=900) Num_Constrains
     
          READ(4,*,END=900) (Constrains(Iread), Iread=1, Num_constrains)

      GO TO 99
C
900   WRITE(6,901)
901   FORMAT(3x,'@Read_constrains, *CONSTRAINS namelist not found or',
     &       ' incomplete.')
      GO TO 999
C
 999  CONTINUE

      CLOSE(UNIT=4,STATUS='KEEP')
      Write(6,*)
      Write(6,"(3x,2a)") "Constrained Cartesian ",
     &                   "optimization is requested but no constraints"
      Write(6,"(3x,a)")  "are specified in the inputfile (ZMAT)."
      Write(6,"(3x,2a)") "Calculation proceeds as a unconstrained ",
     &                   "optimization."
C
      Not_found = .True.

  99  CONTINUE
      Write(6,*)
      Write(6, "(a)") "The constrain table"
      Write(6, "(1x,6I3)") (Constrains(i), i=1, Num_constrains)
      Write(6,*) "Not found any constrains?",  Not_found
      CLOSE(UNIT=4,STATUS='KEEP')

      Return
      End
