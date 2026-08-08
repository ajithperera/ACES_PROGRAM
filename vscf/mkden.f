











      SUBROUTINE MKDEN(COEF,DENS,NOCC,NBF,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION COEF(NBF,NOCC), DENS(NBF,NBF)
      DATA ONE /1.0/
      DATA TWO /2.0/
C
CSSS      Write(6,*) "The coeficient matrix @mkden:" 
CSSS      call  realprt(COEF,NBF*NOCC)

      CALL ZERO(DENS,NBF*NBF)
      CALL XGEMM('n','t',nbf,nbf,NOCC,
     &           ONE,COEF,nbf,
     &               COEF,nbf,
     &           ONE,DENS,nbf)
      IF (IUHF.EQ.0) CALL SSCAL(NBF*NBF,TWO,DENS,1)


      RETURN
      END

