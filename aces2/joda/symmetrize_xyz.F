      SUBROUTINE SYMMETRIZE_XYZ(X, NATOMS, ERRMAX, APPLY)
C
C Project the Cartesian coordinates X onto the totally symmetric
C subspace of the full point group of the structure in /COORD/
C
C      x_P(i) <-- (1/|G|) Sum_g R_g x_i,  P(i) = image of atom i under g
C
C using the FULLSYOP/FULLPERM records that SYMDRV writes for the
C current geometry. X must be in the same frame as Q (the integral
C frame), which holds for the RIC back transformation in PROCESTEP_XYZ.
C The operations are rotated to the Q frame by ORIENTMT as in PROJFCM.
C
C ERRMAX is the largest change |X - SX| of any coordinate (bohr). With
C APPLY = .FALSE. X is not touched (use it to check that a structure is
C symmetric). When FULLORDR is 1 (C1 or SYMMETRY=NONE) nothing is done.
C 10/2026.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
#include "mxatms.par"
#include "machsp.com"
C
      PARAMETER (MXORDR = 120)
      LOGICAL APPLY
      DIMENSION X(3*NATOMS), SYOP(9*MXORDR), ORIENT(9),
     &          XSYM(3*MXATMS), IPTR(MXATMS*MXORDR)
C
      ERRMAX = 0.0D0
      IORDER = 0
      CALL GETREC(-1,'JOBARC','FULLORDR',1,IORDER)
      IF (IORDER .LE. 1) RETURN
      IF (IORDER .GT. MXORDR) THEN
         WRITE(6,"(a,I5)") " @SYMMETRIZE_XYZ: group order too large: ",
     &                       IORDER
         CALL ERREX
      ENDIF
C
      CALL GETREC(20,'JOBARC','FULLSYOP',9*IORDER*IINTFP,SYOP)
      CALL GETREC(20,'JOBARC','FULLPERM',NATOMS*IORDER,IPTR)
      CALL GETREC(20,'JOBARC','ORIENTMT',9*IINTFP,ORIENT)
      CALL TRNOPS(SYOP,ORIENT,IORDER)
C
      CALL ZERO(XSYM,3*NATOMS)
      DO IOP = 1, IORDER
         IOFF = 9*(IOP-1)
         DO IATM = 1, NATOMS
            JATM = IPTR(IATM + NATOMS*(IOP-1))
            IF (JATM .LT. 1 .OR. JATM .GT. NATOMS) THEN
               WRITE(6,"(a,3I5)") " @SYMMETRIZE_XYZ: bad atom image ",
     &                             IOP, IATM, JATM
               CALL ERREX
            ENDIF
            DO I = 1, 3
               S = 0.0D0
               DO J = 1, 3
                  S = S + SYOP(IOFF + I + 3*(J-1))*X(3*IATM-3+J)
               ENDDO
               XSYM(3*JATM-3+I) = XSYM(3*JATM-3+I) + S
            ENDDO
         ENDDO
      ENDDO
C
      FACT = 1.0D0/DFLOAT(IORDER)
      DO I = 1, 3*NATOMS
         XSYM(I) = FACT*XSYM(I)
         ERRMAX  = MAX(ERRMAX, DABS(XSYM(I) - X(I)))
      ENDDO
C
#ifdef _DEBUG_LVL0
      Write(6,"(a,I4,a,ES10.3)") " @SYMMETRIZE_XYZ: order ", IORDER,
     &                            " max change ", ERRMAX
#endif
      IF (APPLY) CALL DCOPY(3*NATOMS, XSYM, 1, X, 1)
C
      RETURN
      END
