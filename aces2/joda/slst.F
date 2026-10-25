      DOUBLE PRECISION FUNCTION SLST(X, X0, R0, NATOMS, W)
C
C The Halgren-Lipscomb LST objective of LST_GUESS,
C S = Sum_i<j (r_ij - r0_ij)^2/r0_ij^4 + W |X - X0|^2. 10/2026.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(3,NATOMS), X0(3,NATOMS), R0(NATOMS,NATOMS)
C
      SLST = 0.0D0
      DO IATM = 1, NATOMS
         DO I = 1, 3
            SLST = SLST + W*(X(I,IATM) - X0(I,IATM))**2
         ENDDO
         DO JATM = 1, IATM-1
            E = (DIST(X(1,IATM), X(1,JATM)) - R0(IATM,JATM))/
     &          R0(IATM,JATM)**2
            SLST = SLST + E*E
         ENDDO
      ENDDO
      RETURN
      END
