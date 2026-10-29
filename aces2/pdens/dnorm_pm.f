










      SUBROUTINE DNORM_PM(DA, DB, NMO, NOCC)
C
C PRINTS NORMS OF DA+DB AND DA-DB FOR THE OO, VV, VO AND OV BLOCKS
C OF THE MO-BASIS PERTURBED DENSITIES (DEBUG ONLY).
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION DA(NMO, NMO), DB(NMO, NMO), SP(4), SM(4), XA(4), XB(4)
C
      CALL ZERO(SP, 4)
      CALL ZERO(SM, 4)
      CALL ZERO(XA, 4)
      CALL ZERO(XB, 4)
      DO J = 1, NMO
         DO I = 1, NMO
            IF (I .LE. NOCC .AND. J .LE. NOCC) K = 1
            IF (I .GT. NOCC .AND. J .GT. NOCC) K = 2
            IF (I .GT. NOCC .AND. J .LE. NOCC) K = 3
            IF (I .LE. NOCC .AND. J .GT. NOCC) K = 4
            SP(K) = SP(K) + (DA(I,J) + DB(I,J))**2
            SM(K) = SM(K) + (DA(I,J) - DB(I,J))**2
            XA(K) = XA(K) + DA(I,J)**2
            XB(K) = XB(K) + DB(I,J)**2
         ENDDO
      ENDDO
      WRITE(6, "(a)") " |DA+DB|, |DA-DB|, |DA|, |DB| for OO,VV,VO,OV:"
      DO K = 1, 4
         WRITE(6, "(1x,I2,4ES14.6)") K, SQRT(SP(K)), SQRT(SM(K)),
     &                               SQRT(XA(K)), SQRT(XB(K))
      ENDDO
C
      RETURN
      END
