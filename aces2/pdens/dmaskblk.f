










      SUBROUTINE DMASKBLK(D, NMO, NOCC)
C
C DEBUG (_BLKCHK): KEEP ONE BLOCK OF THE MO-BASIS DENSITY, SELECTED BY
C IDMASK IN /DBLKMASK/: 0 ALL, 1 OO, 2 VV, 3 VO (VIRTUAL ROW, OCCUPIED
C COLUMN), 4 OV.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION D(NMO, NMO)
      COMMON /DBLKMASK/ IDMASK
C
      IF (IDMASK .EQ. 0) RETURN
      DO J = 1, NMO
         DO I = 1, NMO
            IB = 0
            IF (I .LE. NOCC .AND. J .LE. NOCC) IB = 1
            IF (I .GT. NOCC .AND. J .GT. NOCC) IB = 2
            IF (I .GT. NOCC .AND. J .LE. NOCC) IB = 3
            IF (I .LE. NOCC .AND. J .GT. NOCC) IB = 4
            IF (IB .NE. IDMASK) D(I, J) = 0.0D0
         ENDDO
      ENDDO
C
      RETURN
      END
