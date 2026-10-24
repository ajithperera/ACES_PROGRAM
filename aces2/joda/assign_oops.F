      SUBROUTINE ASSIGN_OOPS(CARTCOORD, IREDUNCO, TOTREDNCO, TOTNOFBND,
     &                       TOTNOFANG, TOTNOFDIH, NOOPS, NRATMS,
     &                       MAXREDUNCO)
C
C Out-of-plane coordinates. At a planar three-coordinated center the
C three angles sum to 360 degrees and their B-rows are dependent, so
C bonds and angles alone do not describe pyramidalization (and there
C are no torsions at all in molecules like H2CO). For every atom C
C with exactly three bonds (to A, B and D) one improper torsion
C A-B-C-D is added at the end of the dihedral block; it is evaluated
C and differentiated as an ordinary dihedral (180 degrees when the
C center is planar). The set does not depend on the geometry, only on
C the bonds, so it does not change during an optimization. Of the
C three choices of B, the one with the angles A-C-B and B-C-D
C farthest from 0 and 180 degrees is taken. An improper that is also
C a proper torsion (three-membered rings) is not added, nor one with
C an angle at C within 10 degrees of 0 or 180 degrees. A file named
C NO_RIC_OOPS in the work directory turns them off. 10/2026.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      INTEGER TOTREDNCO, TOTNOFBND, TOTNOFANG, TOTNOFDIH
      DIMENSION CARTCOORD(3*NRATMS), IREDUNCO(4, MAXREDUNCO),
     &          NBR(3), V1(3), V2(3)
      PARAMETER (SINMIN = 0.17D0)
      LOGICAL NOOOPS
C
      NOOPS = 0
      INQUIRE(FILE='NO_RIC_OOPS', EXIST=NOOOPS)
      IF (NOOOPS) RETURN
      DO ICEN = 1, NRATMS
C
         NNBR = 0
         DO IBND = 1, TOTNOFBND
            JATM = 0
            IF (IREDUNCO(1, IBND) .EQ. ICEN) JATM = IREDUNCO(2, IBND)
            IF (IREDUNCO(2, IBND) .EQ. ICEN) JATM = IREDUNCO(1, IBND)
            IF (JATM .NE. 0) THEN
               NNBR = NNBR + 1
               IF (NNBR .LE. 3) NBR(NNBR) = JATM
            ENDIF
         ENDDO
         IF (NNBR .NE. 3) GOTO 10
C
C Pick B (the middle neighbor); A and D are the other two.
C
         BEST = -1.0D0
         DO IB = 1, 3
            IA = NBR(MOD(IB,   3) + 1)
            ID = NBR(MOD(IB+1, 3) + 1)
            SMIN = MIN(SINANG(IA, NBR(IB)), SINANG(NBR(IB), ID))
            IF (SMIN .GT. BEST) THEN
               BEST = SMIN
               IA_B = IA
               IB_B = NBR(IB)
               ID_B = ID
            ENDIF
         ENDDO
         IF (BEST .LT. SINMIN) GOTO 10
C
         DO IDIH = TOTNOFBND + TOTNOFANG + 1, TOTREDNCO
            IF ((IREDUNCO(1,IDIH) .EQ. IA_B .AND.
     &           IREDUNCO(2,IDIH) .EQ. IB_B .AND.
     &           IREDUNCO(3,IDIH) .EQ. ICEN .AND.
     &           IREDUNCO(4,IDIH) .EQ. ID_B) .OR.
     &          (IREDUNCO(1,IDIH) .EQ. ID_B .AND.
     &           IREDUNCO(2,IDIH) .EQ. ICEN .AND.
     &           IREDUNCO(3,IDIH) .EQ. IB_B .AND.
     &           IREDUNCO(4,IDIH) .EQ. IA_B)) GOTO 10
         ENDDO
C
         NOOPS = NOOPS + 1
         IF (TOTREDNCO + 1 .GE. MAXREDUNCO) THEN
            WRITE(6, 99) TOTREDNCO + 1, MAXREDUNCO
 99         FORMAT(T3,'Maximum No. of redundant coordinates allowed',
     &             ' exceeded.',/,
     &             T3,'Required:',I8,' Current Maximum:',I5)
            CALL ERREX
         ENDIF
         TOTREDNCO = TOTREDNCO + 1
         TOTNOFDIH = TOTNOFDIH + 1
         IREDUNCO(1, TOTREDNCO) = IA_B
         IREDUNCO(2, TOTREDNCO) = IB_B
         IREDUNCO(3, TOTREDNCO) = ICEN
         IREDUNCO(4, TOTREDNCO) = ID_B
C
 10      CONTINUE
      ENDDO
C
#ifdef _DEBUG_LVL0
      Write(6,"(a,I4)") " @ASSIGN_OOPS: out-of-plane coordinates: ",
     &                  NOOPS
#endif
      RETURN
C
      CONTAINS
C
C Sine of the angle I-ICEN-J.
C
      DOUBLE PRECISION FUNCTION SINANG(I, J)
      INTEGER I, J, K
      DOUBLE PRECISION C
      DO K = 1, 3
         V1(K) = CARTCOORD(3*I-3+K) - CARTCOORD(3*ICEN-3+K)
         V2(K) = CARTCOORD(3*J-3+K) - CARTCOORD(3*ICEN-3+K)
      ENDDO
      C = DDOT(3, V1, 1, V2, 1)/DSQRT(DDOT(3, V1, 1, V1, 1)*
     &                                DDOT(3, V2, 1, V2, 1))
      SINANG = DSQRT(MAX(0.0D0, 1.0D0 - C*C))
      END FUNCTION SINANG
      END
