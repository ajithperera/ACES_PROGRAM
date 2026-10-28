










      SUBROUTINE LINDH_ADD(HC, BROW, F, IATM, NAT, NATOMS)
C
C HC <-- HC + F b b^t for a B-matrix row b (BROW) that is non-zero only
C on the NAT atoms IATM; BROW is zeroed on those atoms on return.
C Used by LINDH_HESS_CART. 10/2026.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION HC(3*NATOMS,3*NATOMS), BROW(3*NATOMS), IATM(NAT)
C
      DO IA = 1, NAT
         DO IX = 1, 3
            II = 3*IATM(IA) - 3 + IX
            DO JA = 1, NAT
               DO JX = 1, 3
                  JJ = 3*IATM(JA) - 3 + JX
                  HC(II,JJ) = HC(II,JJ) + F*BROW(II)*BROW(JJ)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      DO IA = 1, NAT
         DO IX = 1, 3
            BROW(3*IATM(IA) - 3 + IX) = 0.0D0
         ENDDO
      ENDDO
C
      RETURN
      END
