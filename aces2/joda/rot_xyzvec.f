










      SUBROUTINE ROT_XYZVEC(V, NATOMS, NVEC, LDV, U, TRANS)
C
C Rotate the Cartesian 3-vectors of NVEC columns of V (each column
C is 3*NATOMS long, leading dimension LDV): v_i <-- U v_i, or
C v_i <-- U^t v_i when TRANS is .TRUE. 10/2026.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL TRANS
      DIMENSION V(LDV,NVEC), U(3,3), W(3)
C
      DO IVEC = 1, NVEC
         DO IATM = 1, NATOMS
            IOFF = 3*(IATM-1)
            DO I = 1, 3
               W(I) = 0.0D0
               DO J = 1, 3
                  IF (TRANS) THEN
                     W(I) = W(I) + U(J,I)*V(IOFF+J,IVEC)
                  ELSE
                     W(I) = W(I) + U(I,J)*V(IOFF+J,IVEC)
                  ENDIF
               ENDDO
            ENDDO
            DO I = 1, 3
               V(IOFF+I,IVEC) = W(I)
            ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END
