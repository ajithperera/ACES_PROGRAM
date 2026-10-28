










      SUBROUTINE ALIGN_FRAME(XREF, XNEW, NATOMS, U, CREF, CNEW, RESID)
C
C Find the orthogonal matrix U (proper or improper) and the centroids
C CREF, CNEW that best map XREF onto XNEW,
C
C        XNEW_i - CNEW = U (XREF_i - CREF)
C
C (orthogonal Procrustes/Kabsch: H = Sum_i (XREF_i-CREF)(XNEW_i-CNEW)^t
C = V S W^t, U = W V^t). RESID is the largest |XNEW_i - CNEW -
C U(XREF_i - CREF)| (bohr).
C
C Used by Cartesian optimizations: the optimizer keeps its coordinates
C (R) in the frame of the first cycle while SYMMETRY reorients the
C geometry used for the integrals (Q) every cycle. U takes vectors
C from the optimizer frame to the integral frame. Improper U is
C allowed since SYMMETRY can change the handedness of the frame.
C 10/2026.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION XREF(3,NATOMS), XNEW(3,NATOMS), U(3,3), CREF(3),
     &          CNEW(3), H(3,3), S(3), V(3,3), WT(3,3), WORK(64)
C
      DO I = 1, 3
         CREF(I) = 0.0D0
         CNEW(I) = 0.0D0
      ENDDO
      DO IATM = 1, NATOMS
         DO I = 1, 3
            CREF(I) = CREF(I) + XREF(I,IATM)
            CNEW(I) = CNEW(I) + XNEW(I,IATM)
         ENDDO
      ENDDO
      DO I = 1, 3
         CREF(I) = CREF(I)/DFLOAT(NATOMS)
         CNEW(I) = CNEW(I)/DFLOAT(NATOMS)
      ENDDO
C
      CALL ZERO(H, 9)
      DO IATM = 1, NATOMS
         DO J = 1, 3
            DO I = 1, 3
               H(I,J) = H(I,J) + (XREF(I,IATM) - CREF(I))*
     &                           (XNEW(J,IATM) - CNEW(J))
            ENDDO
         ENDDO
      ENDDO
C
      LWORK = 64
      CALL DGESVD('A', 'A', 3, 3, H, 3, S, V, 3, WT, 3, WORK, LWORK,
     &            INFO)
      IF (INFO .NE. 0) THEN
         WRITE(6,"(a,I5)") " @ALIGN_FRAME: DGESVD failed, INFO = ", INFO
         CALL ERREX
      ENDIF
C
C U = W V^t = (WT)^t V^t
C
      DO J = 1, 3
         DO I = 1, 3
            U(I,J) = WT(1,I)*V(J,1) + WT(2,I)*V(J,2) + WT(3,I)*V(J,3)
         ENDDO
      ENDDO
C
      RESID = 0.0D0
      DO IATM = 1, NATOMS
         DO I = 1, 3
            D = XNEW(I,IATM) - CNEW(I)
            DO J = 1, 3
               D = D - U(I,J)*(XREF(J,IATM) - CREF(J))
            ENDDO
            RESID = MAX(RESID, DABS(D))
         ENDDO
      ENDDO
C
      Write(6,"(a,ES10.3)") " @ALIGN_FRAME: residual ", RESID
      Write(6,"(3F13.9)") ((U(I,J), J=1,3), I=1,3)
      RETURN
      END
