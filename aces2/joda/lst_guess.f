










      SUBROUTINE LST_GUESS(XR, XP, X, NATOMS, F)
C
C Transition state guess from a reactant XR and a product XP alone
C (QST2): the point at F (0.5) along the linear synchronous transit
C path of T. A. Halgren and W. N. Lipscomb, Chem. Phys. Lett. 49, 225
C (1977). The interatomic distances are interpolated,
C
C    r0_ij = (1-F) r_ij(XR) + F r_ij(XP),
C
C and X minimizes
C
C    S(X) = Sum_i<j (r_ij(X) - r0_ij)^2/r0_ij^4 + W |X - X0|^2,
C
C X0 = (1-F) XR + F XP after XP is moved onto XR by the best proper
C rotation (Kabsch); the small Cartesian term (W = 1.0D-6) fixes the
C overall position and orientation. Gauss-Newton on the normal
C equations, starting from X0, with step halving when S goes up.
C XP is returned rotated onto XR. Coordinates in bohr. 10/2026.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DIMENSION XR(3,NATOMS), XP(3,NATOMS), X(3,NATOMS)
      DIMENSION CR(3), CP(3), H(3,3), S(3), V(3,3), WT(3,3), U(3,3),
     &          WORK(64), D(3)
      DOUBLE PRECISION, ALLOCATABLE :: R0(:,:), X0(:,:), A(:,:),
     &                                 G(:), XT(:,:)
C
      PARAMETER (W = 1.0D-6, MAXIT = 200, TOL = 1.0D-10)
C
      NX = 3*NATOMS
      ALLOCATE(R0(NATOMS,NATOMS), X0(3,NATOMS), A(NX,NX), G(NX),
     &         XT(3,NATOMS))
C
C Rotate XP onto XR: H = Sum (XP-CP)(XR-CR)^t = V S W^t, U = W V^t,
C with the sign of the last singular vector flipped if det U < 0.
C
      DO I = 1, 3
         CR(I) = 0.0D0
         CP(I) = 0.0D0
         DO IATM = 1, NATOMS
            CR(I) = CR(I) + XR(I,IATM)
            CP(I) = CP(I) + XP(I,IATM)
         ENDDO
         CR(I) = CR(I)/DFLOAT(NATOMS)
         CP(I) = CP(I)/DFLOAT(NATOMS)
      ENDDO
      CALL ZERO(H, 9)
      DO IATM = 1, NATOMS
         DO J = 1, 3
            DO I = 1, 3
               H(I,J) = H(I,J) + (XP(I,IATM) - CP(I))*
     &                           (XR(J,IATM) - CR(J))
            ENDDO
         ENDDO
      ENDDO
      CALL DGESVD('A', 'A', 3, 3, H, 3, S, V, 3, WT, 3, WORK, 64,
     &            INFO)
      IF (INFO .NE. 0) THEN
         WRITE(6,"(a,I5)") " @LST_GUESS: DGESVD failed, INFO = ", INFO
         CALL ERREX
      ENDIF
      DO J = 1, 3
         DO I = 1, 3
            U(I,J) = WT(1,I)*V(J,1) + WT(2,I)*V(J,2) + WT(3,I)*V(J,3)
         ENDDO
      ENDDO
      DET = U(1,1)*(U(2,2)*U(3,3) - U(2,3)*U(3,2))
     &    - U(1,2)*(U(2,1)*U(3,3) - U(2,3)*U(3,1))
     &    + U(1,3)*(U(2,1)*U(3,2) - U(2,2)*U(3,1))
      IF (DET .LT. 0.0D0) THEN
         DO J = 1, 3
            DO I = 1, 3
               U(I,J) = U(I,J) - 2.0D0*WT(3,I)*V(J,3)
            ENDDO
         ENDDO
      ENDIF
      DO IATM = 1, NATOMS
         DO I = 1, 3
            D(I) = XP(I,IATM) - CP(I)
         ENDDO
         DO I = 1, 3
            XP(I,IATM) = CR(I) + U(I,1)*D(1) + U(I,2)*D(2) + U(I,3)*D(3)
         ENDDO
      ENDDO
C
C Interpolated distances and Cartesians.
C
      DO IATM = 1, NATOMS
         DO I = 1, 3
            X0(I,IATM) = (1.0D0-F)*XR(I,IATM) + F*XP(I,IATM)
         ENDDO
         DO JATM = 1, NATOMS
            R0(IATM,JATM) = (1.0D0-F)*DIST(XR(1,IATM), XR(1,JATM))
     &                    + F*DIST(XP(1,IATM), XP(1,JATM))
         ENDDO
      ENDDO
      CALL DCOPY(NX, X0, 1, X, 1)
C
      SOLD = SLST(X, X0, R0, NATOMS, W)
      DO ITER = 1, MAXIT
C
C Normal equations A dX = -G, A = J^t J + W, G = J^t e + W (X - X0),
C with e_ij = (r_ij - r0_ij)/r0_ij^2.
C
         CALL ZERO(A, NX*NX)
         DO I = 1, NX
            A(I,I) = W
         ENDDO
         DO IATM = 1, NATOMS
            DO I = 1, 3
               G(3*(IATM-1)+I) = W*(X(I,IATM) - X0(I,IATM))
            ENDDO
         ENDDO
         DO IATM = 2, NATOMS
            DO JATM = 1, IATM-1
               RIJ = DIST(X(1,IATM), X(1,JATM))
               R2  = R0(IATM,JATM)**2
               E   = (RIJ - R0(IATM,JATM))/R2
               DO I = 1, 3
                  D(I) = (X(I,IATM) - X(I,JATM))/(RIJ*R2)
               ENDDO
               II = 3*(IATM-1)
               JJ = 3*(JATM-1)
               DO I = 1, 3
                  G(II+I) = G(II+I) + E*D(I)
                  G(JJ+I) = G(JJ+I) - E*D(I)
                  DO J = 1, 3
                     T = D(I)*D(J)
                     A(II+I,II+J) = A(II+I,II+J) + T
                     A(JJ+I,JJ+J) = A(JJ+I,JJ+J) + T
                     A(II+I,JJ+J) = A(II+I,JJ+J) - T
                     A(JJ+I,II+J) = A(JJ+I,II+J) - T
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         CALL DPOSV('U', NX, 1, A, NX, G, NX, INFO)
         IF (INFO .NE. 0) THEN
            WRITE(6,"(a,I5)") " @LST_GUESS: DPOSV failed, INFO = ",
     &                        INFO
            CALL ERREX
         ENDIF
C
C G now holds -dX. Halve the step until S does not go up.
C
         STEP = 1.0D0
         DO IHALF = 1, 30
            DO IATM = 1, NATOMS
               DO I = 1, 3
                  XT(I,IATM) = X(I,IATM) - STEP*G(3*(IATM-1)+I)
               ENDDO
            ENDDO
            SNEW = SLST(XT, X0, R0, NATOMS, W)
            IF (SNEW .LE. SOLD) GOTO 10
            STEP = 0.5D0*STEP
         ENDDO
 10      CONTINUE
         CALL DCOPY(NX, XT, 1, X, 1)
         DXMAX = STEP*DABS(G(IDAMAX(NX, G, 1)))
         SOLD  = SNEW
         IF (DXMAX .LT. TOL) GOTO 20
      ENDDO
      WRITE(6,"(a,ES10.3)") "  The LST path point did not converge;"//
     &                      " last step ", DXMAX
 20   CONTINUE
C
      WRITE(6,"(a,F5.2,a,I4,a,ES10.3)") "  LST point at f =", F,
     &   " fitted in", MIN(ITER,MAXIT), " iterations, S =", SOLD
C
      DEALLOCATE(R0, X0, A, G, XT)
      RETURN
      END
