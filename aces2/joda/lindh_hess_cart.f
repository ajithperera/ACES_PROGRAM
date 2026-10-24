










      SUBROUTINE LINDH_HESS_CART(X, HC, NATOMS, IATNUM)
C
C Lindh model Hessian in Cartesian coordinates (R. Lindh et al., Chem.
C Phys. Lett. 241, 423 (1995)): H = Sum k_c b_c b_c^t over all atom
C pairs (bonds, k = 0.45 rho_ij), all triples (angles i-j-k, k = 0.15
C rho_ij rho_jk) and all quadruples (dihedrals i-j-k-l, k = 0.005 
C rho_ij rho_jk rho_kl), b_c the Wilson B-matrix row of the coordinate
C (BULT_BNDCRD/BULT_ANGCRD/BULT_DIHANGCRD). Terms with a weight below
C RHOMIN are skipped, and so are angles near 0 degrees and dihedrals 
C with a nearly linear angle; a nearly linear angle gets both bends
C (BULT_LANGCRD for the second). X (bohr) must be in the frame of the optimizer. 10/2026.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
C
      PARAMETER (RHOMIN = 1.0D-8, SIN2MIN = 1.0D-2)
      DIMENSION X(3*NATOMS), HC(3*NATOMS,3*NATOMS), IATNUM(NATOMS),
     &          BROW(3*MXATMS), RHO(MXATMS,MXATMS), IATM(4)
      DOUBLE PRECISION LINDH_RHO
C
      CALL ZERO(HC, 9*NATOMS*NATOMS)
      CALL ZERO(BROW, 3*NATOMS)
      DO IA = 1, NATOMS
         RHO(IA,IA) = 0.0D0
         DO IB = 1, IA-1
            R2 = (X(3*IA-2)-X(3*IB-2))**2 + (X(3*IA-1)-X(3*IB-1))**2 +
     &           (X(3*IA  )-X(3*IB  ))**2
            RHO(IA,IB) = LINDH_RHO(IATNUM(IA), IATNUM(IB), R2)
            RHO(IB,IA) = RHO(IA,IB)
         ENDDO
      ENDDO
C
C Bonds
C
      DO I = 2, NATOMS
         DO J = 1, I-1
            F = 0.45D0*RHO(I,J)
            IF (F .GE. RHOMIN) THEN
               CALL BULT_BNDCRD(X, BROW, D, I, J, 1, 1, NATOMS)
               IATM(1) = I
               IATM(2) = J
               CALL LINDH_ADD(HC, BROW, F, IATM, 2, NATOMS)
            ENDIF
         ENDDO
      ENDDO
C
C Angles I-J-K (vertex J)
C
      DO J = 1, NATOMS
         DO I = 1, NATOMS
            DO K = I+1, NATOMS
               IF (I.NE.J .AND. K.NE.J) THEN
                  F = 0.15D0*RHO(I,J)*RHO(J,K)
                  S2 = SIN2(X,I,J,K)
                  IF (F .GE. RHOMIN .AND. (S2 .GE. SIN2MIN .OR.
     &                COSA(X,I,J,K) .LT. 0.0D0)) THEN
                     CALL BULT_ANGCRD(X, BROW, ANG, I, J, K, 1, 1,
     &                                NATOMS)
                     IATM(1) = I
                     IATM(2) = J
                     IATM(3) = K
                     CALL LINDH_ADD(HC, BROW, F, IATM, 3, NATOMS)
C
C A nearly linear I-J-K bends in two directions; add the second one
C
                     IF (S2 .LT. SIN2MIN) THEN
                        CALL BULT_LANGCRD(X, BROW, ANG, I, J, K, -1, 1,
     &                                    NATOMS, 1, 1.0D-8)
                        CALL LINDH_ADD(HC, BROW, F, IATM, 3, NATOMS)
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
C
C Dihedrals I-J-K-L; (I,J,K,L) and (L,K,J,I) are the same, so J < K
C
      DO J = 1, NATOMS
         DO K = J+1, NATOMS
            IF (RHO(J,K) .GE. RHOMIN) THEN
            DO I = 1, NATOMS
               IF (I.NE.J .AND. I.NE.K) THEN
               DO L = 1, NATOMS
                  IF (L.NE.I .AND. L.NE.J .AND. L.NE.K) THEN
                     F = 0.005D0*RHO(I,J)*RHO(J,K)*RHO(K,L)
                     IF (F .GE. RHOMIN) THEN
                        IF (SIN2(X,I,J,K) .GE. SIN2MIN .AND.
     &                      SIN2(X,J,K,L) .GE. SIN2MIN) THEN
                           CALL BULT_DIHANGCRD(X, BROW, DANG, I, J, K,
     &                                         L, 1, 1, NATOMS)
                           IATM(1) = I
                           IATM(2) = J
                           IATM(3) = K
                           IATM(4) = L
                           CALL LINDH_ADD(HC, BROW, F, IATM, 4, NATOMS)
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
               ENDIF
            ENDDO
            ENDIF
         ENDDO
      ENDDO
C
      RETURN
C
      CONTAINS
C
C Squared sine of the angle A-B-C
C
      DOUBLE PRECISION FUNCTION SIN2(X, IA, IB, IC)
      DOUBLE PRECISION X(*), U(3), V(3), UU, VV, UV
      INTEGER IA, IB, IC, M
      UU = 0.0D0
      VV = 0.0D0
      UV = 0.0D0
      DO M = 1, 3
         U(M) = X(3*IA-3+M) - X(3*IB-3+M)
         V(M) = X(3*IC-3+M) - X(3*IB-3+M)
         UU = UU + U(M)*U(M)
         VV = VV + V(M)*V(M)
         UV = UV + U(M)*V(M)
      ENDDO
      SIN2 = 1.0D0 - UV*UV/(UU*VV)
      END FUNCTION SIN2
C
C Cosine of the angle A-B-C
C
      DOUBLE PRECISION FUNCTION COSA(X, IA, IB, IC)
      DOUBLE PRECISION X(*), UU, VV, UV, U, V
      INTEGER IA, IB, IC, M
      UU = 0.0D0
      VV = 0.0D0
      UV = 0.0D0
      DO M = 1, 3
         U  = X(3*IA-3+M) - X(3*IB-3+M)
         V  = X(3*IC-3+M) - X(3*IB-3+M)
         UU = UU + U*U
         VV = VV + V*V
         UV = UV + U*V
      ENDDO
      COSA = UV/DSQRT(UU*VV)
      END FUNCTION COSA
C
      END
