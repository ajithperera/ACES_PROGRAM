










      SUBROUTINE RIC_FDCHECK(Q, BT, TOTREDNCO, NATOMS, IATNUM)
C
C Consistency check of the RIC B-matrix: compare the analytic
C B(i,k) = dq_i/dx_k (stored transposed, BT(3*NATOMS,TOTREDNCO), as
C left by BUILT_BGMTRX) with central differences of the RIC values
C computed by GEN_NEW_RIC, the routine used in the back
C transformation. Runs only when a file named RIC_FDCHECK exists in
C the work directory. The second component of a linear bend
C (IREDUNCO(4,i) = -1) has a value fixed at 180 degrees in
C GEN_NEW_RIC, so it is reported separately. 10/2026.
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
      INTEGER TOTREDNCO, TOTNOFBND, TOTNOFANG, TOTNOFDIH
      LOGICAL DOIT
      PARAMETER (H = 1.0D-4)
      DIMENSION Q(3*NATOMS), BT(3*NATOMS,TOTREDNCO), IATNUM(NATOMS),
     &          X(3*MXATMS), QP(MAXREDUNCO), QM(MAXREDUNCO),
     &          IREDUNCO(4,MAXREDUNCO), ERRTYP(4), IWORST(2,4)
C
      INQUIRE(FILE='RIC_FDCHECK', EXIST=DOIT)
      IF (.NOT. DOIT) RETURN
C
      PI = DACOS(-1.0D0)
      CALL GETREC(20, 'JOBARC', 'CONTEVIT', 4*TOTREDNCO, IREDUNCO)
      CALL GETREC(20, 'JOBARC', 'TNUMOBND', 1, TOTNOFBND)
      CALL GETREC(20, 'JOBARC', 'TNUMOANG', 1, TOTNOFANG)
      CALL GETREC(20, 'JOBARC', 'TNUMODIH', 1, TOTNOFDIH)
C
      DO I = 1, 4
         ERRTYP(I)   = 0.0D0
         IWORST(1,I) = 0
         IWORST(2,I) = 0
      ENDDO
C
      DO K = 1, 3*NATOMS
         CALL DCOPY(3*NATOMS, Q, 1, X, 1)
         X(K) = Q(K) + H
         CALL GEN_NEW_RIC(X, QP, IATNUM, NATOMS, NB, NW, NA, .FALSE.)
         X(K) = Q(K) - H
         CALL GEN_NEW_RIC(X, QM, IATNUM, NATOMS, NB, NW, NA, .FALSE.)
         DO I = 1, TOTREDNCO
            D = QP(I) - QM(I)
            IF (I .GT. TOTNOFBND+TOTNOFANG) THEN
               IF (D .GT.  PI) D = D - 2.0D0*PI
               IF (D .LT. -PI) D = D + 2.0D0*PI
               ITYP = 3
            ELSE IF (I .GT. TOTNOFBND) THEN
               ITYP = 2
               IF (IREDUNCO(4,I) .EQ. -1) ITYP = 4
            ELSE
               ITYP = 1
            ENDIF
            ERR = DABS(D/(2.0D0*H) - BT(K,I))
            IF (ERR .GT. ERRTYP(ITYP)) THEN
               ERRTYP(ITYP)   = ERR
               IWORST(1,ITYP) = I
               IWORST(2,ITYP) = K
            ENDIF
         ENDDO
      ENDDO
C
      WRITE(6,*)
      WRITE(6,"(a,ES9.2,a)") " @RIC_FDCHECK: max |B - dq/dx| (h =",
     &                       H, " bohr)"
      WRITE(6,"(a,ES10.3,a,2I5)") "   bonds            ", ERRTYP(1),
     &                            "  ric, x:", IWORST(1,1), IWORST(2,1)
      WRITE(6,"(a,ES10.3,a,2I5)") "   angles           ", ERRTYP(2),
     &                            "  ric, x:", IWORST(1,2), IWORST(2,2)
      WRITE(6,"(a,ES10.3,a,2I5)") "   dihedrals        ", ERRTYP(3),
     &                            "  ric, x:", IWORST(1,3), IWORST(2,3)
      WRITE(6,"(a,ES10.3,a,2I5)") "   linear bend (2nd)", ERRTYP(4),
     &                            "  ric, x:", IWORST(1,4), IWORST(2,4)
      DO ITYP = 1, 4
         I = IWORST(1,ITYP)
         IF (I .GT. 0) WRITE(6,"(a,I2,a,4I4)") "   worst type ", ITYP,
     &                 " RIC atoms:", (IREDUNCO(J,I), J=1,4)
      ENDDO
      WRITE(6,*)
C
      RETURN
      END
