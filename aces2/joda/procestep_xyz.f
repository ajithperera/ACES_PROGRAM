










      SUBROUTINE PROCESTEP_XYZ(SCRATCH, AMATRX, TOTREDNCO, STEP)
C
C Back transform the RIC step to Cartesians (Pulay et al., JCP 96, 2856,
C 1992):  x_k+1 = x_k + A(dq0 - (q(x_k) - q0)), A = B^t(BB^t)^-, with A
C fixed at the starting geometry.
C
C On return STEP holds the step actually taken, q(x_new) - q0, which
C is archived and used in the next Hessian update. The requested dq0
C is not reproduced exactly (convergence threshold, best iterate or
C linear fallback, rounding and symmetrization). 10/2026.
C
C Symmetry: the optimizer works with the full (C1) set of RICs and knows
C nothing about symmetry; the electronic structure runs in the point
C group of the starting structure. The point group is preserved by
C projecting every Cartesian iterate onto the totally symmetric
C subspace of the full point group (SYMMETRIZE_XYZ). This replaces the
C equivalence/plus-minus pair logic of PROCESS_REDINTRNLS. 10/2026.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL CONVERGED, CHANGED, NO_ITER, DO_SYM
      INTEGER TOTREDNCO, TOTNOFBND, TOTNOFANG
      DOUBLE PRECISION ORIENT(3,3)
      PARAMETER (MAXITER = 100)
C
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)


c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end







































































































































































































c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
C coord.com : begin
C
      DOUBLE PRECISION Q, R, ATMASS
      INTEGER NCON, NR, ISQUASH, IATNUM, IUNIQUE, NEQ, IEQUIV,
     &        NOPTI, NATOMS
      COMMON /COORD/ Q(3*MXATMS), R(MAXREDUNCO), NCON(MAXREDUNCO),
     &     NR(MXATMS),ISQUASH(MAXREDUNCO),IATNUM(MXATMS),
     &     ATMASS(MXATMS),IUNIQUE(MAXREDUNCO),NEQ(MAXREDUNCO),
     &     IEQUIV(MAXREDUNCO,MAXREDUNCO),
     &     NOPTI(MAXREDUNCO), NATOMS

C coord.com : end


C
      PARAMETER (EPSILON = 1.0D-5, SYMTHRS = 1.0D-5)
      CHARACTER*80 HOLDQ
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      DIMENSION SCRATCH(9*NATOMS*NATOMS), IQFIX(3*MXATMS, 3*MXATMS),
     &          AMATRX(3*NATOMS*TOTREDNCO), STATS(6),
     &          REDUNCO(MAXREDUNCO), SCRATCH_LCL(4*MAXREDUNCO+6*
     &          MXATMS+15*MXATMS), XBEST(3*MXATMS), STEP(TOTREDNCO)
      DIMENSION TMP1(MAXREDUNCO*MAXREDUNCO),TMP2(MAXREDUNCO)
C
      DATA BTOA /0.529177249D0/

      DO I = 1, NUNIQUE
         DO J = 1, NEQ(IUNIQUE(I))
            IQFIX(IUNIQUE(I), J) = IEQUIV(I,J)
         ENDDO
      ENDDO
C
      DO J = 1, NOPT
C
         SCRATCH_LCL(NOPTI(J)) = SCRATCH(J + NOPT)
         SCRATCH_LCL(NOPT + NOPTI(J)) = SCRATCH(J + NOPT)
C
         DO K = 1, NEQ(NOPTI(J))
            SCRATCH_LCL(IQFIX(NOPTI(J), K)) = SCRATCH(J + NOPT)
           SCRATCH_LCL(NOPT + IQFIX(NOPTI(J), K)) = SCRATCH(J + NOPT)
         ENDDO
      ENDDO
C
      Write(6,*)
      WRITE(6,*) "The optimization parameters:nx, nopt, totrednco nxm6",
     &            nx, totrednco, nopt, nxm6
      Write(6,*)
      WRITE(6,*) "The delta q0"
      Write(6,*)
      WRITE(6,10) (SCRATCH_LCL(I), I = 1, TOTREDNCO)
      Write(6,*) "The The rnew (q1)"
      WRITE(6,10) (R(I), I = 1, TOTREDNCO)
   10 Format (5(1X,F10.6))
      Write(6,*)
C
      IOFF4q0  = 3*NATOMS + 3*TOTREDNCO
      IOFF4X0  = IOFF4q0  + TOTREDNCO
C
      CALL DAXPY(TOTREDNCO, -1.0D0, SCRATCH_LCL, 1, R, 1)
      CALL DCOPY(TOTREDNCO, R, 1, SCRATCH_LCL(2*TOTREDNCO + 1), 1)
      CALL DCOPY(TOTREDNCO, SCRATCH_LCL, 1, SCRATCH_LCL(IOFF4q0 + 1),
     &           1)
      CALL DCOPY(3*NATOMS, Q, 1, SCRATCH_LCL(IOFF4X0 + 1), 1)

      CALL GETREC(20,'JOBARC','NUMREDCO', 1, NULLEVAL)
      CALL GETREC(20,'JOBARC','REDEVECS',NULLEVAL*TOTREDNCO*IINTFP,
     &            TMP1)
      CALL XGEMM('T', 'N', NULLEVAL, 1, TOTREDNCO, 1.0D0,
     &               TMP1, TOTREDNCO, SCRATCH_LCL, TOTREDNCO,
     &               0.0D0, TMP2, NULLEVAL)
      Write(6,"(a,I4)") "  Check for the redundencies in delta q0:",
     &                     NULLEVAL
      IF (NULLEVAL .NE. 0) THEN
      Write(6,"(6F13.7)") (tmp2(i), i=1, NULLEVAL)
      ENDIF
C
C
C The starting structure must be symmetric under the operations of
C the point group that was found for it; if it is not (the operations
C records do not belong to this geometry) do not symmetrize.
C
      CALL SYMMETRIZE_XYZ(Q, NATOMS, SYMERR, .FALSE.)
      DO_SYM = (SYMERR .LE. SYMTHRS)
      IF (.NOT. DO_SYM) Write(6,"(a,a,ES10.3,a)") "  The starting ",
     &   "geometry is not symmetric (", SYMERR, " bohr); the step is "
     &   //"not symmetrized."
C
      CONVERGED = .FALSE.
        CHANGED = .FALSE.
        NO_ITER = .FALSE.
          NITER = 0
      RMSBEST   = -1.0D0
C
C Bond angles can not go beyond 180 degrees. If the step asks for that,
C do not iterate and take the linear update.
C
      CALL GETREC(20, 'JOBARC', 'TNUMOBND', 1, TOTNOFBND)
      CALL GETREC(20, 'JOBARC', 'TNUMOANG', 1, TOTNOFANG)
      CALL CHECK4_180(SCRATCH_LCL(2*TOTREDNCO + 1), SCRATCH_LCL,
     &                TOTNOFBND, TOTNOFANG, TOTREDNCO, NO_ITER)
C
      DO WHILE (.NOT. CONVERGED .AND. NITER .LE. MAXITER .AND.
     &          .NOT. CHANGED .AND. .NOT. NO_ITER)
C
         CALL DCOPY(3*NATOMS, Q, 1, SCRATCH_LCL(3*TOTREDNCO + 1), 1)
         Write(6,*)
         WRITE(6,*) "The starting Cartesian during iteration; x_K"
         CALL OUTPUT(Q, 1, 3, 1, Natoms,3, Natoms, 1)
         Write(6,*)
         Write(6,*) "The starting internal increments:delta q_k"
         Write(6,*)
         WRITE(6,10) (SCRATCH_LCL(I), I=1, TOTREDNCO)
         Write(6,*)
         CALL XGEMM('N', 'N', 3*NATOMS, 1, TOTREDNCO, 1.0D0,
     &               AMATRX, 3*NATOMS, SCRATCH_LCL, TOTREDNCO,
     &               1.0D0, Q, 3*NATOMS)
         IF (DO_SYM) CALL SYMMETRIZE_XYZ(Q, NATOMS, SYMERR, .TRUE.)
C
         Print*, "The Cartesian after the update: x_k+1"
         CALL OUTPUT(Q, 1, 3, 1, Natoms,3, Natoms, 1)
         CALL GEN_NEW_RIC(Q, REDUNCO, IATNUM, NATOMS, TOTNOFBND,
     &                    NW_TOTREDNCO, TOTNOFANG, .FALSE.)
C
C In rare cases, it is possible that the # of RIC coordinates
C can change. If that happens leave the iterative procedure and
C and use the linear update.
C
         IF (NW_TOTREDNCO .GT. TOTREDNCO) CHANGED = .TRUE.
C
C
         CALL DCOPY(TOTREDNCO, REDUNCO, 1, R, 1)
         CALL DAXPY(TOTREDNCO, -1.0D0, SCRATCH_LCL(2*TOTREDNCO + 1),
     &              1, R, 1)
         CALL DCOPY(TOTREDNCO, SCRATCH_LCL(IOFF4q0 + 1), 1,
     &              SCRATCH_LCL, 1)
         CALL DAXPY(TOTREDNCO, -1.0D0, R, 1, SCRATCH_LCL, 1)
C
      WRITE(6,*) "The delta delta q_k"
      Write(6,*)
      WRITE(6,10) (SCRATCH_LCL(I), I=1, TOTREDNCO)
      Write(6,*)
         CALL REMOVE_360(SCRATCH_LCL, TOTNOFBND, TOTNOFANG,
     &                   TOTREDNCO)
C
         DO IREDUNCO = 1, TOTREDNCO
            IF (DABS(SCRATCH_LCL(IREDUNCO)) .LT. 1.0D-10) SCRATCH_LCL
     &          (IREDUNCO) = 0.0D0
         ENDDo
C
C The residual dq0 - (q(x_k+1) - q0) measures how well x_k+1
C reproduces the requested step. Keep the best iterate; it is used
C when the iterations do not converge (the first iterate is the
C linear update).
C
         CALL VSTAT(SCRATCH_LCL, STATS, TOTREDNCO)
         IF (RMSBEST .LT. 0.0D0 .OR. STATS(5) .LT. RMSBEST) THEN
            RMSBEST = STATS(5)
            CALL DCOPY(3*NATOMS, Q, 1, XBEST, 1)
         ENDIF
C
      CALL XGEMM('T', 'N', NULLEVAL, 1, TOTREDNCO, 1.0D0,
     &               TMP1, TOTREDNCO, SCRATCH_LCL, TOTREDNCO,
     &               0.0D0, TMP2, NULLEVAL)
      Write(6,*)
      Write(6,"(a,a,I4)") "  Check for the redundencies in delta",
     &                    " delta q_K:",NULLEVAL
      IF (NULLEVAL .NE. 0) THEN
      Write(6,"(6F13.7)") (tmp2(i), i=1, NULLEVAL)
      ENDIF
C
         CALL DAXPY(3*NATOMS, -1.0D0, Q, 1, SCRATCH_LCL
     &             (3*TOTREDNCO + 1), 1)
         CALL VSTAT(SCRATCH_LCL(3*TOTREDNCO+1), STATS, 3*NATOMS)
C
         IF (STATS(5)  .LE. EPSILON) CONVERGED = .TRUE.
         NITER = NITER + 1
C
      ENDDO
      IF (CONVERGED) Write(6, "(a,a,I2,a)") "  Iterative update ",
     &                         "converged in ", NITER, " iterations."
C
C If the iterative procedure fail to converge in maximum iterations
C alowed, take the best iterate (smallest residual). If no iterations
C were done, do a simple linear update.
C
      IF (.NOT. CONVERGED) THEN
         IF (RMSBEST .GE. 0.0D0 .AND. .NOT. CHANGED) THEN
            Write(6, "(a,a,ES10.3,a)") "  Iterative update did not ",
     &            "converge; best iterate (RMS residual ", RMSBEST,
     &            ") is used."
            CALL DCOPY(3*NATOMS, XBEST, 1, Q, 1)
         ELSE
            Write(6, "(a,a)") "  Iterative update did not converged.",
     &                        " linear update is used. "
            CALL DCOPY(3*NATOMS, SCRATCH_LCL(IOFF4X0 +1), 1, Q, 1)
            CALL XGEMM('N', 'N', 3*NATOMS, 1, TOTREDNCO, 1.0D0,
     &                 AMATRX, 3*NATOMS, SCRATCH_LCL(IOFF4q0 + 1),
     &                 TOTREDNCO, 1.0D0, Q, 3*NATOMS)
         ENDIF
      ENDIF
C
C Lets screen the Cartesian coordiantes and filter out changes
C beyond 7th decimal. Rounding does not preserve symmetry for
C operations other than reflections and C2 rotations about the
C axes, so symmetrize again afterwards.
C
      DO IATOM = 1, 3*NATOMS
         WRITE(HOLDQ, "(F12.7)") Q(IATOM)
         READ(HOLDQ,  "(F12.7)") Q(IATOM)
      ENDDO
      IF (DO_SYM) CALL SYMMETRIZE_XYZ(Q, NATOMS, SYMERR, .TRUE.)

      Write(6,*)
      WRITE(6,*) "The new Cartesian Coord., x+delta x0"
      WRITE(6,"(3F13.9)") (Q(I), I = 1, 3*NATOMS)
C
      CALL GEN_NEW_RIC(Q, REDUNCO, IATNUM, NATOMS, TOTNOFBND,
     &                 NW_TOTREDNCO, TOTNOFANG, .TRUE.)
C
      DO I = 1, TOTREDNCO
         STEP(I) = REDUNCO(I) - SCRATCH_LCL(2*TOTREDNCO + I)
      ENDDO
      CALL REMOVE_360(STEP, TOTNOFBND, TOTNOFANG, TOTREDNCO)
      DIFMAX = 0.0D0
      DO I = 1, TOTREDNCO
         DIFMAX = MAX(DIFMAX, DABS(STEP(I) -
     &                SCRATCH_LCL(IOFF4q0 + I)))
      ENDDO
      Write(6,"(a,ES10.3)") "  Max |actual - requested| RIC step: ",
     &                      DIFMAX

C
C Store the Cartesian coordinates with a different record label
C from COORD since finite difference calculations use COORD
C for the grid points, not the current point.
C
      Write(6,*)
      WRITE(6,*) "The new Cartesian Coord"
      WRITE(6,"(3F13.9)") (Q(I), I = 1, 3*NATOMS)

      CALL PUTREC(1,'JOBARC','COORD_OP',IINTFP*3*NATOMS,Q)
C
      RETURN
      END
