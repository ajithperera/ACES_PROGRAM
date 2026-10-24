      SUBROUTINE EVAL_QSTLST_PATH_RIC(QST_TANGENT, LST_TANGENT, WORK,
     &                                PMAT, LST, QST)
C
C The normalized tangent to the synchronous transit path at the current
C point X in redundant internal coordinates (C. Peng and H. B. Schlegel,
C Israel J. Chem. 33, 449 (1993)):
C
C   LST (eq. 1):  T = (P - R)/|P - R|
C   QST (eq. 8):  T ~ (P - X)/|P - X|^2 - (R - X)/|R - X|^2,
C                 the tangent at X of the circle through R, X and P.
C
C R and P are the RIC values of the reactant and product (RX_RICS_,
C PR_RICS_, computed with the coordinate set of the current structure
C in GEN_REDNINTRNLS), X the current values (R in /COORD/). Dihedral
C differences are brought to (-pi, pi]. The tangent is projected onto
C the nonredundant space with PMAT (= G^-G, Peng et al., J. Comput.
C Chem. 17, 49 (1996)) and normalized. When the tangent is not defined
C (X on top of R or P) it is returned as zero and the caller does not
C use it. 10/2026.
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
#include "mxatms.par"
#include "machsp.com"
#include "coord.com"
C
      LOGICAL QST, LST
      DOUBLE PRECISION LST_TANGENT
      INTEGER TOTREDNCO, TOTNOFBND, TOTNOFANG
C
      COMMON /USINT/NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
C
      DIMENSION WORK(NX*NX), LST_TANGENT(*), QST_TANGENT(*),
     &          PMAT(NX*NX)
C
      PARAMETER (TINY = 1.0D-10)
C
      CALL GETREC(20, 'JOBARC', 'REDNCORD', 1, TOTREDNCO)
      CALL GETREC(20, 'JOBARC', 'TNUMOBND', 1, TOTNOFBND)
      CALL GETREC(20, 'JOBARC', 'TNUMOANG', 1, TOTNOFANG)
C
      IREACT = 1
      IPROD  = IREACT + TOTREDNCO
      ITAN   = IPROD  + TOTREDNCO
      IEND   = ITAN   + TOTREDNCO
      IF (IEND .GT. NX*NX) CALL INSMEM("@EVAL_QSTLST_PATH_RIC",
     &                                  IEND, NX*NX)
C
      CALL GETREC(20, 'JOBARC', 'RX_RICS_', TOTREDNCO*IINTFP,
     &            WORK(IREACT))
      CALL GETREC(20, 'JOBARC', 'PR_RICS_', TOTREDNCO*IINTFP,
     &            WORK(IPROD))
C
      IF (LST) THEN
C
C P - R
C
         DO I = 1, TOTREDNCO
            WORK(ITAN-1+I) = WORK(IPROD-1+I) - WORK(IREACT-1+I)
         ENDDO
         CALL REMOVE_360(WORK(ITAN), TOTNOFBND, TOTNOFANG, TOTREDNCO)
C
      ELSE IF (QST) THEN
C
C R - X and P - X
C
         DO I = 1, TOTREDNCO
            WORK(IREACT-1+I) = WORK(IREACT-1+I) - R(I)
            WORK(IPROD-1+I)  = WORK(IPROD-1+I)  - R(I)
         ENDDO
         CALL REMOVE_360(WORK(IREACT), TOTNOFBND, TOTNOFANG, TOTREDNCO)
         CALL REMOVE_360(WORK(IPROD),  TOTNOFBND, TOTNOFANG, TOTREDNCO)
C
         DRX2 = DDOT(TOTREDNCO, WORK(IREACT), 1, WORK(IREACT), 1)
         DPX2 = DDOT(TOTREDNCO, WORK(IPROD),  1, WORK(IPROD),  1)
         CALL ZERO(WORK(ITAN), TOTREDNCO)
         IF (DRX2 .GT. TINY .AND. DPX2 .GT. TINY) THEN
            CALL DAXPY(TOTREDNCO,  1.0D0/DPX2, WORK(IPROD),  1,
     &                 WORK(ITAN), 1)
            CALL DAXPY(TOTREDNCO, -1.0D0/DRX2, WORK(IREACT), 1,
     &                 WORK(ITAN), 1)
         ENDIF
      ENDIF
C
C Project onto the nonredundant space and normalize.
C
      CALL XGEMM("N", "N", TOTREDNCO, 1, TOTREDNCO, 1.0D0,
     &           PMAT, TOTREDNCO, WORK(ITAN), TOTREDNCO,
     &           0.0D0, WORK(IREACT), TOTREDNCO)
      TNORM = DSQRT(DDOT(TOTREDNCO, WORK(IREACT), 1, WORK(IREACT), 1))
      IF (TNORM .GT. TINY) THEN
         CALL DSCAL(TOTREDNCO, 1.0D0/TNORM, WORK(IREACT), 1)
      ELSE
         CALL ZERO(WORK(IREACT), TOTREDNCO)
         WRITE(6,"(a)") "  The LST/QST tangent is not defined at this"//
     &                  " point; it is not used."
      ENDIF
C
      IF (LST) THEN
         CALL DCOPY(TOTREDNCO, WORK(IREACT), 1, LST_TANGENT, 1)
      ELSE
         CALL DCOPY(TOTREDNCO, WORK(IREACT), 1, QST_TANGENT, 1)
      ENDIF
C
#ifdef _DEBUG_LVL0
      Write(6,"(a)") " @EVAL_QSTLST_PATH_RIC: the projected tangent"
      Write(6,"(6(1x,F10.6))") (WORK(IREACT-1+I), I=1, TOTREDNCO)
#endif
      RETURN
      END
