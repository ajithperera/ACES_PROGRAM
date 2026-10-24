      SUBROUTINE LINDH_HESS(HI, NRIC)
C
C Diagonal model Hessian for redundant internal coordinates after
C R. Lindh, A. Bernhardsson, G. Karlstrom, P.-A. Malmqvist, Chem.
C Phys. Lett. 241, 423 (1995):
C
C    bond     k = 0.45  rho_ij
C    angle    k = 0.15  rho_ij rho_jk
C    dihedral k = 0.005 rho_ij rho_jk rho_kl     (rho: LINDH_RHO)
C The out-of-plane coordinates (the last TNUMOOPS dihedrals, improper
C torsions A-B-C-D about the center C, see ASSIGN_OOPS) are not in
C the model; they get 0.045 rho_ac rho_bc rho_cd (the value used for
C out-of-plane bends in geomeTRIC's guess).
C
C A force constant is never smaller than FLOOR so that long (H-bond,
C fragment) coordinates do not make the guess nearly singular.
C 10/2026.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
#include "mxatms.par"
#include "coord.com"
C
      INTEGER TOTNOFBND, TOTNOFANG, TOTNOFDIH
      PARAMETER (FLOOR = 1.0D-3)
      DIMENSION HI(NRIC,NRIC), IREDUNCO(4,MAXREDUNCO)
C
      CALL GETREC(20, 'JOBARC', 'CONTEVIT', 4*NRIC, IREDUNCO)
      CALL GETREC(20, 'JOBARC', 'TNUMOBND', 1, TOTNOFBND)
      CALL GETREC(20, 'JOBARC', 'TNUMOANG', 1, TOTNOFANG)
      CALL GETREC(20, 'JOBARC', 'TNUMODIH', 1, TOTNOFDIH)
C
      CALL ZERO(HI, NRIC*NRIC)
      DO I = 1, NRIC
         IF (I .LE. TOTNOFBND) THEN
            F = 0.45D0*RHOAB(IREDUNCO(1,I), IREDUNCO(2,I))
         ELSE IF (I .LE. TOTNOFBND + TOTNOFANG) THEN
            F = 0.15D0*RHOAB(IREDUNCO(1,I), IREDUNCO(2,I))*
     &                 RHOAB(IREDUNCO(2,I), IREDUNCO(3,I))
         ELSE IF (I .GT. NRIC - NOOPS) THEN
            F = 0.045D0*RHOAB(IREDUNCO(1,I), IREDUNCO(3,I))*
     &                  RHOAB(IREDUNCO(2,I), IREDUNCO(3,I))*
     &                  RHOAB(IREDUNCO(3,I), IREDUNCO(4,I))
         ELSE
            F = 0.005D0*RHOAB(IREDUNCO(1,I), IREDUNCO(2,I))*
     &                  RHOAB(IREDUNCO(2,I), IREDUNCO(3,I))*
     &                  RHOAB(IREDUNCO(3,I), IREDUNCO(4,I))
         ENDIF
         HI(I,I) = MAX(F, FLOOR)
      ENDDO
C
#ifdef _DEBUG_LVL0
      Write(6,"(a)") " @LINDH_HESS: model Hessian diagonal"
      Write(6,"(6(1x,F10.6))") (HI(I,I), I=1, NRIC)
#endif
      RETURN
C
      CONTAINS
C
      DOUBLE PRECISION FUNCTION RHOAB(IA, IB)
      INTEGER IA, IB
      DOUBLE PRECISION R2, LINDH_RHO
      R2  = (Q(3*IA-2)-Q(3*IB-2))**2 + (Q(3*IA-1)-Q(3*IB-1))**2 +
     &      (Q(3*IA  )-Q(3*IB  ))**2
      RHOAB = LINDH_RHO(IATNUM(IA), IATNUM(IB), R2)
      END FUNCTION RHOAB
C
      END
