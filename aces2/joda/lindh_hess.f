










      SUBROUTINE LINDH_HESS(HI, NRIC)
C
C Diagonal model Hessian for redundant internal coordinates after
C R. Lindh, A. Bernhardsson, G. Karlstrom, P.-A. Malmqvist, Chem.
C Phys. Lett. 241, 423 (1995):
C
C    rho_ij  = exp(alpha_ij (r_ref,ij^2 - r_ij^2))
C    bond    k = 0.45  rho_ij
C    angle   k = 0.15  rho_ij rho_jk
C    dihedral k = 0.005 rho_ij rho_jk rho_kl
C
C alpha and r_ref (bohr) depend on the rows of the periodic table of
C the two atoms (rows beyond the third use the third-row values).
C A force constant is never smaller than FLOOR so that long (H-bond,
C fragment) coordinates do not make the guess nearly singular.
C 10/2026.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
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
      INTEGER TOTNOFBND, TOTNOFANG, TOTNOFDIH
      PARAMETER (FLOOR = 1.0D-3)
      DIMENSION HI(NRIC,NRIC), IREDUNCO(4,MAXREDUNCO), ALPHA(3,3),
     &          RREF(3,3)
      DATA ALPHA /1.0000D0, 0.3949D0, 0.3949D0,
     &            0.3949D0, 0.2800D0, 0.2800D0,
     &            0.3949D0, 0.2800D0, 0.2800D0/
      DATA RREF  /1.35D0, 2.10D0, 2.53D0,
     &            2.10D0, 2.87D0, 3.40D0,
     &            2.53D0, 3.40D0, 3.40D0/
C
      CALL GETREC(20, 'JOBARC', 'CONTEVIT', 4*NRIC, IREDUNCO)
      CALL GETREC(20, 'JOBARC', 'TNUMOBND', 1, TOTNOFBND)
      CALL GETREC(20, 'JOBARC', 'TNUMOANG', 1, TOTNOFANG)
      CALL GETREC(20, 'JOBARC', 'TNUMODIH', 1, TOTNOFDIH)
C
      CALL ZERO(HI, NRIC*NRIC)
      DO I = 1, NRIC
         IF (I .LE. TOTNOFBND) THEN
            F = 0.45D0*RHO(IREDUNCO(1,I), IREDUNCO(2,I))
         ELSE IF (I .LE. TOTNOFBND + TOTNOFANG) THEN
            F = 0.15D0*RHO(IREDUNCO(1,I), IREDUNCO(2,I))*
     &                 RHO(IREDUNCO(2,I), IREDUNCO(3,I))
         ELSE
            F = 0.005D0*RHO(IREDUNCO(1,I), IREDUNCO(2,I))*
     &                  RHO(IREDUNCO(2,I), IREDUNCO(3,I))*
     &                  RHO(IREDUNCO(3,I), IREDUNCO(4,I))
         ENDIF
         HI(I,I) = MAX(F, FLOOR)
      ENDDO
C
      Write(6,"(a)") " @LINDH_HESS: model Hessian diagonal"
      Write(6,"(6(1x,F10.6))") (HI(I,I), I=1, NRIC)
      RETURN
C
      CONTAINS
C
      DOUBLE PRECISION FUNCTION RHO(IA, IB)
      INTEGER IA, IB, IRA, IRB
      DOUBLE PRECISION R2
      IRA = IROW(IATNUM(IA))
      IRB = IROW(IATNUM(IB))
      R2  = (Q(3*IA-2)-Q(3*IB-2))**2 + (Q(3*IA-1)-Q(3*IB-1))**2 +
     &      (Q(3*IA  )-Q(3*IB  ))**2
      RHO = DEXP(ALPHA(IRA,IRB)*(RREF(IRA,IRB)**2 - R2))
      END FUNCTION RHO
C
      INTEGER FUNCTION IROW(IZ)
      INTEGER IZ
      IROW = 3
      IF (IZ .LE. 10) IROW = 2
      IF (IZ .LE. 2)  IROW = 1
      END FUNCTION IROW
C
      END
