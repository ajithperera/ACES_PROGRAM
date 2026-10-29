










      SUBROUTINE RSPNS_NO(DAO, SCR, MAXCOR, NAO, NATOM, JPERT, IUHF)
C
C RESPONSE NATURAL ORBITAL (RNO) DECOMPOSITION OF THE FERMI-CONTACT
C RESPONSE. DAO IS THE AO-BASIS FIRST-ORDER SPIN DENSITY OF PERTURBATION
C JPERT (NUCLEUS B). ITS SYMMETRIC PART IS DIAGONALIZED IN THE S METRIC,
C
C    D = SUM(K) N(K) |PHI(K)><PHI(K)|,
C
C WHICH IS UNIQUE (NO LOCALIZATION OR PARTITIONING). FOR THE FC OPERATOR
C AT NUCLEUS A, H(A) = DELTA(R-R(A)), SO
C
C    X(A,B) = TR(D H(A)) = SUM(K) N(K) |PHI(K)(R(A))|**2,
C
C AND EACH MODE CONTRIBUTES ITS WEIGHT TIMES ITS DENSITY AT NUCLEUS A.
C THE CONTRIBUTIONS ARE PRINTED IN RAW UNITS AND AS PARTS OF THE REDUCED
C COUPLING K (10**19 T**2/J, BEFORE THE A<->B SYMMETRIZATION), AND THEIR
C SUM IS CHECKED AGAINST TR(D H(A)).
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MXPRT = 12, PLANCK = 6.62607015D-34)
      PARAMETER (FACTFC = -1.058311758D-14)
C


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



C
      CHARACTER*8 LABEL8
      DIMENSION DAO(NAO, NAO), SCR(MAXCOR), IORD(1000), ICENT(1000),
     &          IANG(1000), FRAC(3, 20)
C
      IF (NAO .GT. 1000) RETURN
      N2 = NAO*NAO
      IS    = 1
      IU    = IS   + N2
      IEV   = IU   + N2
      ISH   = IEV  + NAO
      ISMH  = ISH  + N2
      ID    = ISMH + N2
      IC    = ID   + N2
      ITMP  = IC   + N2
      IH    = ITMP + N2
      ICON  = IH   + N2
      IWRK  = ICON + NAO*NATOM
      LWRK  = 4*NAO
      ITOP  = IWRK + LWRK
      IF (ITOP .GT. MAXCOR) THEN
         WRITE(6,"(a)") " @RSPNS_NO: not enough memory, skipped."
         RETURN
      ENDIF
C
C K PER RAW UNIT: X -> J USES FACTFC*G(A)*G(B) AFTER THE FACTOR OF TWO
C APPLIED TO THE SYMMETRIZED TENSOR IN PERT_DENS_MAIN.
C
      PI    = 4.0D0*ATAN(1.0D0)
      XTOK  = 2.0D0*FACTFC*4.0D0*PI*PI/PLANCK*1.0D-19
C
C S**(1/2) AND S**(-1/2)
C
      CALL GETREC(20, 'JOBARC', 'AOOVRLAP', N2*IINTFP, SCR(IS))
      CALL DCOPY(N2, SCR(IS), 1, SCR(IU), 1)
      CALL DSYEV('V', 'U', NAO, SCR(IU), NAO, SCR(IEV), SCR(IWRK),
     &           LWRK, INFO)
      IF (INFO .NE. 0) RETURN
      CALL ZERO(SCR(ISH), N2)
      CALL ZERO(SCR(ISMH), N2)
      DO K = 1, NAO
         EK = SCR(IEV+K-1)
         DO J = 1, NAO
            UJK = SCR(IU+(K-1)*NAO+J-1)
            DO I = 1, NAO
               UIK = SCR(IU+(K-1)*NAO+I-1)
               IJ = (J-1)*NAO + I - 1
               SCR(ISH+IJ)  = SCR(ISH+IJ)  + UIK*UJK*SQRT(EK)
               SCR(ISMH+IJ) = SCR(ISMH+IJ) + UIK*UJK/SQRT(EK)
            ENDDO
         ENDDO
      ENDDO
C
C SYMMETRIC PART OF D, THEN S**(1/2) D S**(1/2) AND ITS EIGENVECTORS
C
      DO J = 1, NAO
         DO I = 1, NAO
            SCR(ID+(J-1)*NAO+I-1) = 0.5D0*(DAO(I,J) + DAO(J,I))
         ENDDO
      ENDDO
      CALL XGEMM('N', 'N', NAO, NAO, NAO, 1.0D0, SCR(ISH), NAO,
     &           SCR(ID), NAO, 0.0D0, SCR(ITMP), NAO)
      CALL XGEMM('N', 'N', NAO, NAO, NAO, 1.0D0, SCR(ITMP), NAO,
     &           SCR(ISH), NAO, 0.0D0, SCR(IU), NAO)
      CALL DSYEV('V', 'U', NAO, SCR(IU), NAO, SCR(IEV), SCR(IWRK),
     &           LWRK, INFO)
      IF (INFO .NE. 0) RETURN
C
C AO COEFFICIENTS OF THE RNOS: C = S**(-1/2) W
C
      CALL XGEMM('N', 'N', NAO, NAO, NAO, 1.0D0, SCR(ISMH), NAO,
     &           SCR(IU), NAO, 0.0D0, SCR(IC), NAO)
C
C CONTRIBUTION OF MODE K TO X(A,JPERT) FOR EVERY NUCLEUS A
C
      OPEN(UNIT=30, FILE='VPOUT', FORM='UNFORMATTED', STATUS='OLD')
      IRWND = 0
      WRITE(6,*)
      WRITE(6,"(a,I4)") " Response natural orbitals (FC), perturbed"//
     &                  " nucleus ", JPERT
      DO IA = 1, NATOM
         CALL SEEKLB('   DEN  ', IERR, IRWND, 30)
         IF (IERR .NE. 0) RETURN
         CALL LOADINT(SCR(IH), NATOM, N2, NAO, IUHF)
         IRWND = 1
         CALL XGEMM('N', 'N', NAO, NAO, NAO, 1.0D0, SCR(IH), NAO,
     &              SCR(IC), NAO, 0.0D0, SCR(ITMP), NAO)
         TOT = 0.0D0
         TRDH = DDOT(N2, SCR(ID), 1, SCR(IH), 1)
         DO K = 1, NAO
            Q = DDOT(NAO, SCR(IC+(K-1)*NAO), 1,
     &               SCR(ITMP+(K-1)*NAO), 1)
            SCR(ICON+(IA-1)*NAO+K-1) = SCR(IEV+K-1)*Q
            TOT = TOT + SCR(IEV+K-1)*Q
         ENDDO
         WRITE(6,"(a,I3,a,ES16.8,a,ES16.8,a,F10.5)")
     &        "   nucleus ", IA, ": sum over modes ", TOT,
     &        "  Tr(D h) ", TRDH, "  K ", TOT*XTOK
      ENDDO
C
C ORDER THE MODES BY THE LARGEST |CONTRIBUTION| TO ANY OTHER NUCLEUS
C
      DO K = 1, NAO
         IORD(K) = K
      ENDDO
      DO K = 1, NAO
         SCR(ITMP+K-1) = 0.0D0
         DO IA = 1, NATOM
            IF (IA .NE. JPERT) SCR(ITMP+K-1) = MAX(SCR(ITMP+K-1),
     &         ABS(SCR(ICON+(IA-1)*NAO+K-1)))
         ENDDO
      ENDDO
      DO K = 1, NAO-1
         DO L = K+1, NAO
            IF (SCR(ITMP+IORD(L)-1) .GT. SCR(ITMP+IORD(K)-1)) THEN
               ITMPK = IORD(K)
               IORD(K) = IORD(L)
               IORD(L) = ITMPK
            ENDIF
         ENDDO
      ENDDO
C
C LABELS ONLY: MULLIKEN-STYLE FRACTION OF EACH MODE ON EACH ATOM, SPLIT
C INTO S, P AND D-OR-HIGHER FUNCTIONS (THE DECOMPOSITION ITSELF USES NO
C PARTITIONING).
C
      CALL GETREC(20, 'JOBARC', 'CENTERBF', NAO, ICENT)
      CALL GETREC(20, 'JOBARC', 'ANGMOMBF', NAO, IANG)
      CALL XGEMM('N', 'N', NAO, NAO, NAO, 1.0D0, SCR(IS), NAO,
     &           SCR(IC), NAO, 0.0D0, SCR(ITMP), NAO)
C
      WRITE(6,"(a)") "   mode      n(k)      contributions to K"//
     &               " at each nucleus (10**19 T**2/J)"
      DO KK = 1, MIN(MXPRT, NAO)
         K = IORD(KK)
         WRITE(6,"(3x,I4,F12.6,8ES12.4)") K, SCR(IEV+K-1),
     &        (SCR(ICON+(IA-1)*NAO+K-1)*XTOK, IA = 1, MIN(NATOM,8))
         CALL ZERO(FRAC, 60)
         DO MU = 1, NAO
            IA = MIN(ICENT(MU), 20)
            L  = MIN(IANG(MU), 2) + 1
            FRAC(L, IA) = FRAC(L, IA) + SCR(IC+(K-1)*NAO+MU-1)*
     &                    SCR(ITMP+(K-1)*NAO+MU-1)
         ENDDO
         WRITE(6,"(12x,a,8(I3,':',3F6.2))") "s/p/d on atom",
     &        (IA, (FRAC(L, IA), L = 1, 3), IA = 1, MIN(NATOM,8))
      ENDDO
C
C SAVE THE MODES (WEIGHTS AND AO COEFFICIENTS) FOR PLOTTING
C
      WRITE(LABEL8, "(A,I2)") "RNOEVL", JPERT
      CALL PUTREC(20, 'JOBARC', LABEL8, NAO*IINTFP, SCR(IEV))
      WRITE(LABEL8, "(A,I2)") "RNOVEC", JPERT
      CALL PUTREC(20, 'JOBARC', LABEL8, N2*IINTFP, SCR(IC))
C
      RETURN
      END
