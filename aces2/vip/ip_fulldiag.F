#include "flags.h"
      SUBROUTINE IP_FULLDIAG(ICORE, MAXCOR, IUHF, ISPIN, NIP, NROOT)
C
C IP_ROOT_FIND=FULL_DIAG: explicit IP-EOM H-bar build and full
C diagonalization of one symmetry/spin block. Every root is found,
C together with its left eigenvector (rows of R**-1, so <L_k|R_l> =
C delta_kl by construction) and the 1h pole strength
C
C     P_p(k) = R1_p(k) * L1_p(k)
C
C of each occupied orbital p (the Lambda-free spectroscopic factor;
C Sum_k P_p(k) = 1 for every p when all roots are real). Results go to
C the output and to ip_eom_fulldiag.dat. Small systems only: three
C NIP x NIP matrices are held in core.
C
C The right and left H-bar are built independently (HC_MULT with
C ISIDE=1 and 2 on unit vectors); max|H_L - H_R**T| and the residual
C of each left eigenvector against H_L are printed as checks.
C
C UHF reference only (the RHF IP vector is not spin-adapted here).
C
CEND
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR)
C
#include "machsp.com"
C
      N2    = NIP*NIP
      LWORK = 64*NIP
C
      I000 = 1
      I010 = I000 + IINTFP*N2
      I020 = I010 + IINTFP*N2
      I030 = I020 + IINTFP*N2
      I040 = I030 + IINTFP*N2
      I050 = I040 + IINTFP*NIP
      I060 = I050 + IINTFP*NIP
      I070 = I060 + IINTFP*NIP
      I080 = I070 + 2*NIP
      IEND = I080 + IINTFP*LWORK
C
      IF (IEND .GT. MAXCOR) THEN
         WRITE(6,"(a,2I12)") " @IP_FULLDIAG-F, insufficient memory:",
     &                       IEND, MAXCOR
         CALL ERREX
      ENDIF
C
      CALL IP_FULLDIAG_DO(ICORE(I000), ICORE(I010), ICORE(I020),
     &                    ICORE(I030), ICORE(I040), ICORE(I050),
     &                    ICORE(I060), ICORE(I070), ICORE(I080),
     &                    LWORK, ICORE(IEND), MAXCOR-IEND+1,
     &                    IUHF, ISPIN, NIP, NROOT)
C
      RETURN
      END
C
      SUBROUTINE IP_FULLDIAG_DO(HR, HL, RV, RI, WR, WI, V, IORD, WORK,
     &                          LWORK, ICORE, MAXCOR, IUHF, ISPIN,
     &                          NIP, NROOT)
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION HR, HL, RV, RI, WR, WI, V, WORK
      DOUBLE PRECISION EIGVAL, OSCSTR
      DOUBLE PRECISION ONE, ZILCH, FACT, DUM, X, PTOT, R1SQ
      DOUBLE PRECISION HDIFF, HMAX, RESMAX, RESID, EVALIMX, PMAX
      DOUBLE PRECISION DNRM2, SDOT
      LOGICAL CMPLX
      DIMENSION HR(NIP,NIP), HL(NIP,NIP), RV(NIP,NIP), RI(NIP,NIP)
      DIMENSION WR(NIP), WI(NIP), V(NIP), IORD(2*NIP), WORK(LWORK)
      DIMENSION ICORE(MAXCOR)
C
#include "machsp.com"
#include "sym.com"
C
      COMMON /SINFO/ NS(8), SIRREP
      COMMON /SLISTS/ LS1IN, LS1OUT, LS2IN(2,2), LS2OUT(2,2)
      COMMON /EXTINF3/ IROOT, LOCROOT, ITROOT
      COMMON /ROOTS/ EIGVAL(100,8,3), OSCSTR(100,8,3)
C
      DATA ONE, ZILCH /1.0D0, 0.0D0/
      DATA FACT /27.2113957D0/
C
      N2   = NIP*NIP
      NOCC = POP(SIRREP, ISPIN)
C
      WRITE(6,*)
      WRITE(6,"(a,I8)") " @IP_FULLDIAG-I, building IP-EOM H-bar, dim ",
     &                  NIP
C
C BUILD RIGHT AND LEFT H-BAR, ONE UNIT VECTOR AT A TIME
C
      DO K = 1, NIP
         CALL ZERO(V, NIP)
         V(K) = ONE
         CALL PUTS_HHP(V, NIP, ISPIN, IUHF, LS1IN, LS2IN)
         CALL HC_MULT(ICORE, MAXCOR, IUHF, 1, ISPIN)
         CALL LOADS_IPVEC(HR(1,K), NIP, ISPIN, IUHF, LS1OUT, LS2OUT,
     &                    .FALSE., ICORE, MAXCOR)
C
         CALL ZERO(V, NIP)
         V(K) = ONE
         CALL PUTS_HHP(V, NIP, ISPIN, IUHF, LS1IN, LS2IN)
         CALL HC_MULT(ICORE, MAXCOR, IUHF, 2, ISPIN)
         CALL LOADS_IPVEC(HL(1,K), NIP, ISPIN, IUHF, LS1OUT, LS2OUT,
     &                    .FALSE., ICORE, MAXCOR)
      ENDDO
C
      HDIFF = ZILCH
      HMAX  = ZILCH
      DO J = 1, NIP
         DO I = 1, NIP
            HDIFF = MAX(HDIFF, ABS(HL(I,J) - HR(J,I)))
            HMAX  = MAX(HMAX, ABS(HR(I,J)))
         ENDDO
      ENDDO
      WRITE(6,"(a,1PE10.3,a,1PE10.3)")
     &   " @IP_FULLDIAG-I, max|H_L - H_R**T| = ", HDIFF,
     &   "   (max|H_R| = ", HMAX, ")"
C
C DIAGONALIZE H_R (RIGHT EIGENVECTORS ONLY; LEFT = ROWS OF R**-1)
C
      CALL DCOPY(N2, HR, 1, RI, 1)
      CALL DGEEV('N', 'V', NIP, RI, NIP, WR, WI, DUM, 1, RV, NIP,
     &           WORK, LWORK, INFO)
      IF (INFO .NE. 0) THEN
         WRITE(6,*) ' @IP_FULLDIAG-F, DGEEV failed, info = ', INFO
         CALL ERREX
      ENDIF
C
C SORT BY REAL PART (STABLE, SO A COMPLEX PAIR STAYS RE/IM ADJACENT)
C
      DO K = 1, NIP
         IORD(K) = K
      ENDDO
      DO K = 2, NIP
         J = IORD(K)
         L = K - 1
 10      IF (L .GE. 1) THEN
            IF (WR(IORD(L)) .GT. WR(J)) THEN
               IORD(L+1) = IORD(L)
               L = L - 1
               GO TO 10
            ENDIF
         ENDIF
         IORD(L+1) = J
      ENDDO
      DO K = 1, NIP
         CALL DCOPY(NIP, RV(1,IORD(K)), 1, RI(1,K), 1)
         WORK(K)     = WR(IORD(K))
         WORK(NIP+K) = WI(IORD(K))
      ENDDO
      CALL DCOPY(N2, RI, 1, RV, 1)
      CALL DCOPY(NIP, WORK, 1, WR, 1)
      CALL DCOPY(NIP, WORK(NIP+1), 1, WI, 1)
C
      EVALIMX = ZILCH
      DO K = 1, NIP
         EVALIMX = MAX(EVALIMX, ABS(WI(K)))
      ENDDO
      WRITE(6,"(a,1PE10.3)") " @IP_FULLDIAG-I, max|Im(eigenvalue)| = ",
     &                       EVALIMX
C
C LEFT EIGENVECTORS: ROWS OF R**-1, STORED AS COLUMNS OF RI
C
      CALL DCOPY(N2, RV, 1, HR, 1)
      CALL DGETRF(NIP, NIP, HR, NIP, IORD(NIP+1), INFO)
      IF (INFO .NE. 0) THEN
         WRITE(6,*) ' @IP_FULLDIAG-F, DGETRF failed, info = ', INFO
         CALL ERREX
      ENDIF
      CALL DGETRI(NIP, HR, NIP, IORD(NIP+1), WORK, LWORK, INFO)
      IF (INFO .NE. 0) THEN
         WRITE(6,*) ' @IP_FULLDIAG-F, DGETRI failed, info = ', INFO
         CALL ERREX
      ENDIF
      DO J = 1, NIP
         DO I = 1, NIP
            RI(I,J) = HR(J,I)
         ENDDO
      ENDDO
C
C CHECK EACH (REAL) LEFT EIGENVECTOR AGAINST THE INDEPENDENT H_L
C
      RESMAX = ZILCH
      DO K = 1, NIP
         IF (WI(K) .EQ. ZILCH) THEN
            CALL DGEMV('N', NIP, NIP, ONE, HL, NIP, RI(1,K), 1, ZILCH,
     &                 V, 1)
            CALL DAXPY(NIP, -WR(K), RI(1,K), 1, V, 1)
            RESID  = DNRM2(NIP, V, 1)/DNRM2(NIP, RI(1,K), 1)
            RESMAX = MAX(RESMAX, RESID)
         ENDIF
      ENDDO
      WRITE(6,"(a,a,1PE10.3)") " @IP_FULLDIAG-I, max left-eigenvector",
     &   " residual |H_L*l - w*l|/|l| = ", RESMAX
C
C SUMMARY: ENERGIES AND 1H POLE STRENGTHS
C
      OPEN(UNIT=84, FILE='ip_eom_fulldiag.dat', FORM='FORMATTED',
     &     STATUS='UNKNOWN')
      WRITE(84,"(a)") "# IP_ROOT_FIND=FULL_DIAG: NIP NOCC IRREP ISPIN"
      WRITE(84,"(a,4I8)") "# ", NIP, NOCC, SIRREP, ISPIN
      WRITE(84,"(2a)") "# root  IP(Ha)  IP(eV)  Im(Ha)  P_1..P_nocc",
     &                 "  (P_p = R1_p*L1_p; 0 for complex roots)"
C
      WRITE(6,*)
      WRITE(6,"(a)") "  IP-EOM FULL DIAGONALIZATION: ALL ROOTS"
      WRITE(6,"(2a)") "  P_p = R1_p*L1_p (1h pole strength of occ. p),",
     &                " |R1|^2 with |R| = 1"
      WRITE(6,"(2a)") "   root     IP(eV)     Im(eV)    Sum_p P_p",
     &                "   |R1|^2   max P_p  (p)"
C
      DO K = 1, NIP
         CMPLX = WI(K) .NE. ZILCH
         PTOT = ZILCH
         PMAX = ZILCH
         IPMAX = 0
         DO IP = 1, NOCC
            IF (CMPLX) THEN
               WORK(IP) = ZILCH
            ELSE
               WORK(IP) = RV(IP,K)*RI(IP,K)
            ENDIF
            PTOT = PTOT + WORK(IP)
            IF (ABS(WORK(IP)) .GT. ABS(PMAX)) THEN
               PMAX  = WORK(IP)
               IPMAX = IP
            ENDIF
         ENDDO
         R1SQ = SDOT(NOCC, RV(1,K), 1, RV(1,K), 1)
         IF (CMPLX) THEN
            WRITE(6,"(I7,2F11.5,a)") K, WR(K)*FACT, WI(K)*FACT,
     &            "   complex root, pole strength not computed"
         ELSE
            WRITE(6,"(I7,2F11.5,3F10.5,I5)") K, WR(K)*FACT,
     &            WI(K)*FACT, PTOT, R1SQ, PMAX, IPMAX
         ENDIF
         WRITE(84,"(I6,F16.10,F14.6,1PE12.3,0P,100F12.8)") K, WR(K),
     &         WR(K)*FACT, WI(K), (WORK(IP), IP=1,NOCC)
      ENDDO
      CLOSE(UNIT=84)
C
C SUM RULE Sum_k P_p(k) = 1 (EXACT FOR AN ALL-REAL SPECTRUM)
C
      WRITE(6,*)
      WRITE(6,"(a)") "  Sum over roots of P_p (should be 1):"
      DO IP = 1, NOCC
         X = ZILCH
         DO K = 1, NIP
            IF (WI(K) .EQ. ZILCH) X = X + RV(IP,K)*RI(IP,K)
         ENDDO
         WRITE(6,"(a,I4,F14.10)") "   p = ", IP, X
      ENDDO
      WRITE(6,*)
C
C LOWEST NROOT ROOTS GO TO THE USUAL IP SUMMARY TABLE
C
      IROOT = 0
      DO K = 1, NIP
         IF (IROOT .LT. MIN(NROOT,100)) THEN
            IROOT = IROOT + 1
            EIGVAL(IROOT, SIRREP, ISPIN) = WR(K)
         ENDIF
      ENDDO
      RETURN
      END
