










      SUBROUTINE CNTINN(TWOIN,TWOOUT,ISTRCF,IOFFVC,WORK,
     &                  NSET3,NSET4,NPRIM3,NPRIM4,NCONT3,NCONT4,
     &                  NINT12,TPRI34,TCON34)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PURPOSE: TRANSFORMATION OF TWO INNERMOST INDICES
C              THIRD INDEX IS NEXT INNTERMOST
C              FOURTH INDEX IS INNERMOST
C
C     NOTE: AFTER TRANSFORMATION ORDER OF TWO INNERMOST INDICES
C           IS REVERSED
C
C     IN:  TWOIN(NINT12*NPRIM3*NPRIM4)
C          COEF3(NCONT3*NPRIM3)
C          COEF4(NCONT4*NPRIM4)
C
C     OUT: TWOOUT(NINT12*NCONT4*NCONT3)
C
C     SCRATCH: WORK(NINT12*NPRIM3*NCONT4)
C

c These parameters are gathered from vmol and vdint and are used by ecp
c as well. It just so happens that the vmol parameters do not exist in
c vdint and vice versa. LET'S TRY TO KEEP IT THAT WAY!

c VMOL PARAMETERS ------------------------------------------------------

C     MAXPRIM - Maximum number of primitives for a given shell.
      INTEGER    MAXPRIM
      PARAMETER (MAXPRIM=72)

C     MAXFNC  - Maximum number of contracted functions for a given shell.
C               (vmol/readin requires this to be the same as MAXPRIM)
      INTEGER    MAXFNC
      PARAMETER (MAXFNC=MAXPRIM)

C     NHT     - Maximum angular momentum
      INTEGER    NHT
      PARAMETER (NHT=7)

C     MAXATM  - Maximum number of atoms
      INTEGER    MAXATM
      PARAMETER (MAXATM=100)

C     MXTNPR  - Maximum total number of primitives for all symmetry
C               inequivalent centers.
      INTEGER    MXTNPR
      PARAMETER (MXTNPR=MAXPRIM*MAXPRIM)

C     MXTNCC  - Maximum total number of contraction coefficients for
C               all symmetry inequivalent centers.
      INTEGER    MXTNCC
      PARAMETER (MXTNCC=180000)

C     MXTNSH  - Maximum total number of shells for all symmetry
C               inequivalent centers.
      INTEGER    MXTNSH
      PARAMETER (MXTNSH=200)

C     MXCBF   - Maximum number of Cartesian basis functions for the
C               whole system (NOT the number of contracted functions).
c mxcbf.par : begin

c MXCBF := the maximum number of Cartesian basis functions (limited by vmol)

c This parameter is the same as MAXBASFN. Do NOT change this without changing
c maxbasfn.par as well.

      INTEGER MXCBF
      PARAMETER (MXCBF=1000)
c mxcbf.par : end

c VDINT PARAMETERS -----------------------------------------------------

C     MXPRIM - Maximum number of primitives for all symmetry
C              inequivalent centers.
      INTEGER    MXPRIM
      PARAMETER (MXPRIM=MXTNPR)

C     MXSHEL - Maximum number of shells for all symmetry inequivalent centers.
      INTEGER    MXSHEL
      PARAMETER (MXSHEL=MXTNSH)

C     MXCORB - Maximum number of contracted basis functions.
      INTEGER    MXCORB
      PARAMETER (MXCORB=MXCBF)

C     MXORBT - Length of the upper or lower triangle length of MXCORB.
      INTEGER    MXORBT
      PARAMETER (MXORBT=MXCORB*(MXCORB+1)/2)

C     MXAOVC - Maximum number of subshells per center.
      INTEGER    MXAOVC,    MXAOSQ
      PARAMETER (MXAOVC=32, MXAOSQ=MXAOVC*MXAOVC)

c     MXCONT - ???
      INTEGER    MXCONT
      PARAMETER (MXCONT=MXAOVC)

      COMMON /CCFCOM/ CONT1 (MXCONT*MXAOVC), CONT2 (MXCONT*MXAOVC),
     &                CONT3 (MXCONT*MXAOVC), CONT4 (MXCONT*MXAOVC),
     &                CONTT1(MXCONT*MXAOVC), CONTT2(MXCONT*MXAOVC),
     &                CONTT3(MXCONT*MXAOVC), CONTT4(MXCONT*MXAOVC),
     &                NUC1X (MXAOVC),        NUC2X (MXAOVC),
     &                NUC3X (MXAOVC),        NUC4X (MXAOVC),
     &                NRC1X (MXAOVC),        NRC2X (MXAOVC),
     &                NRC3X (MXAOVC),        NRC4X (MXAOVC)
      LOGICAL TPRI34, TCON34
      DIMENSION TWOIN(*), TWOOUT(*), WORK(*)
C
      DATA ONE,AZERO /1.D0,0.D0/
C
C     TRANSFORM FOURTH INDEX
C
      IF (TPRI34) THEN
       N34D = NPRIM3*(NPRIM3 + 1)/2
       N34Q = NPRIM3*NPRIM4
       IADR = 1
       IOFF = NPRIM3*NINT12*NCONT4
       DO 100 I = 1, NPRIM3
        DO 110 J = 1, I
         IJ = IOFF + (I - 1)*NPRIM3 + J
         JI = IOFF + (J - 1)*NPRIM3 + I
         CALL SCOPY(NINT12,TWOIN(IADR),N34D,WORK(IJ),N34Q)
         CALL SCOPY(NINT12,TWOIN(IADR),N34D,WORK(JI),N34Q)
         IADR = IADR + 1
  110   CONTINUE
  100  CONTINUE
       IF (NSET4 .EQ. 1) THEN
        CALL XGEMM('N','N',NCONT4,NPRIM3*NINT12,
     &              NPRIM4,ONE,CONT3(ISTRCF),NCONT4,
     &              WORK(IOFF+1),NPRIM4,AZERO,
     &              WORK(1),NCONT4)
       ELSE
        ICONT  = ISTRCF
        IWORK1 = IOFF + 1
        IWORK2 = 1
        DO 200 I = 1, NSET4
         NRC4 = NRC3X(IOFFVC + I)
         NUC4 = NUC3X(IOFFVC + I)
         CALL XGEMM('N','N',NRC4,NPRIM3*NINT12,
     &               NUC4,ONE,CONT3(ICONT),NRC4,
     &               WORK(IWORK1),NPRIM4,AZERO,
     &               WORK(IWORK2),NCONT4)
         ICONT  = ICONT  + NRC4*NUC4
         IWORK1 = IWORK1 + NUC4
         IWORK2 = IWORK2 + NRC4
  200   CONTINUE
       END IF
      ELSE
       IF (NSET4 .EQ. 1) THEN
        CALL XGEMM('N','N',NCONT4,NPRIM3*NINT12,
     &              NPRIM4,ONE,CONT3(ISTRCF),NCONT4,
     &              TWOIN(1),NPRIM4,AZERO,
     &              WORK(1),NCONT4)
       ELSE
        ICONT  = ISTRCF
        IWORK1 = 1
        IWORK2 = 1
        DO 300 I = 1, NSET4
         NRC4 = NRC3X(IOFFVC + I)
         NUC4 = NUC3X(IOFFVC + I)
C
C 03/09/97. In the original code NUC4 argument was missing.
C This probably would not have caused problems because NSET4 
C is 1, AP. 
C
         CALL XGEMM('N','N',NRC4,NPRIM3*NINT12,
     &                   NUC4, ONE,CONT3(ICONT),NRC4,
     &                    TWOIN(IWORK1),NPRIM4,AZERO,
     &                    WORK(IWORK2),NCONT4)
         ICONT  = ICONT  + NRC4*NUC4
         IWORK1 = IWORK1 + NUC4
         IWORK2 = IWORK2 + NRC4
  300   CONTINUE
       END IF
      END IF
C
C     CHANGE ORDER OF THIRD AND FOURTH INDICES
C
      NPR34 = NPRIM3*NCONT4
      IADR10 = 0
      IADR20 = 0
      DO 400 I = 1, NPRIM3
       DO 410 J = 1, NCONT4
        IADR1 = IADR10 + J - 1
        IADR2 = IADR20 + (J - 1)*NPRIM3
        DO 420 K = 1, (NINT12 - 1)*NPR34 + 1, NPR34
         TWOIN(IADR2 + K) = WORK(IADR1 + K)
  420   CONTINUE
  410  CONTINUE
       IADR10 = IADR10 + NCONT4
       IADR20 = IADR20 + 1
  400 CONTINUE
C
C     TRANSFORM THIRD INDEX
C
      NCT412 = NCONT4*NINT12
      IF (TCON34) THEN
       IF (NSET3 .EQ. 1) THEN
        CALL XGEMM('N','N',NCONT3,NCT412,NPRIM3,ONE,CONT4(ISTRCF),
     &             NCONT3,TWOIN(1),NPRIM3,AZERO,WORK(1),NCONT3)
       ELSE
        ICONT  = ISTRCF
        IWORK1 = 1
        IWORK2 = 1
        DO 500 I = 1, NSET3
         NRC3 = NRC4X(IOFFVC + I)
         NUC3 = NUC4X(IOFFVC + I)
         CALL XGEMM('N','N',NRC3,NCT412,NUC3,ONE,CONT4(ICONT),NRC3,
     &              TWOIN(IWORK1),NPRIM3,AZERO,WORK(IWORK2),NCONT3)
         ICONT  = ICONT  + NRC3*NUC3
         IWORK1 = IWORK1 + NUC3
         IWORK2 = IWORK2 + NRC3
  500   CONTINUE
       END IF
       NCTSQ = NCONT3*NCONT3
       NCTTR = NCONT3*(NCONT3 + 1)/2
       IADR2 = 0
       DO 600 I = 1, NCONT3
        IADR1 = (I - 1)*NCONT3
        DO 610 J = 1, I
         CALL SCOPY(NINT12,WORK(IADR1+J),NCTSQ,TWOOUT(IADR2+J),
     &              NCTTR)
  610   CONTINUE
        IADR2 = IADR2 + I
  600  CONTINUE
      ELSE
       IF (NSET3 .EQ. 1) THEN
        CALL XGEMM('N','N',NCONT3,NCT412,NPRIM3,ONE,CONT4(ISTRCF),
     &             NCONT3,TWOIN(1),NPRIM3,AZERO,TWOOUT(1),NCONT3)
       ELSE
        ICONT  = ISTRCF
        IWORK1 = 1
        IWORK2 = 1
        DO 700 I = 1, NSET3
         NRC3 = NRC4X(IOFFVC + I)
         NUC3 = NUC4X(IOFFVC + I)
         CALL XGEMM('N','N',NRC3,NCT412,NUC3,ONE,CONT4(ICONT),NRC3,
     &              TWOIN(IWORK1),NPRIM3,AZERO,TWOOUT(IWORK2),NCONT3)
         ICONT  = ICONT  + NRC3*NUC3
         IWORK1 = IWORK1 + NUC3
         IWORK2 = IWORK2 + NRC3
  700   CONTINUE
       END IF
      END IF
      RETURN
      END
