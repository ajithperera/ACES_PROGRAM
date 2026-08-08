










      SUBROUTINE CNTOUT(TWOIN,TWOOUT,WORK1,WORK2,LWORK1,LWORK2,
     &                  ISTRCF,IOFFVC,NSET1,NSET2,
     &                  NPRIM1,NPRIM2,NCONT1,NCONT2,
     &                  NINT34,TPRI12,TCON12)
C
C
C     Purpose: Transformation of two outermost indices
C              Index 1 is outermost index
C              Index 2 is next outermost index
C
C     Note: After transformation the order of two outermost indices
C           is reversed
C
C     In:  TWOIN(NPRIM1*NPRIM2*NINT34)
C
C     Out: TWOOUT(NCONT2*NCONT1*NINT34)
C
C     Scratch: WORK1(NCONT1*NPRIM2*NINT34)
C
C              WORK2(NPRIM1*NPRIM2*NINT34)              (TPRI12)
C              WORK2(NCONT1*NPRIM2*NINT34)         (.NOT.TPRI12)
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

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
C
      LOGICAL TPRI12, TCON12
      DIMENSION TWOIN(NPRIM1*NPRIM2*NINT34),
     &          TWOOUT(NCONT1*NCONT2*NINT34),
     &          WORK1(LWORK1),WORK2(LWORK2)
C
      DATA AZERO,ONE /0.D0,1.D0/
C
      NPR234 = NPRIM2*NINT34
      NCT134 = NCONT1*NINT34
C
C TRANSFORM FIRST INDEX
C
      IF (TPRI12) THEN
       DO 100 I = 1, NPRIM1
        IJ = (I - 1)*NPRIM1*NINT34
        IADR = ((I-1)*I)/2*NINT34
        DO 110 K = 1, I*NINT34
          WORK2(IJ + K) = TWOIN(IADR + K)
  110    CONTINUE
  100  CONTINUE
       IADR = NINT34
       DO 200 I = 2, NPRIM1
        DO 210 J = 1, I - 1
         JI = ((J - 1)*NPRIM1 + I - 1)*NINT34
         DO 220 K = 1, NINT34
          WORK2(JI + K) = TWOIN(IADR + K)
  220    CONTINUE
         IADR = IADR + NINT34
  210   CONTINUE
        IADR = IADR + NINT34
  200  CONTINUE
       IWORK1 = 1
       IWORK2 = 1
       ICONT  = ISTRCF
       DO 300 I = 1, NSET1
        NUC1 = NUC2X(IOFFVC + I)
        NRC1 = NRC2X(IOFFVC + I)
        CALL XGEMM('N','N',NPR234,NRC1,NUC1,ONE,WORK2(IWORK2),NPR234,
     &             CONTT2(ICONT),NUC1,AZERO,WORK1(IWORK1),NPR234)
        IWORK1 = IWORK1 + NRC1*NPR234
        IWORK2 = IWORK2 + NUC1*NPR234
        ICONT  = ICONT  + NRC1*NUC1
  300  CONTINUE
      ELSE
       IWORK1 = 1
       IWORK2 = 1
       ICONT  = ISTRCF
       DO 400 I = 1, NSET1
        NUC1 = NUC2X(IOFFVC + I)
        NRC1 = NRC2X(IOFFVC + I)
        CALL XGEMM('N','N',NPR234,NRC1,NUC1,ONE,TWOIN(IWORK2),NPR234,
     &             CONTT2(ICONT),NUC1,AZERO,WORK1(IWORK1),NPR234)
        IWORK1 = IWORK1 + NRC1*NPR234
        IWORK2 = IWORK2 + NUC1*NPR234
        ICONT  = ICONT  + NRC1*NUC1
  400  CONTINUE
      END IF
C
C CHANGE ORDER OF FIRST AND SECOND INDICES
C
      IADR1  = 0
      IADR20 = 0
      DO 500 I = 1, NCONT1
       IADR2 = IADR20
       DO 510 J = 1, NPRIM2
        DO 520 K = 1, NINT34
         WORK2(IADR2+K) = WORK1(IADR1 + K)
  520   CONTINUE
        IADR1 = IADR1 + NINT34
        IADR2 = IADR2 + NCT134
  510  CONTINUE
       IADR20 = IADR20 + NINT34
  500 CONTINUE
C
C TRANSFORM SECOND INDEX
C
      IF (TCON12) THEN
       IWORK1 = 1
       IWORK2 = 1
       ICONT  = ISTRCF
       DO 600 I = 1, NSET2
        NUC2 = NUC1X(IOFFVC + I)
        NRC2 = NRC1X(IOFFVC + I)
        CALL XGEMM('N','N',NCT134,NRC2,NUC2,ONE,WORK2(IWORK2),NCT134,
     &             CONTT1(ICONT),NUC2,AZERO,WORK1(IWORK1),NCT134)
        IWORK1 = IWORK1 + NRC2*NCT134
        IWORK2 = IWORK2 + NUC2*NCT134
        ICONT  = ICONT  + NRC2*NUC2
  600  CONTINUE
       DO 700 I = 1, NCONT1
        IADR1 = ((I - 1)*I)/2*NINT34
        IADR2 = (I - 1)*NCT134
        DO 710 J = 1, I*NINT34
         TWOOUT(IADR1 + J) = WORK1(IADR2 + J)
  710   CONTINUE
  700  CONTINUE
      ELSE
       IWORK1 = 1
       IWORK2 = 1
       ICONT  = ISTRCF
       DO 800 I = 1, NSET2
        NUC2 = NUC1X(IOFFVC + I)
        NRC2 = NRC1X(IOFFVC + I)
        CALL XGEMM('N','N',NCT134,NRC2,NUC2,ONE,WORK2(IWORK2),NCT134,
     &             CONTT1(ICONT),NUC2,AZERO,TWOOUT(IWORK1),NCT134)
        IWORK1 = IWORK1 + NRC2*NCT134
        IWORK2 = IWORK2 + NUC2*NCT134
        ICONT  = ICONT  + NRC2*NUC2
  800  CONTINUE
      END IF
      RETURN
      END
