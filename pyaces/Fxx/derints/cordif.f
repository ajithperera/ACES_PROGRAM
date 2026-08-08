










      SUBROUTINE CORDIF(NORBA,NORBB,THRESH,D0X,D0Y,D0Z,
     *                  IPRINT,NUCOA,NUCOB,JSTRA,JSTRB)
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

      LOGICAL DIFX, DIFY, DIFZ, D0X, D0Y, D0Z
      DIMENSION NUCOA(MXAOVC), NUCOB(MXAOVC),
     *          JSTRA(MXAOVC), JSTRB(MXAOVC)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
      COMMON /PRIMIT/ PRIEXP(MXPRIM), PRICCF(MXPRIM,MXCONT),
     *                PRICRX(MXPRIM), PRICRY(MXPRIM), PRICRZ(MXPRIM)
C
      DIFX = .FALSE.
      DIFY = .FALSE.
      DIFZ = .FALSE.
C
C     A - A
C
      DO 100 I = 1, NORBA
      DO 100 J = 1, NUCOA(I)
         IJ  = JSTRA(I) + J
         CRX = PRICRX(IJ)
         CRY = PRICRY(IJ)
         CRZ = PRICRZ(IJ)
         DO 200 K = 1, I
         DO 200 L = 1, NUCOA(K)
            KL   = JSTRA(K) + L
            DIFX = DIFX .OR. ABS(PRICRX(KL)-CRX) .GT. THRESH
            DIFY = DIFY .OR. ABS(PRICRY(KL)-CRY) .GT. THRESH
            DIFZ = DIFZ .OR. ABS(PRICRZ(KL)-CRZ) .GT. THRESH
  200    CONTINUE
  100 CONTINUE
C
C     B - B
C
      IF (.NOT.(DIFX .AND. DIFY .AND. DIFZ)) THEN
         DO 300 I = 1, NORBB
         DO 300 J = 1, NUCOB(I)
            IJ  = JSTRB(I) + J
            CRX = PRICRX(IJ)
            CRY = PRICRY(IJ)
            CRZ = PRICRZ(IJ)
            DO 400 K = 1, I
            DO 400 L = 1, NUCOB(K)
               KL   = JSTRB(K) + L
               DIFX = DIFX .OR. ABS(PRICRX(KL)-CRX) .GT. THRESH
               DIFY = DIFY .OR. ABS(PRICRY(KL)-CRY) .GT. THRESH
               DIFZ = DIFZ .OR. ABS(PRICRZ(KL)-CRZ) .GT. THRESH
  400       CONTINUE
  300    CONTINUE
      END IF
C
C     A - B
C
      IF (.NOT.(DIFX .AND. DIFY .AND. DIFZ)) THEN
         DO 500 I = 1, NORBA
         DO 500 J = 1, NUCOA(I)
            IJ  = JSTRA(I) + J
            CRX = PRICRX(IJ)
            CRY = PRICRY(IJ)
            CRZ = PRICRZ(IJ)
            DO 600 K = 1, NORBB
            DO 600 L = 1, NUCOB(K)
               KL   = JSTRB(K) + L
               DIFX = DIFX .OR. ABS(PRICRX(KL)-CRX) .GT. THRESH
               DIFY = DIFY .OR. ABS(PRICRY(KL)-CRY) .GT. THRESH
               DIFZ = DIFZ .OR. ABS(PRICRZ(KL)-CRZ) .GT. THRESH
  600       CONTINUE
  500    CONTINUE
      END IF
      D0X = .NOT.DIFX
      D0Y = .NOT.DIFY
      D0Z = .NOT.DIFZ
C
      IF (IPRINT .LT. 05) RETURN
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      CALL HEADER('SUBROUTINE CORDIF',-1)
      WRITE (LUPRI, 1010) NORBA, NORBB
      WRITE (LUPRI, 1020) (NUCOA(I), I = 1, NORBA)
      WRITE (LUPRI, 1030) (NUCOB(I), I = 1, NORBB)
      WRITE (LUPRI, 1040) (JSTRA(I), I = 1, NORBA)
      WRITE (LUPRI, 1050) (JSTRB(I), I = 1, NORBB)
      WRITE (LUPRI, 1060) D0X, D0Y, D0Z
 1010 FORMAT(  '  NORB     ',2I7)
 1020 FORMAT(  '  NUCOA:   ',15I7)
 1030 FORMAT(  '  NUCOB:   ',15I7)
 1040 FORMAT(  '  JSTRA:   ',15I7)
 1050 FORMAT(  '  JSTRB:   ',15I7)
 1060 FORMAT(  '  D0X/Y/Z: ',3L5)
      RETURN
      END
