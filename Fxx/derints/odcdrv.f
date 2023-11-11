










      SUBROUTINE ODCDRV(NHKT1,NHKT2,NSET1,NSET2,ISTEP1,ISTEP2,
     *                  NUC1,NUC2,NUC12,NORB1,NORB2,
     *                  NUCO1,NUCO2,NRCO1,NRCO2,JSTR1,JSTR2,
     *                  SIGN1X,SIGN1Y,SIGN1Z,SIGN2X,SIGN2Y,SIGN2Z,
     *                  COR1X,COR1Y,COR1Z,COR2X,COR2Y,COR2Z,
     *                  SAMEX,SAMEY,SAMEZ,I0X,I0Y,I0Z,TPRIAB,
     *                  ONECEN,DO1,DO2,BIGVEC,SEGMEN,DTEST,
     *                  ITYPE,THRESH,MAXDER,IPRINT)
C
C     TUH Apr 11 1988
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)

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

      PARAMETER (ONE = 1.00 D00)
      LOGICAL SAMEX, SAMEY, SAMEZ, ONECEN, DO1, DO2, TPRIAB,
     *        BIGVEC, SEGMEN, DTEST
      DIMENSION NUCO1(MXAOVC), NUCO2(MXAOVC),
     *          NRCO1(MXCONT), NRCO2(MXCONT),
     *          JSTR1(MXAOVC), JSTR2(MXAOVC)
      LOGICAL TKTIME
      IF (BIGVEC) THEN
         CALL CORDIF(NSET1,NSET2,THRESH,SAMEX,SAMEY,SAMEZ,IPRINT,
     *               NUCO1,NUCO2,JSTR1,JSTR2)
      ELSE
         DIFX  = SIGN1X*COR1X - SIGN2X*COR2X
         DIFY  = SIGN1Y*COR1Y - SIGN2Y*COR2Y
         DIFZ  = SIGN1Z*COR1Z - SIGN2Z*COR2Z
         SAMEX = ABS(DIFX) .LT. THRESH
         SAMEY = ABS(DIFY) .LT. THRESH
         SAMEZ = ABS(DIFZ) .LT. THRESH
      END IF
      I0X = 0
      I0Y = 0
      I0Z = 0
      IF (SAMEX) I0X = 1
      IF (SAMEY) I0Y = 1
      IF (SAMEZ) I0Z = 1
C
C     ****************************************
C     ***** Overlap Distribution Vectors *****
C     ****************************************
C
      CALL ODCVEC(NUC1,NUC2,NORB1,NORB2,NSET1,NSET2,NUCO1,NUCO2,
     *            NRCO1,NRCO2,JSTR1,JSTR2,NUC12,TPRIAB,
     *            SIGN1X,SIGN1Y,SIGN1Z,
     *            SIGN2X,SIGN2Y,SIGN2Z,
     *            THRESH,ITYPE,IPRINT)

      IF (NUC12 .EQ. 0) RETURN
C
C    *********************************************
C    ***** Overlap Distribution Coefficients *****
C    *********************************************
C
      CALL EXCOEF_DERINTS
     *            (BIGVEC,ONECEN,SEGMEN,NHKT1,NHKT2,ISTEP1,ISTEP2,NUC12,
     *            I0X,I0Y,I0Z,DIFX,DIFY,DIFZ,SIGN1X,SIGN1Y,SIGN1Z,
     *            SIGN2X,SIGN2Y,SIGN2Z,THRESH,ITYPE,MAXDER,
     *            DO1,DO2,DTEST,IPRINT)
      RETURN
      END
