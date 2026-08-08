










      SUBROUTINE SET1(NRSS,I2,A2,I1,A1)
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

      parameter (nh4=4*nht-3)
      DIMENSION I2(*),A2(*),I1(50),A1(*)
      COMMON /TAB/ NNNGAM
      COMMON/NULL_COM/FACT(nh4),RFACT(nh4),FACTM(nh4),RFACTM(nh4)
     & ,MAA,IFD1,IFD2,KCD,KBCD,NHCD,NHBCD,NNC,IFPL(3),
     & NN1,NN2,NHKTA,NHKTB,NHKTC,NHKTD,KHKTA,KHKTB,KHKTC,KHKTD,NNB
      COMMON /REP/ NEWIND(MXCBF) , MSTOLD(8)
      COMMON /CONI/ NUCR,NUCS,NUCA,NUCQ,
     & NRCR,NRCS,NRCA,NRCQ,JSTR,JSTS,JSTA,JSTQ,
     & IFD1X,IFD2X,NUCAQ,NUCRS,NAQRS,NCQRS,KQ,KR,KS,IVA
      NSB=NNNGAM + 3*(NNB+NNNGAM-1)
      NSB=MAX0(NSB,8)
      NRSR = (I1(50)-I1(7))/(NSB*NUCAQ)
      IF(NRSR .LE. 0) GO TO 666
      NPASS = 1 + NUCRS/(NRSR+1)
C      FOR TESTING PURPOSES
C      NPASS=MIN0(NUCRS,NUCRS)
C
  667 NRSS=(NUCRS-1)/NPASS + 1
      KQ = 1
      KR = 1
      KS = 1
      IVA= 1
      RETURN
      ENTRY SET2(NRSS,ITAG,I2,A2,I1,A1)
      NSB=NNNGAM + 3*(NNB+NNNGAM-1)
      NSB=MAX0(NSB,8)
      NRSR = (I1(50)-I1(7))/(NSB*NUCAQ)
      IF(NRSR .LE. 0) GO TO 666
      NPASS = 1 + NUCRS/(NRSR+1)
C
      ITAG=1
      IF(NPASS .GT.NUCR) THEN
      ITAG=3
      NPASS = NUCRS
      ELSE IF(NPASS .GT. 1) THEN
      ITAG=2
      NPASS = NUCR
      ENDIF
      GO TO 667
  666 WRITE (6,*) ' COMMON/W1/ IS TOO SMALL. INCREASE, AND TRY AGAIN'
      WRITE (6,*) ' NRSR,I1(7),I1(50), NSB,NUCAQ'
      WRITE (6,6666) NRSR,I1(7),I1(50),NSB,NUCAQ
 6666 FORMAT (5I10)
      CALL ERREX
      END
