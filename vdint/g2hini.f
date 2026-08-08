










      SUBROUTINE G2HINI(WORK1,LWORK1)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
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

      PARAMETER (AZERO = 0.00 D00, ONE = 1.00 D00)
      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN)

c Dimensions of the WORK2 and WORK3 arrays of the VDINT AME. For
c machines that have small stack space, these values must be small.
c Marshall Cory, who has access to CRAY machines with large memory, uses
c 900000 and 500000 is recommended for the others. Small values for
c these arrays will crash the program for basis sets that have large
c primitive spaces. For example the ANO basis sets. A. Perera, 03/2005. 

      PARAMETER (LWORK2 = 300 0000)
      PARAMETER (LWORK3 = 900 0000)
      LOGICAL PATH1, CROSS, DORDER(2), DHOVEC, DHTVEC, DMXVEC,
     *        PQSYM, DTEST, SEGCON, TPRI34, TCON34
      DIMENSION IXDER(6), IYDER(6), IZDER(6), JODD1(3),
     *          ITRI1(0:18), ITRI2(0:18),
     *          JODDIF(27), ITADD(27), IUADD(27), IVADD(27),
     *          FACTOR(27), IOFFHC(27), IOFFCC(27), NRCCPP(27),
     *          IDHO(2,2), IDHT(3,3,2), IDMX(3,3,2,2)
      DIMENSION WORK1(LWORK1)
      DIMENSION IEX(3),IEY(3),IEZ(3),IXADD(3),IYADD(3),IZADD(3)
      COMMON /CWORK2/ WK2LOW, WORK2(LWORK2), WK2HGH
      COMMON /CWORK3/ WK3LOW, WORK3(LWORK3), WK3HGH
      LOGICAL TKTIME
      COMMON /CC2INF/ KHKT12, KHKT34, KH1234,
     *                NORB3, NORB4, NORB12, NORB34, NO1234,
     *                TPRI34, TCON34,
     *                IPQ0X, IPQ0Y, IPQ0Z, INCRMT, INCRMU, INCRMV,
     *                ISTRET, ISTREU, ISTREV, LVAL34, MVAL34, NVAL34,
     *                NCCPP, NUC3, NUC4, NUC34, NSET3, NSET4,
     *                IADDPV, ISTRCF, IOFFVC, IODALL(MXAQN**2), PQSYM,
     *                ICMP34, IPATH, PATH1, MAXDER, CROSS, DTEST,
     *                SEGCON, IPRINT
      COMMON /PRIVEC/ NPRM12(MXAOSQ), NPRM34(MXAOSQ),
     *                LSTP12(MXAOSQ), LSTP34(MXAOSQ),
     *                LSTO12(MXAOSQ), LSTO34(MXAOSQ),
     *                LASTP(MXAOSQ),  LASTO(MXAOSQ)
      COMMON /INTADR/ IWKAO, IWKSO, IWKHHS, IWK1HH, IWK1HC, IWKLST
      COMMON /CCFCOM/ CONT1 (MXCONT*MXAOVC), CONT2 (MXCONT*MXAOVC),
     *                CONT3 (MXCONT*MXAOVC), CONT4 (MXCONT*MXAOVC),
     *                CONTT1(MXCONT*MXAOVC), CONTT2(MXCONT*MXAOVC),
     *                CONTT3(MXCONT*MXAOVC), CONTT4(MXCONT*MXAOVC),
     *                NUC1X (MXAOVC),        NUC2X (MXAOVC),
     *                NUC3X (MXAOVC),        NUC4X (MXAOVC),
     *                NRC1X (MXAOVC),        NRC2X (MXAOVC),
     *                NRC3X (MXAOVC),        NRC4X (MXAOVC)
      COMMON /ODCADR/ IE1X00, IE1X10, IE1X01, IE1X20, IE1X11, IE1X02,
     *                IE1Y00, IE1Y10, IE1Y01, IE1Y20, IE1Y11, IE1Y02,
     *                IE1Z00, IE1Z10, IE1Z01, IE1Z20, IE1Z11, IE1Z02,
     *                IE2X00, IE2X10, IE2X01, IE2X20, IE2X11, IE2X02,
     *                IE2Y00, IE2Y10, IE2Y01, IE2Y20, IE2Y11, IE2Y02,
     *                IE2Z00, IE2Z10, IE2Z01, IE2Z20, IE2Z11, IE2Z02,
     *                IELAST, ILST12
      LOGICAL         DZER
      COMMON /DERZER/ FZERO, DZER, IZERO
      COMMON /DHCINF/ IDHC(10), NDHC(10)
      COMMON /DHODIR/ DHOVEC(18)
      COMMON /DHOADR/ IHOVEC(18)
      COMMON /DHOFAC/ FHOVEC(18)
      COMMON /DHTDIR/ DHTVEC(9)
      COMMON /DHTADR/ IHTVEC(9)
      COMMON /DHTFAC/ FHTVEC(9)
      COMMON /DMXDIR/ DMXVEC(36)
      COMMON /DMXADR/ IMXVEC(36)
      COMMON /DMXFAC/ FMXVEC(36)
      COMMON/INT2H/NCOOR,JODDIF,ITADD,IUADD,IVADD,IOFFCC,FACTOR,
     *     IEXADR, IEYADR, IEZADR, NWORK, IADDCC, NRCCPP,
     *     IOFFHC
      COMMON/INT2H1/IEX,IEY,IEZ
      DATA ITRI1 /0,1,3,6,10,15,21,28,36,45,55,66,78,91,105,
     *            120,136,153,171/
     *     ITRI2 /0,1,4,10,20,35,56,84,120,165,220,286,364,455,
     *            560,680,816,969,1140/
      DATA IXDER /0,1,1,2,2,2/
     *     IYDER /0,1,0,2,1,0/
     *     IZDER /0,0,1,0,1,2/
      IHCADR(I,J,K) = NCCPP*(ITRI2(I+J+K) + ITRI1(J+K) + K)
C
C   HOVEC(18) 
C
      DATA IDHO /3,0,  12,6/
C
C     ARRANGEMENT OF VECTORS HOVEC(18)
C
C     1   XP00 YP00 ZP00
C     4   XQ00 YQ00 ZQ00
C     7   XXPP XYPP XZPP YYPP YZPP ZZPP
C     13  XXQQ XYQQ XZQQ YYQQ YZQQ ZZQQ
C
C                ***** HTVEC(9) *****
C
      DATA IDHT / 1, 2, 3, 4, 5, 6, 7, 8, 9,
     *            1, 4, 7, 2, 5, 8, 3, 6, 9/
C
C     ARRANGEMENT OF VECTORS HTVEC(9)
C
C     1   XXPQ XYPQ XZPQ
C     4   XYQP YYPQ YZPQ
C     7   XZQP YZQP ZZPQ
C
C                ***** MXVEC(36) *****
C
      DATA IDMX /19, 22, 25, 20, 23, 26, 21, 24, 27,
     *           28, 31, 34, 29, 32, 35, 30, 33, 36,
     *            1,  2,  3,  4,  5,  6,  7,  8,  9,
     *           10, 11, 12, 13, 14, 15, 16, 17, 18/
C
C     ARRANGEMENT OF VECTORS MXVEC(36)
C
C     1   XXPC XYPC XZPC
C     4   XYCP YYPC YZPC
C     7   XZCP YZCP ZZPC
C
C     10  XXPD XYPD XZPD
C     13  XYDP YYPD YZPD
C     16  XZDP YZDP ZZPD
C
C     19  XXAQ XYAQ XZAQ
C     22  XYQA YYAQ YZAQ
C     25  XZQA YZQA ZZAQ
C
C     28  XXBQ XYBQ XZBQ
C     31  XYQB YYBQ YZBQ
C     34  XZQB YZQB ZZBQ
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
C  INITIALIZATION - ENTRY G2HINI 
C
C     IEX, IEY, AND IEZ
C
      IF (PATH1) THEN
       IEXADR=IE2X00
       IEYADR=IE2Y00
       IEZADR=IE2Z00
       IEX(1)=IE2X10
       IEX(2)=IEXADR
       IEX(3)=IEXADR
       IEY(1)=IEYADR
       IEY(2)=IE2Y10
       IEY(3)=IEYADR
       IEZ(1)=IEZADR
       IEZ(2)=IEZADR
       IEZ(3)=IE2Z10
      ELSE
       IEXADR=IE1X00
       IEYADR=IE1Y00
       IEZADR=IE1Z00
       IEX(1)=IE1X10
       IEX(2)=IEXADR
       IEX(3)=IEXADR
       IEY(1)=IEYADR
       IEY(2)=IE1Y10
       IEY(3)=IEYADR
       IEZ(1)=IEZADR
       IEZ(2)=IEZADR
       IEZ(3)=IE1Z10
      END IF
C
C  ITADD, IUADD, AND IVADD

      ITADD(1)=1
      ITADD(2)=0
      ITADD(3)=0
      IUADD(1)=0
      IUADD(2)=1
      IUADD(3)=0
      IVADD(1)=0
      IVADD(2)=0
      IVADD(3)=1
C
C     NCOOR, ITADD, IUADD, IVADD, JODDIF, IOFFCC & FACTOR
C
      ICOOR = 0
      JDHO=IDHO(2,1)
      DO 100 IORDER=1,6
       IADRHO=JDHO+IORDER
       IADRHC=1
       ICOOR = ICOOR + 1
       NRCCPP(ICOOR) = NDHC(IADRHC)
       IOFFHC(ICOOR) = IDHC(IADRHC)
       JODDIF(ICOOR) = 0
       IOFFCC(ICOOR) = IWKAO + IHOVEC(IADRHO)
100   CONTINUE
      NCOOR = 6 
      NWORK = NCOOR*KHKT12*NCCPP
      IWKMAX = IWKLST + NWORK
      IF (IWKMAX .GT. LWORK1) THEN
       WRITE (LUPRI, 1000) IWKMAX, IWKLST
       CALL ERREX
      END IF
      IADDCC = NORB12*KH1234
C
      MCOOR=6
      RETURN
 1000 FORMAT (//,1X,' Work space requirement ',I8,' exceeds ',
     *        ' current limit ',I8,' of WORK1.')
      END
