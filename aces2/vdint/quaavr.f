










      SUBROUTINE QUAAVR(KHKTA,KHKTB,IDENA,IDENB,NCENTA,NCENTB,
     *                  ISYMOP,SBX,SBY,SBZ,LDIAG,ONECEN,DIFDIP,
     *                  WORK1,QUADME,DQUADE,NCOORD)
C
C     CALCULATES THE VARIOUS EXPECTATION VALUES FOR QUADRUPOLE MOMENT
C     AND QUADRUPOLE MOMENT DERIVATIVES
C
CEND
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
      PARAMETER (LUCMD = 5, LUPRI = 6)

c IMPORTANT: The ECP library contains three routines (doitgr, ecpabl, ecpder)
c that use this parameter but do not include this file. Changing this value
c requires updating those files as well. [Yau: If we linked them together,
c then compiling ecp would require the vdint directory. Alternatively, we
c could put the mxcent.par file into the top-level include directory.]

C MXCENT : Maximum number of atoms currently allowed
C MXCOOR : The number of Cartesian that correspond to MXCENT

      INTEGER MXCENT, MXCOOR
      PARAMETER (MXCENT=200, MXCOOR=3*MXCENT)

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

      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN*
     *           MXCONT*MXCONT)
      PARAMETER (D0 = 0.00 D00)
      LOGICAL DIFDIP, ONECEN, LDIAG
      DIMENSION WORK1(1),QUADME(6),DQUADE(6,NCOORD),DDIPE(1,1)
      COMMON/DF/ DSHELL(MXAQNS),FSHELL(MXAQNS)
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     *                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     *                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     *                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     *                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     *                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     *                IPTXYZ(3,0:7)
      COMMON /ADER/ ADER0 (MXAQNS)
      COMMON /POINTER/IS0000, IS000X, IS000Y, IS000Z,
     *                IS00XX, IS00XY, IS00XZ, IS00YY,
     *                IS00YZ, IS00ZZ, IT0000, IT000X,
     *                IT000Y, IT000Z, IT00XX, IT00XY,
     *                IT00XZ, IT00YY, IT00YZ, IT00ZZ,
     *                ID0000, ID000X, ID000Y, ID000Z,
     *                ID00XX, ID00XY, ID00XZ, ID00YX,
     *                ID00YY, ID00YZ, ID00ZX, ID00ZY,
     *                ID00ZZ,
     *                IA0000, IA0X00, IA0Y00, IA0Z00,
     *                IAXX00, IAXY00, IAXZ00, IAYY00,
     *                IAYZ00, IAZZ00, IA000X, IA000Y,
     *                IA000Z, IA00XX, IA00XY, IA00XZ,
     *                IA00YY, IA00YZ, IA00ZZ, IA0X0X,
     *                IA0X0Y, IA0X0Z, IA0Y0X, IA0Y0Y,
     *                IA0Y0Z, IA0Z0X, IA0Z0Y, IA0Z0Z
      COMMON/POINT2/IQ00XX,IQ00XY,IQ00XZ,IQ00YY,IQ00YZ,IQ00ZZ
      COMMON/DENPOIN/KDEN,KFOC
      COMMON /GENCON/ NRCA,NRCB,CONTA(MXCONT),CONTB(MXCONT)
C
      DATA AZERO /0.0D0/
C
      IBTAND(I,J) = AND(I,J)
      ITRI(I,J) = MAX(I,J)*(MAX(I,J) - 1)/2 + MIN(I,J)
C
      NAX = 3*NCENTA - 2
      NAY = 3*NCENTA - 1
      NAZ = 3*NCENTA
      NBX = 3*NCENTB - 2
      NBY = 3*NCENTB - 1
      NBZ = 3*NCENTB
      IX  = IPTAX(1)
      IY  = IPTAX(2)
      IZ  = IPTAX(3)
C
      IOFF=0
      ISKIP=NRCA*NRCB
      IF(LDIAG) ISKIP=NRCA*(NRCA+1)/2
      MAXCMP = ISKIP*KHKTA*KHKTB
C
C     Expectation value of undifferentiated elements
C
      XYAVR0=AZERO
      XZAVR0=AZERO
      YZAVR0=AZERO
      XXAVR0=SDOT(MAXCMP,DSHELL,1,WORK1(IQ00XX+1),1)
      YYAVR0=SDOT(MAXCMP,DSHELL,1,WORK1(IQ00YY+1),1)
      ZZAVR0=SDOT(MAXCMP,DSHELL,1,WORK1(IQ00ZZ+1),1)
      IF(ISYTYP(1).EQ.ISYTYP(2))
     *    XYAVR0=SDOT(MAXCMP,DSHELL,1,WORK1(IQ00XY+1),1)
      IF(ISYTYP(1).EQ.ISYTYP(3))
     *    XZAVR0=SDOT(MAXCMP,DSHELL,1,WORK1(IQ00XZ+1),1)
      IF(ISYTYP(3).EQ.ISYTYP(2))
     *    YZAVR0=SDOT(MAXCMP,DSHELL,1,WORK1(IQ00YZ+1),1) 
C
      QUADME(1)=QUADME(1)+XXAVR0
      QUADME(2)=QUADME(2)+XYAVR0
      QUADME(3)=QUADME(3)+XZAVR0
      QUADME(4)=QUADME(4)+YYAVR0
      QUADME(5)=QUADME(5)+YZAVR0
      QUADME(6)=QUADME(6)+ZZAVR0
C
c      IF (DIFDIP) THEN
c       IF (ONECEN) THEN
c        DO 300 IREP = 0, MAXLOP
c         IAX = IPTCNT(NAX,IREP)
c         IAY = IPTCNT(NAY,IREP)
c         IAZ = IPTCNT(NAZ,IREP)
c         IF (ISYTYP(1) .EQ. IREP .AND. IAX .GT. 0) THEN
c          DDIPE(IX,IAX) = DDIPE(IX,IAX) - SAVR0
c         END IF
c         IF (ISYTYP(2) .EQ. IREP .AND. IAY .GT. 0) THEN
c          DDIPE(IY,IAY) = DDIPE(IY,IAY) - SAVR0
c         END IF
c         IF (ISYTYP(3) .EQ. IREP .AND. IAZ .GT. 0) THEN
c          DDIPE(IZ,IAZ) = DDIPE(IZ,IAZ) - SAVR0
c         END IF
c  300   CONTINUE
c       ELSE
c        XAVRX=-SDOT(MAXCMP,DSHELL,1,WORK1(ID00XX+1),1)
c        XAVRY=-SDOT(MAXCMP,DSHELL,1,WORK1(ID00YX+1),1)
c        XAVRZ=-SDOT(MAXCMP,DSHELL,1,WORK1(ID00ZX+1),1)
c        YAVRX=-SDOT(MAXCMP,DSHELL,1,WORK1(ID00XY+1),1)
c        YAVRY=-SDOT(MAXCMP,DSHELL,1,WORK1(ID00YY+1),1)
c        YAVRZ=-SDOT(MAXCMP,DSHELL,1,WORK1(ID00ZY+1),1)
c        ZAVRX=-SDOT(MAXCMP,DSHELL,1,WORK1(ID00XZ+1),1)
c        ZAVRY=-SDOT(MAXCMP,DSHELL,1,WORK1(ID00YZ+1),1)
c        ZAVRZ=-SDOT(MAXCMP,DSHELL,1,WORK1(ID00ZZ+1),1)
c        DO 500 IREP = 0, MAXLOP
c         CHI  = PT(IBTAND(ISYMOP,IREP))
c         CSBX = CHI*SBX
c         CSBY = CHI*SBY
c         CSBZ = CHI*SBZ
c         IAX = IPTCNT(NAX,IREP)
c         IAY = IPTCNT(NAY,IREP)
c         IAZ = IPTCNT(NAZ,IREP)
c         IBX = IPTCNT(NBX,IREP)
c         IBY = IPTCNT(NBY,IREP)
c         IBZ = IPTCNT(NBZ,IREP)
c         IF (ISYTYP(1) .EQ. IREP) THEN
c          IF (IAX.GT.0)
c     *    DDIPE(IX,IAX) = DDIPE(IX,IAX) + XAVRX
c          IF (IAY.GT.0)
c     *    DDIPE(IX,IAY) = DDIPE(IX,IAY) + XAVRY
c          IF (IAZ.GT.0)
c     *    DDIPE(IX,IAZ) = DDIPE(IX,IAZ) + XAVRZ
c          IF (IBX.GT.0)
c     *    DDIPE(IX,IBX) = DDIPE(IX,IBX) - CSBX*XAVRX - SAVR0
c          IF (IBY.GT.0)
c     *    DDIPE(IX,IBY) = DDIPE(IX,IBY) - CSBY*XAVRY
c          IF (IBZ.GT.0)
c     *    DDIPE(IX,IBZ) = DDIPE(IX,IBZ) - CSBZ*XAVRZ
c         END IF
c         IF (ISYTYP(2) .EQ. IREP) THEN
c          IF (IAX.GT.0)
c     *    DDIPE(IY,IAX) = DDIPE(IY,IAX) + YAVRX
c          IF (IAY.GT.0)
c     *    DDIPE(IY,IAY) = DDIPE(IY,IAY) + YAVRY
c          IF (IAZ.GT.0)
c     *    DDIPE(IY,IAZ) = DDIPE(IY,IAZ) + YAVRZ
c          IF (IBX.GT.0)
c     *    DDIPE(IY,IBX) = DDIPE(IY,IBX) - CSBX*YAVRX
c          IF (IBY.GT.0)
c     *    DDIPE(IY,IBY) = DDIPE(IY,IBY) - CSBY*YAVRY - SAVR0
c          IF (IBZ.GT.0)
c     *    DDIPE(IY,IBZ) = DDIPE(IY,IBZ) - CSBZ*YAVRZ
c         END IF
c         IF (ISYTYP(3) .EQ. IREP) THEN
c          IF (IAX.GT.0)
c     *    DDIPE(IZ,IAX) = DDIPE(IZ,IAX) + ZAVRX
c          IF (IAY.GT.0)
c     *    DDIPE(IZ,IAY) = DDIPE(IZ,IAY) + ZAVRY
c          IF (IAZ.GT.0)
c     *    DDIPE(IZ,IAZ) = DDIPE(IZ,IAZ) + ZAVRZ
c          IF (IBX.GT.0)
c     *    DDIPE(IZ,IBX) = DDIPE(IZ,IBX) - CSBX*ZAVRX
c          IF (IBY.GT.0)
c     *    DDIPE(IZ,IBY) = DDIPE(IZ,IBY) - CSBY*ZAVRY
c          IF (IBZ.GT.0)
c     *    DDIPE(IZ,IBZ) = DDIPE(IZ,IBZ) - CSBZ*ZAVRZ - SAVR0
c         END IF
c  500   CONTINUE
c       END IF
c      END IF
      RETURN
      END
