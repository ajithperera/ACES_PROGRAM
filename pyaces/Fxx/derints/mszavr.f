










      SUBROUTINE MSZAVR(KHKTA,KHKTB,IDENA,IDENB,
     &                  ISYMOP,LDIAG,WORK1,MSZ,MSZSCF)
C
C  CALCULATES THE EXPECTATION VALUES FOR THE
C  MAGNETIC SUSCEPTIBILITY
C
CEND
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
      DOUBLE PRECISION MSZ,MSZSCF
C

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
     &           MXCONT*MXCONT)
      PARAMETER (D0 = 0.00 D00)
      LOGICAL LDIAG,SCF,NONHF
      DIMENSION WORK1(1),MSZ(6),MSZSCF(6)
      COMMON/DF/ DSHELL(MXAQNS),DSHELL2(MXAQNS)
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     &                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     &                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     &                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     &                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     &                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     &                IPTXYZ(3,0:7)
      COMMON /ADER/ ADER0 (MXAQNS)
      COMMON /POINTER/IS0000, IS000X, IS000Y, IS000Z,
     &                IS00XX, IS00XY, IS00XZ, IS00YY,
     &                IS00YZ, IS00ZZ, IT0000, IT000X,
     &                IT000Y, IT000Z, IT00XX, IT00XY,
     &                IT00XZ, IT00YY, IT00YZ, IT00ZZ,
     &                IL0000, IL000X, IL000Y, IL000Z,
     &                IL00XX, IL00XY, IL00XZ, IL00YX,
     &                IL00YY, IL00YZ, IL00ZX, IL00ZY,
     &                IL00ZZ,
     &                IA0000, IA0X00, IA0Y00, IA0Z00,
     &                IAXX00, IAXY00, IAXZ00, IAYY00,
     &                IAYZ00, IAZZ00, IA000X, IA000Y,
     &                IA000Z, IA00XX, IA00XY, IA00XZ,
     &                IA00YY, IA00YZ, IA00ZZ, IA0X0X,
     &                IA0X0Y, IA0X0Z, IA0Y0X, IA0Y0Y,
     &                IA0Y0Z, IA0Z0X, IA0Z0Y, IA0Z0Z
      COMMON/POINT2/IQ00XX,IQ00XY,IQ00XZ,IQ00YY,IQ00YZ,IQ00ZZ
      COMMON/DENPOIN/KDEN,KFOC
      COMMON/METHOD/IUHF,SCF,NONHF
      COMMON /GENCON/ NRCA,NRCB,CONTA(MXCONT),CONTB(MXCONT)
C
      DATA AZERO /0.0D0/
C
      IBTAND(I,J) = AND(I,J)
      ITRI(I,J) = MAX(I,J)*(MAX(I,J) - 1)/2 + MIN(I,J)
C
      IOFF=0
      ISKIP=NRCA*NRCB
      IF(LDIAG) ISKIP=NRCA*(NRCA+1)/2
      MAXCMP = ISKIP*KHKTA*KHKTB
C
C EXPECTATION VALUE OF MAGNETIC SUSCEPTIBILITY
C
      XYAVR0=AZERO
      XZAVR0=AZERO
      YZAVR0=AZERO
      XXAVR0=SDOT(MAXCMP,DSHELL,1,WORK1(IL00XX+1),1)
      YYAVR0=SDOT(MAXCMP,DSHELL,1,WORK1(IL00YY+1),1)
      ZZAVR0=SDOT(MAXCMP,DSHELL,1,WORK1(IL00ZZ+1),1)
      IF(ISYTYP(1).EQ.ISYTYP(2))
     &    XYAVR0=SDOT(MAXCMP,DSHELL,1,WORK1(IL00XY+1),1)
      IF(ISYTYP(1).EQ.ISYTYP(3))
     &    XZAVR0=SDOT(MAXCMP,DSHELL,1,WORK1(IL00XZ+1),1)
      IF(ISYTYP(3).EQ.ISYTYP(2))
     &    YZAVR0=SDOT(MAXCMP,DSHELL,1,WORK1(IL00YZ+1),1) 
C
      MSZ(1)=MSZ(1)+XXAVR0
      MSZ(2)=MSZ(2)+XYAVR0
      MSZ(3)=MSZ(3)+XZAVR0
      MSZ(4)=MSZ(4)+YYAVR0
      MSZ(5)=MSZ(5)+YZAVR0
      MSZ(6)=MSZ(6)+ZZAVR0
C
      IF(.NOT.SCF) THEN
C
       IOFF=0
       ISKIP=NRCA*NRCB
       IF(LDIAG) ISKIP=NRCA*(NRCA+1)/2
       MAXCMP = ISKIP*KHKTA*KHKTB
C
C EXPECTATION VALUE OF MAGNETIC SUSCEPTIBILITY
C
       XYAVR0=AZERO
       XZAVR0=AZERO
       YZAVR0=AZERO
       XXAVR0=SDOT(MAXCMP,DSHELL2,1,WORK1(IL00XX+1),1)
       YYAVR0=SDOT(MAXCMP,DSHELL2,1,WORK1(IL00YY+1),1)
       ZZAVR0=SDOT(MAXCMP,DSHELL2,1,WORK1(IL00ZZ+1),1)
       IF(ISYTYP(1).EQ.ISYTYP(2))
     &     XYAVR0=SDOT(MAXCMP,DSHELL2,1,WORK1(IL00XY+1),1)
       IF(ISYTYP(1).EQ.ISYTYP(3))
     &     XZAVR0=SDOT(MAXCMP,DSHELL2,1,WORK1(IL00XZ+1),1)
       IF(ISYTYP(3).EQ.ISYTYP(2))
     &     YZAVR0=SDOT(MAXCMP,DSHELL2,1,WORK1(IL00YZ+1),1) 
C
       MSZSCF(1)=MSZSCF(1)+XXAVR0
       MSZSCF(2)=MSZSCF(2)+XYAVR0
       MSZSCF(3)=MSZSCF(3)+XZAVR0
       MSZSCF(4)=MSZSCF(4)+YYAVR0
       MSZSCF(5)=MSZSCF(5)+YZAVR0
       MSZSCF(6)=MSZSCF(6)+ZZAVR0
C
      ENDIF
C
C ALL DONE, RETURN
      RETURN
      END
