










      SUBROUTINE FSYM3B(IRREPX,FAOA,FAOB,NBAST,IPRINT,FSOA,
     &                  FSOB,NIR,ANTI)
C
C TRANSFORMS THE NON-SYMMETRIC FOCK MATRIX DERIVATIVES
C FROM THE PRIMITVE AO TO THE SYMMETRY-ADAPTED AO BASIS
C
CEND
C
C  MAY/91 JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SHARE
      LOGICAL SCF,NONHF
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


c IMPORTANT: The ECP library contains three routines (doitgr, ecpabl, ecpder)
c that use this parameter but do not include this file. Changing this value
c requires updating those files as well. [Yau: If we linked them together,
c then compiling ecp would require the vdint directory. Alternatively, we
c could put the mxcent.par file into the top-level include directory.]

C MXCENT : Maximum number of atoms currently allowed
C MXCOOR : The number of Cartesian that correspond to MXCENT

      INTEGER MXCENT, MXCOOR
      PARAMETER (MXCENT=200, MXCOOR=3*MXCENT)
      DIMENSION FSOA(1),FSOB(1)
      COMMON /SHELLSi/ KMAX,
     &                NHKT(MXSHEL),   KHKT(MXSHEL), MHKT(MXSHEL),
     &                ISTBAO(MXSHEL), NUCO(MXSHEL), JSTRT(MXSHEL),
     &                NSTRT(MXSHEL),  MST(MXSHEL),  NCENT(MXSHEL),
     &                NRCO(MXSHEL), NUMCF(MXSHEL),
     &                NBCH(MXSHEL),   KSTRT(MXSHEL)
      COMMON /SHELLS/ CENT(MXSHEL,3), SHARE(MXSHEL)
      COMMON /PINCOM/ IPIND(MXCORB), IBLOCK(MXCORB), INDGEN(MXCORB)
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     &                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     &                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     &                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     &                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     &                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     &                IPTXYZ(3,0:7)
      DIMENSION FAOA(NBAST*(NBAST+1)/2,0:NIR),
     &          FAOB(NBAST*(NBAST+1)/2,0:NIR)
      COMMON/METHOD/IUHF,SCF,NONHF
      COMMON/BASSYM/NBAS(8),NBASIS,NBASSQ,NBASTT
      COMMON/LSYM/NLENQ(8),NLENT(8)
C
      DATA AZERO,ONE,TWO /0.D0,1.0D0,2.D0/
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
      CALL ZERO(FSOA,NLENT(IRREPX))
      IF(IUHF.NE.0) CALL ZERO(FSOB,NLENT(IRREPX))
C
C  LOOP OVER ALL IRREPS 
C
      ISOFF=0
      ISTR=1
      DO 100 IREPA = 0, MAXLOP
       IREPB=IBTXOR(IRREPX-1,IREPA)
       NORBI=NAOS(IREPA+1)
       NORBJ=NAOS(IREPB+1)
       IF(IREPA.LE.IREPB) GO TO 110
       IF(MIN(NORBI,NORBJ).EQ.0) GOTO 110
       JSTR=1
       DO 99 IRREP=1,IREPB
        JSTR=JSTR+NAOS(IRREP)
99     CONTINUE
       DO 200 I=ISTR,ISTR+NORBI-1
        IA=IBTSHR(IPIND(I),16)
        NA=IBTSHR(IPIND(I),8) - IA*256
        IOFF=KSTRT(IA)
        MULA=ISTBAO(IA)
        INDA=IOFF + NA
        DO 300 J=JSTR,JSTR+NORBJ-1
         IB=IBTSHR(IPIND(J),16)
         NB=IBTSHR(IPIND(J),8) - IB*256
         JOFF=KSTRT(IB)
         NHKTB=NHKT(IB)
         KHKTB=KHKT(IB)
         MULB=ISTBAO(IB)
         MAB=IBTOR(MULA,MULB)
         ISOFF=ISOFF + 1
         FAIJ=AZERO
         FBIJ=AZERO
         FACT=ONE
         DO 400 ISYMOP = 0, MAXLOT
          FAC=PT(IBTAND(ISYMOP,IREPB))*ANTI
          INDB=JOFF + NB
          IF(INDA.GE.INDB) THEN
           FAC=FAC*PT(IBTAND(IRREPX-1,ISYMOP))*ANTI
          ENDIF
          INDM=MAX(INDA,INDB)
          IND=(INDM*(INDM-3))/2+INDA+INDB
          FAIJ=FAIJ+FAC*FAOA(IND,ISYMOP)
          IF(IUHF.NE.0) THEN
           FBIJ=FBIJ+FAC*FAOB(IND,ISYMOP)
          ENDIF
400      CONTINUE
         FSOA(ISOFF)=FAIJ*FACT
         IF(IUHF.NE.0) FSOB(ISOFF)=FBIJ*FACT
300     CONTINUE
200    CONTINUE
110   CONTINUE
      ISTR=ISTR+NORBI
100   CONTINUE
C
      RETURN
      END
