










      SUBROUTINE DASYM3(CAODER,WORK1,KNDAO,NATOMC,
     &                  ONECEN,ICENTA,ICENTB,ICENTC,
     &                  ISYMOP,JSYMOP,MULA,MULB,NHKTA,NHKTB,
     &                  KHKTA,KHKTB,HKAB,LDIAG,LAEQB,
     &                  THRESH,IPRINT)
C
C ARRANGE CALCULATION OF SYMMETRY-ADAPTED INTEGRAL DERIVATIVES
C FROM DISTINCT AO NUCLEAR ATTRACTION INTEGRAL DERIVATIVES
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
      LOGICAL LDIAG, ONECEN, LAEQB
      LOGICAL FULMAT
      LOGICAL DCORD, DCORGD, NOORBT, DOPERT
      CHARACTER NAMEX*6
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

      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN)
C
      DIMENSION CAODER(1),WORK1(1),KNDAO(8)
C
      COMMON /NUCLEIi/ NOORBT(MXCENT),
     &                NUCIND, NUCDEP, NUCPRE(MXCENT), NUCNUM(MXCENT,8),
     &                NUCDEG(MXCENT), ISTBNU(MXCENT), NDCORD,
     &                NDCOOR(MXCOOR), NTRACO, NROTCO, ITRACO(3),
     &                IROTCO(3),
     &                NATOMS, NFLOAT,
     &                IPTGDV(3*MXCENT),
     &                NGDVEC(8), IGDVEC(8)
      COMMON /NUCLEI/ CHARGE(MXCENT), CORD(MXCENT,3),
     &                DCORD(MXCENT,3),DCORGD(MXCENT,3),
     &                DOPERT(0:3*MXCENT)
      COMMON /NUCLEC/ NAMEX(MXCOOR)
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     &                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     &                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     &                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     &                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     &                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     &                IPTXYZ(3,0:7)
      COMMON/GENCON/NRCA,NRCB,CONTA(MXCONT),CONTB(MXCONT)
      COMMON/LSYM/NLENQ(8),NLENT(8)
      COMMON/PERT/NTPERT,NPERT(8),IPERT(8),IXPERT,IYPERT,IZPERT,
     &            IBPERT(3)
C
      DATA ONEM /-1.D0/
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C 
      FULMAT = .TRUE.
      KHKTAB = KHKTA*KHKTB
      NRCAB=NRCA*NRCB
      IF(LDIAG) NRCAB=NRCA*(NRCA+1)/2
      NINTS=KHKTAB*NRCAB*NATOMC
      NMATS = 3*NUCDEP*(MAXLOP+1)
      ICLMX=1
      DO 10 ICL = 1,ICLMX
C
       JCENT = ICENTC
       MULJ = ISTBNU(JCENT)
       IOFF = 1
C
C RUN OVER CARTESIAN DIRECTIONS
C
       DO 20 JCDIR = 1,3
        IDIR = JCDIR - 1
        ISYTYJ=IBPERT(JCDIR)-1
C
C RUN OVER IRREPS OF THE DIFFERENTIATION OPERATOR
C
        DO 30 IREPD = 0,MAXLOP
         IF (IBTAND(MULJ,IBTXOR(IREPD,ISYTYJ)) .EQ. 0) THEN
          FAC=HKAB*PT(IBTAND(ISYTYJ,JSYMOP))
     &        *PT(IBTAND(IREPD,JSYMOP))
          IMAT0 = NMATS
     &          + (IPTCNT(3*(JCENT-1)+IDIR+1,IREPD)-1)*(MAXLOP+1)+1
          IOFFSET=(IPTCNT(3*(JCENT-1)+IDIR+1,IREPD)
     &            -IPERT(IREPD+1)-1)*NLENT(IREPD+1)
          IF (IREPD .EQ. 0) THEN
           CALL SYM1S(CAODER(IOFF),WORK1(KNDAO(IREPD+1)
     &                +IOFFSET),ISYMOP,MULA,MULB,
     &                NHKTA,NHKTB,KHKTA,KHKTB,FAC,LDIAG,LAEQB,
     &                FULMAT,
     &                THRESH,IMAT0,NATOMC,IPRINT,ONEM)
          ELSE
           CALL SYM1N(CAODER(IOFF),WORK1(KNDAO(IREPD+1)
     &                +IOFFSET),IREPD,ISYMOP,MULA,
     &                MULB,NHKTA,NHKTB,KHKTA,KHKTB,-FAC,LDIAG,
     &                LAEQB,
     &                FULMAT,THRESH,IMAT0,NATOMC,IPRINT,ONEM)
          END IF
         END IF
30      CONTINUE
        IOFF = IOFF + NINTS
20     CONTINUE
10    CONTINUE
      RETURN
      END
