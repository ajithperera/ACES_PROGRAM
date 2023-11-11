










      SUBROUTINE SAOORD(REORD,NBASIS,NIRREP)
      IMPLICIT INTEGER (A-Z)
      INTEGER AND,OR,XOR
      DOUBLE PRECISION PT,FMULT,CENTSH
      LOGICAL BIGVEC,SEGMEN
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


c IMPORTANT: The ECP library contains three routines (doitgr, ecpabl, ecpder)
c that use this parameter but do not include this file. Changing this value
c requires updating those files as well. [Yau: If we linked them together,
c then compiling ecp would require the vdint directory. Alternatively, we
c could put the mxcent.par file into the top-level include directory.]

C MXCENT : Maximum number of atoms currently allowed
C MXCOOR : The number of Cartesian that correspond to MXCENT

      INTEGER MXCENT, MXCOOR
      PARAMETER (MXCENT=200, MXCOOR=3*MXCENT)
C
      DIMENSION REORD(NBASIS,NIRREP)
      DIMENSION ICORB(MXAOVC)
C
      COMMON/BLOCKS/CENTSH(MXSHEL,3),
     &              MAXSHL,BIGVEC,SEGMEN,
     &              NHKTSH(MXSHEL),KHKTSH(MXSHEL),MHKTSH(MXSHEL),
     &              ISTBSH(MXSHEL),NUCOSH(MXSHEL),NORBSH(MXSHEL),
     &              NSTRSH(MXSHEL),NCNTSH(MXSHEL),NSETSH(MXSHEL),
     &              JSTRSH(MXSHEL,MXAOVC),
     &              NPRIMS(MXSHEL,MXAOVC),
     &              NCONTS(MXSHEL,MXAOVC),
     &              IORBSH(MXSHEL,MXAOVC),
     &              IORBSB(MXCORB),NRCSH(MXSHEL)
      COMMON/SYMMET/FMULT(0:7),PT(0:7),
     &              MAXLOP,MAXLOT,MULT(0:7),ISYTYP(3),
     &              ITYPE(8,36),NPARSU(8),NPAR(8),NAOS(8),
     &              NPARNU(8,8),IPTSYM(MXCORB,0:7),
     &              IPTCNT(3*MXCENT,0:7),NCRREP(0:7),
     &              IPTCOR(MXCENT*3),NAXREP(0:7),IPTAX(3),
     &              IPTXYZ(3,0:7)
C
      IBTAND(I,J)=AND(I,J)
      IBTOR(I,J)=OR(I,J)
      IBTXOR(I,J)=XOR(I,J)
C
C  LOOP OVER ALL IRREPS
C
      IND=0
      DO 1000 IRREP=1,NIRREP
C
C  LOOP OVER SHELLS
C
       DO 100 ISHEL=1,MAXSHL 
C
        NORB=NORBSH(ISHEL)
        DO 10 I=1,NORB
         ICORB(I)=IORBSH(ISHEL,I)-1
10      CONTINUE
C
        KHKT=KHKTSH(ISHEL)
        NHKT=NHKTSH(ISHEL)
        MUL=ISTBSH(ISHEL)
C
        DO 200 ICOMP=1,KHKT
C
         ITYN=ITYPE(NHKT,ICOMP)
C
C CHECK WHETHER THERE IS A CONTRIBUTION OF THIS BASIS FUNCTIONS 
C
         IF(IBTAND(MUL,IBTXOR(IRREP-1,ITYN)).EQ.0) THEN 
C
          DO 300 IORB=1,NORB
C
           IND=IND+1
           INDX=ICORB(IORB)+IORB+(ICOMP-1)*NORB
           REORD(INDX,IRREP)=IND
C
300       CONTINUE
C
         ENDIF
C
200     CONTINUE
C
100    CONTINUE
C
1000  CONTINUE
C
      RETURN
      END
