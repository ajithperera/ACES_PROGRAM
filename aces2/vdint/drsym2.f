










      SUBROUTINE DRSYM2(SO,AO,HKABCD,ISYMR,ISYMT,ISYMTS,MULATM,MULE,
     &                  SHABAB,NOABCD,LWRKAO,IPRINT)
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL TRANS
      LOGICAL TWOCEN, THRCEN, FOUCEN, DERONE, DERTWO
      LOGICAL         AEQB,   CEQD,   DIAGAB, DIAGCD, DIACAC,
     &                ONECEN, PQSYM,  DTEST,
     &                TPRIAB, TPRICD, TCONAB, TCONCD
      LOGICAL SHABAB
      INTEGER AND,OR,XOR
      DOUBLE PRECISION
     &        SIGN1X, SIGN1Y, SIGN1Z, SIGN2X, SIGN2Y, SIGN2Z,
     &        SIGN3X, SIGN3Y, SIGN3Z, SIGN4X, SIGN4Y, SIGN4Z
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
C
      DIMENSION AO(LWRKAO), SO(*)
C
      COMMON /DIRPRT/ SIGNDR(4,3), TRANS, NATOMS, IATOMS(4), ISOPDR(4)
      COMMON /EXPCOM/ SIGN1X, SIGN1Y, SIGN1Z, SIGN2X, SIGN2Y, SIGN2Z,
     &                SIGN3X, SIGN3Y, SIGN3Z, SIGN4X, SIGN4Y, SIGN4Z,
     &                NCENT1, NCENT2, NCENT3, NCENT4,
     &                ISO1,   ISO2,   ISO3,   ISO4,
     &                DERONE, DERTWO, TWOCEN, THRCEN, FOUCEN,
     &                NINTYP, NCCINT
      COMMON /INTINF/ THRESH,
     &                NHKTA,  NHKTB,  NHKTC,  NHKTD,
     &                MAXAB,  MAXCD,  JMAX0,
     &                KHKTA,  KHKTB,  KHKTC,  KHKTD,
     &                KHKTAB, KHKTCD, KHABCD,
     &                MHKTA,  MHKTB,  MHKTC,  MHKTD,
     &                MULA,   MULB,   MULC,   MULD,
     &                NORBA,  NORBB,  NORBC,  NORBD, NORBAB, NORBCD,
     &                NUCA,   NUCB,   NUCC,   NUCD,  NUCAB,  NUCCD,
     &                NSETA,  NSETB,  NSETC,  NSETD,
     &                ISTEPA, ISTEPB, ISTEPC, ISTEPD,
     &                NSTRA,  NSTRB,  NSTRC,  NSTRD,
     &                AEQB,   CEQD,
     &                DIAGAB, IAB0X,  IAB0Y,  IAB0Z,
     &                DIAGCD, ICD0X,  ICD0Y,  ICD0Z,
     &                DIACAC, ONECEN, PQSYM,  IPQ0X, IPQ0Y, IPQ0Z,
     &                TPRIAB, TPRICD, TCONAB, TCONCD,
     &                MAXDER, DTEST
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     &                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     &                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     &                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     &                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     &                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     &                IPTXYZ(3,0:7)
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
      IF (IPRINT .GT. 10) THEN
       CALL HEADER('Subroutine DRSYM2',-1)
       WRITE (LUPRI,'(A,3I5)') ' ISYMR, ISYMT, ISYMTS ',
     &                             ISYMR, ISYMT, ISYMTS
      END IF
      IF (NATOMS .GT. 1) THEN
       NWRKAO = 3*(NATOMS + 1)*NCCINT
       IF (NWRKAO .GT. LWRKAO) THEN
        WRITE (LUPRI,'(/,A,2(/,A,I10))')
     &       ' Work space exceeded in DRSYM2.',
     &       ' NWRKAO:', NWRKAO, ' LWRKAO:', LWRKAO
        WRITE (LUPRI,'(A)') ' Increase dimension of WORK1.'
        CALL ERREX 
       END IF
      END IF
      ISOSTR = 1
      DO 100 ICOOR = 1, 3
       ISYTYE = ISYTYP(ICOOR)
       DO 200 IREPE = 0, MAXLOP
        IF (IBTAND(MULE,IBTXOR(IREPE,ISYTYE)) .EQ. 0) THEN
         IF (NATOMS .EQ. 1) THEN
          FACSYM = PT(IBTAND(ISOPDR(1),IREPE))
     &             *SIGNDR(1,ICOOR)*HKABCD
          IAOSTR = (ICOOR - 1)*NCCINT
         ELSE
          FACSYM = HKABCD
          IAOSTR = 3*NATOMS*NCCINT
          DO 300 ICENT = 1, NATOMS
           FACTOR = PT(IBTAND(ISOPDR(ICENT),IREPE))
     &              *SIGNDR(ICENT,ICOOR)
           IOFF   = (3*ICENT + ICOOR - 4)*NCCINT
           IF (ICENT .EQ. 1) THEN
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
           DO 400 INT = 1, NCCINT
            AO(IAOSTR + INT) = FACTOR*AO(IOFF + INT)
  400      CONTINUE
          ELSE
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
          DO 450 INT = 1, NCCINT
           AO(IAOSTR + INT) = AO(IAOSTR + INT)
     &                      + FACTOR*AO(IOFF + INT)
  450     CONTINUE
         END IF
  300   CONTINUE
       END IF
       IF (IPRINT .GT. 20) THEN
        WRITE (LUPRI,'(A,3I5)') ' ICOOR, IREPE, ISYTYE ',
     &                            ICOOR, IREPE, ISYTYE
        WRITE (LUPRI,'(A,F12.6)') ' FACSYM ', FACSYM
       END IF
       CALL SYM2(SO(ISOSTR),AO(IAOSTR+1),ISYMR,ISYMT,ISYMTS,
     &           SHABAB,NOABCD,FACSYM,IREPE,NSOINT,IPRINT)
       ISOSTR = ISOSTR + NSOINT
      END IF
  200 CONTINUE
  100 CONTINUE
      RETURN
      END
