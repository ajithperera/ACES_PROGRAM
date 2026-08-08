











      SUBROUTINE GSYM3(SO,AO,ISYMR,ISYMT,ISYMTS,SHABAB,NOABCD,
     &                 HKABCD,IPRINT)
C
C     TAKE BLOCK OF DISTINT AO TWO-ELECTRON GIAO INTEGRALS AND 
C     GENERATE SYMMETRIZED CONTRIBUTIONS TO SO GIAO INTEGRALS. 
C
CEND
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


c IMPORTANT: The ECP library contains three routines (doitgr, ecpabl, ecpder)
c that use this parameter but do not include this file. Changing this value
c requires updating those files as well. [Yau: If we linked them together,
c then compiling ecp would require the vdint directory. Alternatively, we
c could put the mxcent.par file into the top-level include directory.]

C MXCENT : Maximum number of atoms currently allowed
C MXCOOR : The number of Cartesian that correspond to MXCENT

      INTEGER MXCENT, MXCOOR
      PARAMETER (MXCENT=200, MXCOOR=3*MXCENT)
      LOGICAL SHABAB
      LOGICAL         AEQB,   CEQD,   DIAGAB, DIAGCD, DIACAC,
     &                ONECEN, PQSYM,  DTEST,
     &                TPRIAB, TPRICD, TCONAB, TCONCD
C
      LOGICAL TWOCEN,THRCEN,FOUCEN,DERONE,DERTWO
      DIMENSION AO(1), SO(1)
C
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
      COMMON/PERT/NTPERT,NPERT(8),IPERT(8),IXPERT,IYPERT,IZPERT,
     &            IXYPERT(3)
      COMMON/EXPCOM/SIGNXYZ(12),NCENT1,NCENT2,NCENT3,NCENT4,
     &              ISON(4),DERONE,DERTWO,TWOCEN,THRCEN,FOUCEN,
     &              NINTYP,NCCINT
      COMMON/GOFF/IGOFF
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
C     LOOP OVER COMPONENTS
C
      IAOFF = 1
      ISOFF = 1 
      NCCINT2=NCCINT*IGOFF  
      DO 100 ICOMPA = 1,KHKTA
       ITYNA = ITYPE(NHKTA,ICOMPA)
       KHKTBB = KHKTB
       IF (DIAGAB) KHKTBB = ICOMPA
       DO 200 ICOMPB = 1,KHKTBB
        ITYNB = ITYPE(NHKTB,ICOMPB)
        DO 300 ICOMPC = 1,KHKTC
         ITYNC = ITYPE(NHKTC,ICOMPC)
         KHKTDD = KHKTD
         IF (DIAGCD) KHKTDD = ICOMPC
         DO 400 ICOMPD = 1,KHKTDD
          ITYND = ITYPE(NHKTD,ICOMPD)
C
C WHAT IS THE MEANING OF THIS STATEMENT ?
C
c          IF (SHABAB .AND. (ICOMPA .LT. ICOMPC .OR.
c     &       (ICOMPA .EQ. ICOMPC .AND. ICOMPB .LT. ICOMPD))) GO TO 240
C
C LOOP OVER IRREPS
C
          DO 500 IREPA = 0, MAXLOP
           IF (IBTAND(MULA,IBTXOR(IREPA,ITYNA)) .NE. 0) GOTO 500
           DO 600 IREPB = 0, MAXLOP
            IF (IBTAND(MULB,IBTXOR(IREPB,ITYNB)) .NE. 0) GOTO 600
            IRPAB = IBTXOR(IREPB,IREPA)
            FACB   = PT(IBTAND(ISYMR,IBTXOR(IREPB,ITYNB)))*HKABCD
            DO 700 IREPC = 0, MAXLOP
             IF (IBTAND(MULC,IBTXOR(IREPC,ITYNC)) .NE. 0) GOTO 700
             IRPABC = IBTXOR(IREPC,IRPAB)
             FACBC=PT(IBTAND(ISYMT,IBTXOR(IREPC,ITYNC)))*FACB
             DO 800 IREPD = 0,MAXLOP
              IF (IBTAND(MULD,IBTXOR(IREPD,ITYND)) .NE. 0) GOTO 800
              IRPABCD = IBTXOR(IREPC,IRPABC)
              FAC=PT(IBTAND(ISYMTS,IBTXOR(IREPD,ITYND)))*FACBC
C             IF (IPRINT .GT. 25) THEN
C               WRITE (LUPRI,'(A,4I5,5x,4I5)') ' comps, IREPE/A/B/C ',
C    *                  ICOMPA, ICOMPB, ICOMPC, ICOMPD,
C    *                  IREPE, IREPA, IREPB, IREPC
C               WRITE(LUPRI,'(A,6F8.4)') 'FACTOR, FAC', FACTOR, FAC
C             END IF
C
C  Compute SO integral contributions from this block
C  of contracted AO integrals.
C
              ISOADR=ISOFF
              IAOADR=IAOFF
C
CDIR$ IVDEP
              DO 900 INT = 1,NOABCD
               INT1=ISOADR+INT-1
               INT11=IAOADR+INT-1
               SO(INT1)=FAC*AO(INT11)
               SO(INT1+NCCINT2)=FAC*AO(INT11+NCCINT)
               SO(INT1+2*NCCINT2)=FAC*AO(INT11+2*NCCINT)
               SO(INT1+3*NCCINT2)=FAC*AO(INT11+3*NCCINT)
               SO(INT1+4*NCCINT2)=FAC*AO(INT11+4*NCCINT)
               SO(INT1+5*NCCINT2)=FAC*AO(INT11+5*NCCINT)
900           CONTINUE
              ISOADR=ISOFF+NOABCD
              IAOADR=IAOFF+NOABCD
              ISOFF=ISOFF+NOABCD
800          CONTINUE
700         CONTINUE
600        CONTINUE
500       CONTINUE
240       CONTINUE
          IAOFF=IAOFF+NOABCD
400      CONTINUE
300     CONTINUE
200    CONTINUE
100   CONTINUE
      RETURN
      END
