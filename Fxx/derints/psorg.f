











C     ORGANIZE SORTING OF SYMMETRY-ADAPTED SECOND-ORDER REDUCED
C     DENSITY MATRIX TO MATCH CALCULATION OF DISTINCT TWO-ELECTRON
C     INTEGRAL DERIVATIVES

      SUBROUTINE PSORG(IWORK1,LWORK1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*80 FNAME
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
      PARAMETER ( MAXBUF2=10240 )
C
      LOGICAL BIGVEC, SEGMEN
      LOGICAL REDUCE1,REDUCE2
      DIMENSION IWORK1(LWORK1)
      COMMON /PTRFIL/ JOBIN, JOBOUT, LUMC, LUSCR, LUDA, LUPSO, LUPAO
      COMMON /CCOM/ THRS, NHTYP, IBX
      COMMON /BLOCKS/ CENTSH(MXSHEL,3),
     &                MAXSHL, BIGVEC, SEGMEN,
     &                NHKTSH(MXSHEL), KHKTSH(MXSHEL), MHKTSH(MXSHEL),
     &                ISTBSH(MXSHEL), NUCOSH(MXSHEL), NORBSH(MXSHEL),
     &                NSTRSH(MXSHEL), NCNTSH(MXSHEL), NSETSH(MXSHEL),
     &                JSTRSH(MXSHEL,MXAOVC),
     &                NPRIMS(MXSHEL,MXAOVC),
     &                NCONTS(MXSHEL,MXAOVC),
     &                IORBSH(MXSHEL,MXAOVC),
     &                IORBSB(MXCORB), NRCSH(MXSHEL)
      COMMON /PINCOM/ IPIND(MXCORB), IBLOCK(MXCORB), INDGEN(MXCORB)
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     &                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     &                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     &                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     &                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     &                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     &                IPTXYZ(3,0:7)
      COMMON /PTRTOT/ NPTOTR
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/IPRI/IPRINT
      COMMON/MVDINT/REDUCE1,REDUCE2
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
c 
c YAU - not used in this file
c#include "icdacc.com"
 
C
      DATA IONE /1/
C
      CALL TIMER(1)
C
C MAXIMUM AMOUNT OF CORE MEMORY AVAILABLE FOR SORT
C
      MXCOR=LWORK1*IINTFP
      IF(REDUCE2) MXCOR=MXCOR/2 
C
C ASSIGN UNITS TO THE FILES
C
C  LUDA : SORT FILE
C  LUPAO : SORTED GAMMA FILE
C  LUPSO : UNSORTED GAMMA FILE
C 
      LUDA=24
      LUPAO=25 
C
C LBUF :  BUFFER LENGTH FOR UNSORTED GAMMA FILE 
C         COMPRESSED AND EXPANDED ARRAY
C LBUF2 : BUFFER LENGTH FOR SORT FILE
C         VALUES, INDICES, NUMBER OF VALUES, NEXt RECORD
C
      CALL GETREC(20,'JOBARC','MAXAODSZ',IONE,LBUF)
      LBUF2=MAXBUF2+IINTFP*MAXBUF2+2
C
C  SUBTRACT THIS FROM THE SORTING AREA

      MXCOR = MXCOR-MAX(LBUF2,3*LBUF*IINTFP)
C
C  SUBTRACT ALSO MEMORY REQUIRED FOR THE INDEXING ARRAYS
C
      MXCOR = MXCOR-3*MAXSHL*(MAXSHL+1)/2-MAXSHL*16
C
C  OPEN THE SORTED GAMMA FILE
C
      CALL GFNAME('AOGAM   ',FNAME,ILENGTH)
      OPEN (UNIT=LUPAO,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     1      FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED')
C
C  I001 :   ARRAY IJADD  LENGTH MAXSHL*(MAXSHL+1)/2
      I001=1
C  I002 :   ARRAY NNIJ   LENGTH MAXSHL*(MAXSHL+1)/2
      I002=I001+MAXSHL*(MAXSHL+1)/2+1
C  I003 :   ARRAY ISHLD LENGTH MAXSHL*(MAXSHL+1)/2
      I003=I002+MAXSHL*(MAXSHL+1)/2
C  I004 :   ARRAY NSTSH LENGTH 8*MAXSHL
      I004=I003+MAXSHL*(MAXSHL+1)/2
C  I005 :   ARRAY NENDSH LENGTH 8*MAXSHL
      I005=I004+8*MAXSHL
C  I010 : OFFSET FOR BUFFERS
      I010=I005+8*MAXSHL
      I010=I010+MOD(I010+1,IINTFP)
C
C  STARTING ADDRESS FOR ARRAY LASTAD WHICH CONTAINS THE
C  NUMBER OF THE LAST RECORD FOR EACH BIN
C
      I030=I010+MAXBUF2
C
C  TOTAL AMOUNT OF MEMORY AVAILABLE FOR THE SORT
      LSORT=MXCOR/IINTFP
C
C  TOTAL AMOUNT OF MEMORY AVAILABLE FOR THE SORT
      LBUCK=MXCOR/(1+IINTFP)
C
C  SORT AREA OFFSET AND FOR THE BUCKETS
C
      I040=I010+MAX(3*LBUF*IINTFP,LBUF2)
      I040=I040+MOD(I040+1,IINTFP)
C
C  OFFSETS FOR THE PART OF THE BUCKETS HOLDING INDICES
C
      I050=I040+IINTFP*LBUCK
C
C  OFFSETS FOR THE EXPANDED ARRAY IN THE READ OF THE UNSORTED GAMMAS
C
       I060=I010+LBUF*IINTFP
C
      CALL GAMSRT(IWORK1(I010),IWORK1(I010),LBUF,IWORK1(I060),
     &          IWORK1(I010),IWORK1(I030),
     &          IWORK1(I040),LSORT,IWORK1(I040),IWORK1(I050),
     &          LBUCK,IWORK1(I001),IWORK1(I002),IWORK1(I003),
     &          IWORK1(I004),IWORK1(I005))
C
C  WRITE END FILE TO LUPAO
C
      END FILE LUPAO
C
C  DELETE SORT AND SCR FILE
C
      CLOSE(LUDA ,STATUS='DELETE')
C
      CALL TIMER(1)
cmn      write(*,9000) TIMENEW
9000  FORMAT(' Sort of Gammas required ',F7.2,' seconds.')
      RETURN
      END
