










      SUBROUTINE PTPRNT
C
C....    PRINT DATA READ FROM INPUT AND FROM INTERFACE, TOGETHER
C....    WITH DATA INITIALIZED IN COMMON
C
C
C....    ******      VERSION 1.0      RELEASE 820423      ******
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -------------------------------------------------------
C     IRAT  = (real word length) / (integer word length)
C     IRAT2 = (real word length) / (half-integer word length)
C             if available and used, otherwise IRAT2 = IRAT
      PARAMETER (IRAT = 2, IRAT2 = 2)

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

C
      COMMON /PTRFIL/ JOBIN, JOBOUT, LUMC, LUSCR, LUDA, LUPSO, LUPAO
      COMMON /PTRBUF/ MX1BUF, L1BUF, MEMS, MEMT, MAXCHN, MXABUF, LABUF,
     *                LDAMAX, LASTAD(200), MX2BUF, L2BUF,  MAXSCM,
     *                MAXLCM, IADR(200),   NCHAIN, MAXADR, NBLOCK
      COMMON /PTRORB/ NISH(8), NASH(8), NORB(8), NAC, NORBT, NAPRE(8),
     1                NORPRE(8), IASYM(MXCORB), IOSYM(MXCORB), LSMAX,
     2                MOFF(MXCORB), NORBM, NORBL1, NORBL2
      COMMON /PTRINF/ THR1, THR2, DAWEIG, SQWEIG, IPRINT
      CHARACTER*6 TITLE, RUNTP
      COMMON /PTRCHR/ TITLE(12), RUNTP
      COMMON /PTRSYM/ NSYM, MUL(8,8), NENT(8,8), ISYOFF(8,8)
C
C....    PRINT TITLE FRAME ONLY IF A TITLE WAS SPECIFIED
C
      IF (TITLE(1) .NE. '      '  .OR.
     1    TITLE(2) .NE. '      '  .OR.
     2    TITLE(3) .NE. '      ' ) THEN
         WRITE(JOBOUT,2000) (TITLE(I),I = 1,12)
      END IF
C
      MXAB = (LDAMAX - 2)/(IRAT+1)
      WRITE(JOBOUT,2050) RUNTP
      WRITE(JOBOUT,2100) JOBIN, JOBOUT, LUMC, LUSCR, LUDA, LUPSO
      WRITE(JOBOUT,2200) LUSCR, MX1BUF, LUDA, MXAB, LUPSO, MX2BUF
      WRITE(JOBOUT,2300) MAXCHN
      WRITE(JOBOUT,2400) MAXSCM, MAXLCM
      WRITE(JOBOUT,2500) SQWEIG, DAWEIG
      WRITE(JOBOUT,2600) THR1, THR2
      WRITE(JOBOUT,2700) IPRINT
      WRITE(JOBOUT,2800) NSYM
      DO 30 I = 1,NSYM
        WRITE(JOBOUT,2850) (MUL(I,J),J = 1,NSYM)
   30 CONTINUE
      WRITE(JOBOUT,2900) (I,I = 1,NSYM)
      WRITE(JOBOUT,2910) (NISH(I),I = 1,NSYM)
      WRITE(JOBOUT,2920) (NASH(I),I = 1,NSYM)
      WRITE(JOBOUT,2930) (NORB(I),I = 1,NSYM)
      IF (IPRINT .LT. 5) GOTO 20
      WRITE(JOBOUT,3000) (I,IASYM(I),I = 1,NAC)
      WRITE(JOBOUT,3100) (I,IOSYM(I),I = 1,NORBT)
      WRITE(JOBOUT,3200) (NAPRE(I),I = 1,NSYM)
      WRITE(JOBOUT,3300) (NORPRE(I),I = 1,NSYM)
   20 CONTINUE
      RETURN
 2000 FORMAT(//,1X,71('*'),/,' *',69X,'*',/,' *',2X,11A6,' *',
     1        /,' *',69X,'*',/,1X,71('*'))
 2050 FORMAT(//,'    TRANSFORMATION OF A ',A6,' 2-MATRIX')
 2100 FORMAT(//,'    UNIT NUMBERS:',
     1       //,'    FUNCTION',23X,'UNIT',
     2       //,4X,'INPUT',26X,I4,/,
     3       4X,'OUTPUT',25X,I4,/,
     4       4X,'INTERFACE',22X,I4,/,
     5       4X,'HALF-TRANSFORMED 2-MATRIX',6X,I4,/,
     6       4X,'DIRECT ACCESS SORT',13X,I4,/,
     7       4X,'TRANSFORMED 2-MATRIX',11X,I4)
 2200 FORMAT(//,'    NUMBER OF 2-MATRIX ELEMENTS PER BUFFER:',//,
     2        6X,I5,11X,I5,/,
     3        6X,I5,11X,I5,' (DEFAULT MAXIMUM ONLY)',/,
     4        6X,I5,11X,I5)
 2300 FORMAT(//,'    MAXIMUM NUMBER OF BIN SORT CHAINS',I5)
 2400 FORMAT(//,'    WORKING STORAGE:',//,
     1        4X,I8,' WORDS SCM AVAILABLE',/,
     2        4X,I8,' WORDS LCM AVAILABLE')
 2500 FORMAT(//,'    I/O WEIGHT FACTORS ',
     1          '(EXPRESSED RELATIVE TO DA WRITE):',//,
     2        6X,'SEQUENTIAL READ',2X,F10.4,/,
     3        6X,'DA READ',10X,F10.4)
 2600 FORMAT(//,'    THRESHOLDS FOR 2-MATRIX ELEMENTS:',//,
     1        6X,'DURING FIRST  HALF-TRANSFORMATION',2X,E16.8,/,
     2        6X,'DURING SECOND HALF-TRANSFORMATION',2X,E16.8)
 2700 FORMAT(//,'    PRINT LEVEL IS',I5)
 2800 FORMAT(//,'    SYMMETRY DATA:',//,
     1        6X,'ORDER OF GROUP',I4,//,
     2        10X,'MULTIPLICATION TABLE',/)
 2850 FORMAT(15X,8I4)
 2900 FORMAT(//,'    ORBITAL DATA:',//,4X,'IRREP',15X,8I5)
 2910 FORMAT(4X,'INACTIVE MOS',8X,8I5)
 2920 FORMAT(4X,'ACTIVE MOS',10X,8I5)
 2930 FORMAT(4X,'SYMMETRY ORBITALS',3X,8I5)
 3000 FORMAT(//,'    ACTIVE ORBITAL SYMMETRIES:',//,
     1        10X,'ORBITAL',4X,'SYMMETRY',//,
     2        (10X,I5,6X,I5))
 3100 FORMAT(//,'    SYMMETRY ORBITALS:',//,
     1        10X,'ORBITAL',4X,'SYMMETRY',//,
     2        (10X,I5,6X,I5))
 3200 FORMAT(//,'    ACTIVE ORBITAL OFFSET VECTOR:',//,4X,8I5)
 3300 FORMAT(//,'    SYMMETRY ORBITAL OFFSET VECTOR',/,4X,8I5)
      END
