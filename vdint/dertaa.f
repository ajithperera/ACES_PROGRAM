










      SUBROUTINE DERTAA(MAXA,MAXB,ISTEPA,ISTEPB,NUCAB,IOFF,IOFF0,ITYPE,
     *                  IAB0,WORD,IPRINT)
C
C     TUH 84
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
C     -------------------------------------------------------
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

      DIMENSION STEXP1(MXAOSQ)
      CHARACTER WORD*4

c Dimensions of the WORK2 and WORK3 arrays of the VDINT AME. For
c machines that have small stack space, these values must be small.
c Marshall Cory, who has access to CRAY machines with large memory, uses
c 900000 and 500000 is recommended for the others. Small values for
c these arrays will crash the program for basis sets that have large
c primitive spaces. For example the ANO basis sets. A. Perera, 03/2005. 

      PARAMETER (LWORK2 = 300 0000)
      PARAMETER (LWORK3 = 900 0000)
      COMMON /CWORK2/ WK2LOW, WORK2(LWORK2), WK2HGH
      COMMON /TWOVEC/ DIFPAX(MXAOSQ), DIFPAY(MXAOSQ), DIFPAZ(MXAOSQ),
     *                DIFPBX(MXAOSQ), DIFPBY(MXAOSQ), DIFPBZ(MXAOSQ),
     *                TEXP1(MXAOSQ),  TEXP2(MXAOSQ),  HEXPPI(MXAOSQ)
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
      IF (ITYPE .EQ. 1) THEN
         MAX1 = MAXA
         MAX2 = MAXB
         ISTEP1 = ISTEPA
         ISTEP2 = ISTEPB
         IADDTE = 0
      ELSE
         MAX1 = MAXB
         MAX2 = MAXA
         ISTEP1 = ISTEPB
         ISTEP2 = ISTEPA
         IADDTE = MXAOSQ
      END IF
      DO 90 I = 1, NUCAB
         TEXP = TEXP1(IADDTE + I)
         STEXP1(I) = TEXP*TEXP
   90 CONTINUE
      IAB0P = IAB0 + 1
      NUCADD = IAB0P*NUCAB
      IADD0 = - IOFF + IOFF0
      IADDM = IADD0 - ISTEP1 - ISTEP1
      IADDP = IADD0 + ISTEP1 + ISTEP1
C
C     *******************************
C     ********** I1 = 0, 1 **********
C     *******************************
C
      ISTRT1 = IOFF
      DO 100 I1 = 0, MIN(MAX1,1)
         FAC0 = - FLOAT(I1 + I1 + 1)
         ISTR12 = ISTRT1
         DO 110 I2 = 0, MAX2
            I12 = I1 + I2
            MINT = IBTAND(I12,IAB0)
            IADR0 = ISTR12 + MINT*NUCAB
C
C           ***** IT = 0, I12 *****
C
            DO 200 IT = MINT, I12, IAB0P
               IODC00 = IADR0 + IADD0
               IODCP0 = IADR0 + IADDP
               DO 210 I = 1, NUCAB
                  TEXP = TEXP1(IADDTE + I)
                  WORK2(IADR0 + I) = TEXP*(TEXP*WORK2(IODCP0 + I)
     *                                  + FAC0*WORK2(IODC00 + I))
  210          CONTINUE
               IADR0 = IADR0 + NUCADD
  200       CONTINUE
C
C           ***** IT = I12 + 1, I12 + 2 *****
C
            IF (IAB0 .EQ. 0) THEN
               IADR1 = IADR0 + NUCADD
               IODCP0 = IADR0 + IADDP
               IODCP1 = IADR1 + IADDP
               DO 220 I = 1, NUCAB
                  STEXP = STEXP1(I)
                  WORK2(IADR0 + I) = STEXP*WORK2(IODCP0 + I)
                  WORK2(IADR1 + I) = STEXP*WORK2(IODCP1 + I)
  220          CONTINUE
            ELSE
               IODCP0 = IADR0 + IADDP
               DO 230 I = 1, NUCAB
                  WORK2(IADR0 + I) = STEXP1(I)*WORK2(IODCP0 + I)
  230          CONTINUE
            END IF
            ISTR12 = ISTR12 + ISTEP2
  110    CONTINUE
         ISTRT1 = ISTRT1 + ISTEP1
  100 CONTINUE
C
C     ****************************
C     ********** I1 > 1 **********
C     ****************************
C
      DO 300 I1 = 2, MAX1
         FACM = FLOAT(I1*(I1 - 1))
         FAC0 = - FLOAT(I1 + I1 + 1)
         ISTR12 = ISTRT1
         DO 310 I2 = 0, MAX2
            I12 = I1 + I2
            MINT = IBTAND(I12,IAB0)
            IADR0 = ISTR12 + MINT*NUCAB
C
C           ***** IT =  0, I12 - 2 *****
C
            DO 420 IT = MINT, I12 - 2, IAB0P
               IODCM0 = IADR0 + IADDM
               IODC00 = IADR0 + IADD0
               IODCP0 = IADR0 + IADDP
               DO 430 I = 1, NUCAB
                  TEXP = TEXP1(IADDTE + I)
                  WORK2(IADR0 + I) = FACM*WORK2(IODCM0 + I)
     *                             + TEXP*(TEXP*WORK2(IODCP0 + I)
     *                             + FAC0*WORK2(IODC00 + I))
  430          CONTINUE
               IADR0 = IADR0 + NUCADD
  420       CONTINUE
C
C           ***** IT = I12 - 1, I12 + 2 *****
C
            IF (IAB0 .EQ. 0) THEN
               IADR1 = IADR0 + NUCADD
               IADR2 = IADR1 + NUCADD
               IADR3 = IADR2 + NUCADD
               IODC00 = IADR0 + IADD0
               IODC01 = IADR1 + IADD0
               IODCP0 = IADR0 + IADDP
               IODCP1 = IADR1 + IADDP
               IODCP2 = IADR2 + IADDP
               IODCP3 = IADR3 + IADDP
               DO 540 I = 1, NUCAB
                  TEXP = TEXP1(IADDTE + I)
                  STEXP = TEXP*TEXP
                  TFAC0 = TEXP*FAC0
                  WORK2(IADR0 + I) = STEXP*WORK2(IODCP0 + I)
     *                             + TFAC0*WORK2(IODC00 + I)
                  WORK2(IADR1 + I) = STEXP*WORK2(IODCP1 + I)
     *                             + TFAC0*WORK2(IODC01 + I)
                  WORK2(IADR2 + I) = STEXP*WORK2(IODCP2 + I)
                  WORK2(IADR3 + I) = STEXP*WORK2(IODCP3 + I)
  540          CONTINUE
            ELSE
               IADR1 = IADR0 + NUCADD
               IODC00 = IADR0 + IADD0
               IODCP0 = IADR0 + IADDP
               IODCP1 = IADR1 + IADDP
               DO 550 I = 1, NUCAB
                  TEXP = TEXP1(IADDTE + I)
                  WORK2(IADR0 + I)  = TEXP*(TEXP*WORK2(IODCP0 + I)
     *                                    + FAC0*WORK2(IODC00 + I))
                  WORK2(IADR1 + I) = TEXP*TEXP*WORK2(IODCP1 + I)
  550          CONTINUE
            END IF
            ISTR12 = ISTR12 + ISTEP2
  310    CONTINUE
         ISTRT1 = ISTRT1 + ISTEP1
  300 CONTINUE
      IF (IPRINT .LT. 10) RETURN
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      WRITE (LUPRI, 2000)
      WRITE (LUPRI, 2010) MAX1, MAX2
      WRITE (LUPRI, 2020) ISTEPA, ISTEPB
      WRITE (LUPRI, 2030) NUCAB
      WRITE (LUPRI, 2040) IOFF0
      WRITE (LUPRI, 2050) IOFF
      WRITE (LUPRI, 2060) ITYPE
      IF (IPRINT .LT. 20) RETURN
      ISTR12 = IOFF
      DO 1000 IA = 0, MAXA
         IADR0 = ISTR12
         DO 1100 IB = 0, MAXB
            IADRT = IADR0
            DO 1200 IT = 0, IA + IB + 2
               IODD = IBTAND(IA + IB + 2 - IT,IAB0)
               IF (IODD .EQ. 0) THEN
                  WRITE (LUPRI, 2070) WORD, IA, IB, IT
                  WRITE (LUPRI, 2080) (WORK2(IADRT + I), I = 1, NUCAB)
               END IF
               IADRT = IADRT + NUCAB
 1200       CONTINUE
            IADR0 = IADR0 + ISTEPB
 1100    CONTINUE
         ISTR12 = ISTR12 + ISTEPA
 1000 CONTINUE
      RETURN
 2000 FORMAT (//,'  <<<<<<<<<< SUBROUTINE DERTAA >>>>>>>>>> ',/)
 2010 FORMAT ('  MAXA/B:   ',2I7)
 2020 FORMAT ('  ISTEPA/B: ',2I7)
 2030 FORMAT ('  NUCAB:    ',I7)
 2040 FORMAT ('  IOFF0:    ',I7)
 2050 FORMAT ('  IOFF:     ',I7)
 2060 FORMAT ('  ITYPE:    ',I7)
 2070 FORMAT (/,1X,A4,'(',I1,',',I1,';',I1,')',/)
 2080 FORMAT (1X,6F12.8)
      END
