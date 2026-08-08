










      SUBROUTINE SYM1N(AO,SO,IREPO,KB,MULA,MULB,NHKTA,NHKTB,
     &                 KHKTA,KHKTB,HKAB,LDIAG,LAEQB,FULMAT,
     &                 THRESH,IMAT0,JSKIP,IPRINT,ASYM)
C
C  Take block of distinct AO integral (derivatives) and
C  generate symmetrized contributions to SO integral
C  (derivatives) over non-symmetric operators
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
      LOGICAL LDIAG, FULMAT, LAEQB
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
      PARAMETER (LUTEMP = 48)
C
      DIMENSION AO(*), SO(*)
      COMMON /CSYM1/ BUF(600), IBUF(600), LENGTH, INDMAX
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     &                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     &                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     &                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     &                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     &                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     &                IPTXYZ(3,0:7)
      COMMON /SYMIND/ INDFA(8,10), INDFB(8,10), ISOFRA(8), ISOFRB(8)
      COMMON /GENCON/ NRCA,NRCB,CONTA(MXCONT), CONTB(MXCONT)
C
      IBTAND(I,J) = AND(I,J)
      IBTXOR(I,J) = XOR(I,J)
C
      IF (IPRINT .GT. 10) THEN
       CALL HEADER('Subroutine SYM1N',-1)
       WRITE (LUPRI, '(A,2I5)') ' IREPO, KB ', IREPO, KB
       WRITE (LUPRI, '(A,2I5)') ' NHKTA/B ', NHKTA, NHKTB
      END IF
C
C  Loop over irreps for first basis function - those for second
C  are obtained from operator symmetry IREPO
C
      DO 100 IREPA = 0, MAXLOP
       IREPB = IBTXOR(IREPO,IREPA)
       IF (FULMAT) INDOFF = NPARNU(IREPO+1,MAX(IREPA,IREPB)+1)
C
C  Loop over AOs which are of symmetry IREPA in stabilizer MULA
C
       ISKIP=NRCA*NRCB
       IF(LDIAG) ISKIP=NRCA*(NRCA+1)/2
C
       DO 200 NA = 1, KHKTA
        IF (IBTAND(MULA,IBTXOR(IREPA,ITYPE(NHKTA,NA))).EQ.0) THEN
         NAT = KHKTB*(NA - 1)
         DO 300 NB = 1,KHKTB
          IF (IBTAND(MULB,IBTXOR(IREPB,ITYPE(NHKTB,NB))).EQ.0) THEN
C
C  Weight and parity factor
C
           FAC = HKAB*PT(IBTAND(KB,IBTXOR(IREPB,ITYPE(NHKTB,NB))))
C
C  Locate SO integrals to which AOs contribute
C
           INT=0
           DO 160 IRCA=1,NRCA
            MAXB=NRCB
            IF(LDIAG) MAXB=IRCA
            DO 150 IRCB=1,MAXB
             INT=INT+1
             IF(LDIAG.AND.(IRCA.EQ.IRCB).AND.NA.LT.NB) GO TO 150
             IF(LAEQB.AND.NA.EQ.NB.AND.IRCA.EQ.IRCB.AND.
     &          IREPA.LT.IREPB) GO TO  150
             IF(LAEQB.AND.IRCA.LT.IRCB) GO TO 150
             IF(LAEQB.AND.(IRCA.EQ.IRCB).AND.NA.LT.NB) GO TO 150
             INDA = INDFA(IREPA + 1,NA)-NRCA+IRCA
             INDB = INDFB(IREPB + 1,NB)-NRCB+IRCB
             RINT = FAC*AO((INT+(NAT+NB-1)*ISKIP-1)*JSKIP+1)
             IF (FULMAT) THEN
              IF (IREPA .GE. IREPB) THEN
               IND  = INDOFF + NAOS(IREPB+1)*(INDA-1) + INDB
               RINT=RINT*ASYM
              ELSE
               IND  = INDOFF + NAOS(IREPA+1)*(INDB-1) + INDA
              ENDIF
              SO(IND) = SO(IND) + RINT
             ELSE
              IF (ABS(RINT) .GT. THRESH) THEN
               IF (IREPA .GE. IREPB) THEN
                IND  = NAOS(IREPB+1)*(INDA-1) + INDB
               ELSE
                IND  = NAOS(IREPA+1)*(INDB-1) + INDA
               ENDIF
               INDMAX = MAX(IND,INDMAX)
               LABEL  = IND*2**16 + IMAT0 + MAX(IREPA,IREPB)
               IF (IPRINT .GT. 20) THEN
                WRITE (LUPRI,'(A,F12.6,2I3,2I2,I4,I2,I5)')
     &                 'SYM1N - NA/B,IREPA/B,IND,IREPO,IMAT',
     &                 RINT, NA, NB, IREPA, IREPB, IND, IREPO,
     &                 IMAT0 + MAX(IREPA,IREPB)
               END IF
               IF (LENGTH .EQ. 600) THEN
                WRITE (LUTEMP) BUF, IBUF, LENGTH
                IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A,I4,A)')
     &                             ' Buffer of length',LENGTH,
     &                             ' has been written in SYM1N.'
                LENGTH = 0
               ENDIF
               LENGTH = LENGTH + 1
               BUF(LENGTH)  = RINT
               IBUF(LENGTH) = LABEL
              ENDIF
             ENDIF
150         CONTINUE
160        CONTINUE
          END IF
300      CONTINUE
        END IF
200    CONTINUE
100   CONTINUE
      RETURN
      END
