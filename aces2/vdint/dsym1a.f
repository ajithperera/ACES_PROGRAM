










      SUBROUTINE DSYM1A(DEN,DSOA,DSOB,NBAST,IPRINT)
C
C     Take density matrix in symmetry orbital basis and generate
C     density matrix over distinct pairs of AOs
C
CEND
C                                          880418  PRT
C  
C  ADAPTED TO HE ACES II ENVIRONEMENT   OCT/90/JG 
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


c IMPORTANT: The ECP library contains three routines (doitgr, ecpabl, ecpder)
c that use this parameter but do not include this file. Changing this value
c requires updating those files as well. [Yau: If we linked them together,
c then compiling ecp would require the vdint directory. Alternatively, we
c could put the mxcent.par file into the top-level include directory.]

C MXCENT : Maximum number of atoms currently allowed
C MXCOOR : The number of Cartesian that correspond to MXCENT

      INTEGER MXCENT, MXCOOR
      PARAMETER (MXCENT=200, MXCOOR=3*MXCENT)
      DIMENSION DSOA(1), DSOB(1),DEN(1)
      LOGICAL SHARE
      LOGICAL SCF,NONHF
      COMMON /SHELLSi/ KMAX,
     &                NHKT(MXSHEL),   KHKT(MXSHEL), MHKT(MXSHEL),
     &                ISTBAO(MXSHEL), NUCO(MXSHEL), JSTRT(MXSHEL),
     &                NSTRT(MXSHEL),  MST(MXSHEL),  NCENT(MXSHEL),
     &                NRCO(MXSHEL), NUMCF(MXSHEL),
     &                NBCH(MXSHEL),   KSTRT(MXSHEL)
      COMMON /SHELLS/ CENT(MXSHEL,3), SHARE(MXSHEL)
      COMMON /PINCOM/ IPIND(MXCORB), IBLOCK(MXCORB), INDGEN(MXCORB)
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     *                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     *                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     *                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     *                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     *                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     *                IPTXYZ(3,0:7)
      COMMON/METHOD/IUHF,SCF,NONHF
      COMMON/BASSYM/NBAS(8),NBASIS,NBASSQ,NBASTT
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Subroutine DSYM1',-1)
      END IF
C
C     Loop over all irreps in molecule
C
      ISOFF = 0
      ISTR = 1
      NNBASX = NBAST*(NBAST + 1)/2
      CALL ZERO(DSOA,NNBASX)
      IOFFB=NBASTT*IUHF
      IF(IUHF.NE.0) CALL ZERO(DSOB,NNBASX)
      DO 100 IREP = 0, MAXLOP
         NORBI = NAOS(IREP+1)
         IF (NORBI .EQ. 0) GOTO 110
         DO 200 I = ISTR,ISTR + NORBI - 1
            IA   = IBTSHR(IPIND(I),16)
            NA   = IBTSHR(IPIND(I),8) - IA*256
            IOFF = KSTRT(IA)
            MULA = ISTBAO(IA)
            INDA = IOFF + NA
            DO 300 J = ISTR,I
               IB     = IBTSHR(IPIND(J),16)
               NB     = IBTSHR(IPIND(J),8) - IB*256
               JOFF   = KSTRT(IB)
               NHKTB  = NHKT(IB)
               KHKTB  = KHKT(IB)
               NRCB   = NRCO(IB)
               NBTYP  = (NB-1)/NRCB+1
               MULB   = ISTBAO(IB)
               MAB    = IBTOR(MULA,MULB)
               KAB    = IBTAND(MULA,MULB)
               HKAB   = FMULT(KAB)
               ISOFF  = ISOFF + 1
               DSYMIJA = DEN(ISOFF)
               DSYMIJB = DEN(ISOFF+IOFFB)
               INDB   = JOFF + NB - KHKTB*NRCB
               DO 400 ISYMOP = 0, MAXLOT
                  IF (IBTAND(ISYMOP,MAB) .NE. 0) GOTO 400
                  INDB = INDB + KHKTB*NRCB
C
C                 Weight and parity factor
C
                  FAC = HKAB*
     *                   PT(IBTAND(ISYMOP,IBTXOR(IREP,
     *                       ITYPE(NHKTB,NBTYP))))
                  INDM = MAX(INDA,INDB)
                  IND  = (INDM*(INDM - 3))/2 + INDA + INDB
                  DSOA(IND) = DSOA(IND) + FAC*DSYMIJA
                  IF(IUHF.NE.0) DSOB(IND) = DSOB(IND) + FAC*DSYMIJB
400            CONTINUE
300         CONTINUE
200      CONTINUE
110      CONTINUE
         ISTR = ISTR + NORBI
100   CONTINUE
      IF (IPRINT .GT. 10) THEN
         CALL PRITRI(DSOA,NBAST,
     *    'Folded alpha density matrix (symmetry distinct AO basis)')
        IF(IUHF.NE.0) THEN
         CALL PRITRI(DSOB,NBAST,
     *    'Folded beta density matrix (symmetry distinct AO basis)')
       ENDIF
      ENDIF
      RETURN
      END
