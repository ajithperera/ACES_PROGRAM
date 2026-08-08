










      SUBROUTINE DSYM2(DHFSO,DRELSO,DNHFSO,NBAST,IPRINT,
     &                 DHFA,DHFB,DRELA,DRELB,DNHFA,DNHFB,
     &                 NIR)
C
C     Take total and active density matrices in symmetry
C     orbital basis and generate density matrices (non-folded)
C     over distinct pairs of AOs
C
CEND
C                                          880418  PRT
C
C  ADAPTED TO THE ACES II ENVIRONMENT AND EXTENDED FOR UHF AND NON-HF
C  OCT/90/JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
      LOGICAL SHARE
c&line modif
      LOGICAL SCF,NONHF
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
      DIMENSION DHFSO(1),DRELSO(1),DNHFSO(1)
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
      DIMENSION DHFA(NBAST,NBAST,0:NIR),DHFB(NBAST,NBAST,0:NIR),
     &              DRELA(NBAST,NBAST,0:NIR),DRELB(NBAST,NBAST,0:NIR),
     &              DNHFA(NBAST,NBAST,0:NIR),DNHFB(NBAST,NBAST,0:NIR)
      COMMON/METHOD/IUHF,SCF,NONHF
c&line del      COMMON/QRHF/QRHFP
c&new lines
      INTEGER POPRHF,VRTRHF,POPDOC,VRTDOC
      COMMON/QRHFINF/POPRHF(8),VRTRHF(8),NOSH1(8),NOSH2(8),
     &               POPDOC(8),VRTDOC(8),NAI,N1I,N2A,
     &               NUMISCF,NUMASCF,ISPINP,ISPINM,IQRHF
c&new lines end
      COMMON/BASSYM/NBAS(8),NBASIS,NBASSQ,NBASTT
C
      DATA ONE,TWO /1.0D0,2.D0/
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Subroutine DSYM2',-1)
      END IF
C
      IF(IUHF.EQ.0) THEN
       CALL ZERO(DHFA(1,1,0),NBAST*NBAST*(NIR+1))
       IF(.NOT.SCF) THEN
        CALL ZERO(DRELA(1,1,0),NBAST*NBAST*(NIR+1))
       ENDIF
      ELSE
       CALL ZERO(DHFA(1,1,0),NBAST*NBAST*(NIR+1))
       CALL ZERO(DHFB(1,1,0),NBAST*NBAST*(NIR+1))
       IF(.NOT.SCF) THEN
        IF(NONHF) THEN
         CALL ZERO(DNHFA(1,1,0),NBAST*NBAST*(NIR+1))
         CALL ZERO(DNHFB(1,1,0),NBAST*NBAST*(NIR+1))
        ENDIF
        CALL ZERO(DRELA(1,1,0),NBAST*NBAST*(NIR+1))
        CALL ZERO(DRELB(1,1,0),NBAST*NBAST*(NIR+1))
       ENDIF
      ENDIF
C
      IOFFA=0
      IOFFB=IOFFA+NBASTT*IUHF
C
C     Loop over all irreps in molecule
C
      ISOFF = 0
      ISTR = 1
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
               IB    = IBTSHR(IPIND(J),16)
               NB    = IBTSHR(IPIND(J),8) - IB*256
               JOFF  = KSTRT(IB)
               NHKTB = NHKT(IB)
               KHKTB = KHKT(IB)
               MULB  = ISTBAO(IB)
               MAB   = IBTOR(MULA,MULB)
               ISOFF = ISOFF + 1
               DHFAIJ = DHFSO(ISOFF)
               DHFBIJ = DHFSO(ISOFF+IOFFB)
               DRELAIJ=DRELSO(ISOFF)
               DRELBIJ=DRELSO(ISOFF+IOFFB)
               DNHFAIJ=DNHFSO(ISOFF)
               DNHFBIJ=DNHFSO(ISOFF+IOFFB)
               DO 400 ISYMOP = 0, MAXLOT
                  FAC = PT(IBTAND(ISYMOP,IREP))
                  INDB = JOFF + NB
                  DHFA(INDA,INDB,ISYMOP) = DHFA(INDA,INDB,ISYMOP)
     *                                     + FAC*DHFAIJ
                  DHFA(INDB,INDA,ISYMOP) = DHFA(INDA,INDB,ISYMOP)
                IF(.NOT.SCF) THEN
                  DRELA(INDA,INDB,ISYMOP) = DRELA(INDA,INDB,ISYMOP)
     *                                     + FAC*DRELAIJ
                  DRELA(INDB,INDA,ISYMOP) = DRELA(INDA,INDB,ISYMOP)
                 IF(NONHF) THEN
                  DNHFA(INDA,INDB,ISYMOP) = DNHFA(INDA,INDB,ISYMOP)
     *                                     + FAC*DNHFAIJ
                  DNHFA(INDB,INDA,ISYMOP) = DNHFA(INDA,INDB,ISYMOP)
                 ENDIF
                ENDIF
                IF(IUHF.NE.0) THEN
                  DHFB(INDA,INDB,ISYMOP) = DHFB(INDA,INDB,ISYMOP)
     *                                     + FAC*DHFBIJ
                  DHFB(INDB,INDA,ISYMOP) = DHFB(INDA,INDB,ISYMOP)
                 IF(.NOT.SCF) THEN
                  DRELB(INDA,INDB,ISYMOP) = DRELB(INDA,INDB,ISYMOP)
     *                                     + FAC*DRELBIJ
                  DRELB(INDB,INDA,ISYMOP) = DRELB(INDA,INDB,ISYMOP)
                  IF(NONHF) THEN
                  DNHFB(INDA,INDB,ISYMOP) = DNHFB(INDA,INDB,ISYMOP)
     *                                     + FAC*DNHFBIJ
                  DNHFB(INDB,INDA,ISYMOP) = DNHFB(INDA,INDB,ISYMOP)
                  ENDIF
                 ENDIF
                ENDIF
400            CONTINUE
300         CONTINUE
200      CONTINUE
110      CONTINUE
         ISTR = ISTR + NORBI
100   CONTINUE
CSSS      IPRINT = 11
      IF (IPRINT .GT. 10) THEN
         DO 500 ISYMOP = 0, MAXLOP
         CALL HEADER
     *        ('Total density matrix (symmetry distinct AO basis)',-1)
         CALL OUTPUT(DHFA(1,1,ISYMOP),1,NBAST,1,NBAST,
     *        NBAST,NBAST,1,LUPRI)
         CALL OUTPUT(DHFB(1,1,ISYMOP),1,NBAST,1,NBAST,
     *        NBAST,NBAST,1,LUPRI)
         IF(.NOT.SCF) THEN
         CALL HEADER
     *        ('Active density matrix (symmetry distinct AO basis)',-1)
          CALL OUTPUT(DRELA(1,1,ISYMOP),1,NBAST,1,NBAST,
     *                NBAST,NBAST,1,LUPRI)
          CALL OUTPUT(DRELB(1,1,ISYMOP),1,NBAST,1,NBAST,
     *                NBAST,NBAST,1,LUPRI)
          IF(NONHF) THEN
           CALL OUTPUT(DNHFA(1,1,ISYMOP),1,NBAST,1,NBAST,
     *                 NBAST,NBAST,1,LUPRI)
           CALL OUTPUT(DNHFB(1,1,ISYMOP),1,NBAST,1,NBAST,
     *                 NBAST,NBAST,1,LUPRI)
          ENDIF
         ENDIF
  500    CONTINUE
      END IF
      IPRINT = 0
C
C IN THE NONHF CASE ADD ALPHA AND BETA PART OF NONHF DENSITIES
C AND STORE REFERENCE DENSITY IN THE BETA BLOCK
C
      IF(NONHF) THEN 
       CALL SAXPY(NBAST*NBAST*(NIR+1),ONE,
     &            DNHFB(1,1,0),1,DNHFA(1,1,0),1)
c&line modif
       IF(IQRHF.EQ.1) THEN
        CALL SCOPY(NBAST*NBAST*(NIR+1),DHFA(1,1,0),1,
     &             DNHFB(1,1,0),1)
        CALL SSCAL(NBAST*NBAST*(NIR+1),TWO,DNHFB(1,1,0),1)
       ELSE
        CALL SCOPY(NBAST*NBAST*(NIR+1),DHFB(1,1,0),1,
     &             DNHFB(1,1,0),1)
        CALL SSCAL(NBAST*NBAST*(NIR+1),TWO,DNHFB(1,1,0),1)
       ENDIF
      ENDIF
      RETURN
      END
