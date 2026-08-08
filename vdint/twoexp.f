











C ROUTINE WHICH CONTROLS THE EVALUATION OF VARIOUS DERIVATIVES
C OF THE TWO-ELECTRON INTEGRALS
C
C  ADAPETD TO THE CRAPS ENVIRONMENT AND EXTENDED FOR UHF AND NON-HF
C  OCT/90/JG
C
C  EXTENDED FOR GIAO CALCULATIONS ON MAGNETIC PROPERTIES, SEP/91/JG  

       SUBROUTINE TWOEXP(WORK1,LWORK1,ENEREE,SCF_GRADEE,GRADEE,
     &                   SCF_HESSEE,HESSEE,MSZ,NCOORD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SCF,NONHF
      LOGICAL ELECT,GEOM,MAGNET,MAGNET2,GIAO
      LOGICAL         SKIP,   RUNPTR, RUNSOR, RUNINT, DTEST, TKTIME,
     *                RETUR,  NODC,   NODV,   NOPV
      LOGICAL MOLGRD,MOLHES,DIPDER,POLAR,INPTES,VIB,RESTAR,
     &        DOWALK,GDALL,FOCK,H2MO
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
      PARAMETER (LUPAO = 25)

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
      DIMENSION WORK1(LWORK1),GRADEE(NCOORD),HESSEE(NCOORD,NCOORD),
     &          MSZ(3,3),SCF_GRADEE(NCOORD),SCF_HESSEE(NCOORD,NCOORD)
      COMMON /CBITWO/ SKIP,   RUNPTR, RUNSOR, RUNINT, DTEST, TKTIME,
     &                IPRALL, IPRPTR, IPRSOR,
     &                JPRINT, IPRNTA, IPRNTB, IPRNTC, IPRNTD, RETUR,
     &                NODC,   NODV,   NOPV,
     &                MAXDIF, MAXVEC
      COMMON/ABAINF/IPRDEF,MOLGRD,MOLHES,DIPDER,POLAR,INPTES,
     &              VIB,RESTAR,DOWALK,GDALL,FOCK,H2MO
      COMMON/SYMMET/FMULT(0:7),PT(0:7),MAXLOP,MAXLOT,
     &              MULT(0:7),ISYTYP(3),
     &              ITYPE(8,36),NPARSU(8),NPAR(8),NAOS(8),
     &              NPARNU(8,8),IPTSYM(MXCORB,0:7),
     &              IPTCNT(3*MXCENT,0:7),NCRREO(0:7),
     &              IPTCOR(MXCENT*3),NAXREP(0:7),IPTAX(3),
     &              IPTXYZ(3,0:7)
C
C ADDED COMMON BLOCKS IN THE CRAPS ENVIRONMENT
C
      COMMON/PROP/ELECT,MAGNET,GEOM,MAGNET2,GIAO
      COMMON/METHOD/IUHF,SCF,NONHF
      COMMON/BASSYM/NBAS(8),NBASIS,NBASSQ,NBASTT
      COMMON/IPRI/IPRINT
 
c YAU - not used in this file
c#include "icdacc.com"
 
C
      DATA TWO /2.D0/ 
C
C WE HAVE OUR OWN BACK TRANSFORMATION SO WE DON'T NEED PTRAN
C
      RUNPTR = .FALSE.
      NODV=.TRUE.
      NOPV=.TRUE.
C
      IF(SCF.OR.FOCK) THEN
C
C  SCF : NO SORT REQUIRED
C
       RUNSOR = .FALSE.
C
      ELSE
C
C  CORRELATED METHODS : SORT IS REQUIRED FOR GRADIENTS
C  BUT NOT FOR CHEMICAL SHIFTS
C
       IF(.NOT.MAGNET) THEN
        RUNSOR=.TRUE.
       ELSE
        RUNSOR=.FALSE.
       ENDIF
C
      ENDIF
C
C  SET UP COMMON /BLOCKS/ FOR PSORT AND TWOINT 
C
      CALL PAOVEC(MAXVEC,0,IPRALL)
C
C  SET UP ONE-ELECTRON DENSITY MATRICES 
C
      IF (RUNINT) THEN
C
C  CORE MEMORY ALLOCATION IN THE CRAPS ENVIRONMENT
C
       KHFA=1
       KEND=KHFA+NBASIS*NBASIS*(MAXLOP+1)
       IF(IUHF.NE.0) THEN
        KHFB=KEND
        KEND=KHFB+NBASIS*NBASIS*(MAXLOP+1)
       ENDIF
       IF(.NOT.SCF) THEN
        KRELA=KEND
        KEND=KRELA+NBASIS*NBASIS*(MAXLOP+1)
        IF(IUHF.NE.0) THEN
         KRELB=KEND
         KEND=KRELB+NBASIS*NBASIS*(MAXLOP+1)
         IF(NONHF) THEN
          KNHFA=KEND
          KNHFB=KNHFA+NBASIS*NBASIS*(MAXLOP+1)
          KEND=KNHFB+NBASIS*NBASIS*(MAXLOP+1)
         ENDIF
        ENDIF
       ENDIF
       IF(FOCK) THEN
        KHFFA=KEND
        KHFFB=KHFFA+NBASIS*(NBASIS+1)/2*IUHF*(MAXLOP+1)
        KEND=KHFFB+NBASTT*(NBASIS+1)/2*(MAXLOP+1)
        IF(.NOT.SCF) THEN
         KRELFA=KEND
         KRELFB=KRELFA+NBASIS*(NBASIS+1)/2*(MAXLOP+1)*IUHF
         KEND=KRELFB+NBASIS*(NBASIS+1)/2*(MAXLOP+1)
        ENDIF
       ENDIF
       KCMO=KEND
       KDV=KCMO+NBASIS*NBASIS
       KDNH=KDV+NBASIS*NBASIS
       KSCR=KDNH+NBASIS*NBASIS
       KHFSO=KSCR+2*NBASIS*NBASIS
       KRELSO=KHFSO+NBASTT*(IUHF+1)
       KNHFSO=KRELSO+NBASTT*(IUHF+1)
       KEND2=KNHFSO+NBASTT*(IUHF+1)
       IF(KEND2.GE.LWORK1) CALL ERREX
C
C CONSTRUCT HERE DENSITY MATRICES IN SO BASIS
C
       CALL ONEDSO(WORK1(KCMO),WORK1(KDV),WORK1(KDNH),WORK1(KSCR),
     &             WORK1(KHFSO),WORK1(KRELSO),WORK1(KNHFSO))
C
C CONSTRUCT HERE DENSITY MATRICES IN AO BASIS NEEDED FOR THE TWO-ELECTRON 
C CONTRIBUTION TO THE GRADIENT
C
       CALL DSYM2(WORK1(KHFSO),WORK1(KRELSO),WORK1(KNHFSO),NBASIS,
     &            IPRALL,WORK1(KHFA),WORK1(KHFB),WORK1(KRELA),
     &            WORK1(KRELB),WORK1(KNHFA),WORK1(KNHFB),MAXLOP)
C
C CONTRUCT HERE DENSITY MATRICES IN AO BASIS NEEDED FOR FOCK-MATRIX
C CONSTRUCTION
C
       IF(FOCK) THEN
C
C FIRST CALCULATE FOLDED DENSITY MATRIX
C
        CALL SSCAL((1+IUHF)*NBASTT,TWO,WORK1(KHFSO),1)
C
C THEN TRANSFORM
C
        CALL DSYM2A(WORK1(KHFSO),NBASIS,IPRALL,WORK1(KHFFA),
     &              WORK1(KHFFB),MAXLOP)
C
C FOR CORRELATED SECOND DERIVATIVES WE NEED ALSO THE FOLDED
C RELAXED DENSITY
C
        IF(.NOT.SCF) THEN
         CALL SSCAL((1+IUHF)*NBASTT,TWO,WORK1(KRELSO),1)
         CALL DSYM2A(WORK1(KRELSO),NBASIS,IPRALL,WORK1(KRELFA),
     &               WORK1(KRELFB),MAXLOP)
        ENDIF
       ENDIF
      ENDIF
C
C  SORT TWO-PARTICLE DENSITY MATRIX
C
      IF (RUNSOR) THEN
       CALL PSORG(WORK1(KEND),LWORK1+1-KEND)
      END IF
C
C  CALCULATE EXPECTAION VALUE, FOCK MATRICES, ETC.
C
      IF (RUNINT) THEN
       IF(.NOT.MAGNET) THEN
         CALL TWOINT(2,MAXDIF,FOCK,0,NODV,NOPV,TKTIME,
     &               IPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,RETUR,
     &               WORK1(KEND),LWORK1+1-KEND,ENEREE,SCF_GRADEE,
     &               GRADEE,SCF_HESSEE,HESSEE,NCOORD,WORK1(KHFA),
     &               WORK1(KHFB),
     &               WORK1(KRELA),WORK1(KRELB),WORK1(KNHFA),
     &               WORK1(KNHFB),WORK1(KHFFA),WORK1(KHFFB),
     &               WORK1(KRELFA),WORK1(KRELFB),NBASIS,MAXLOP)
       ELSE
        CALL NMR2DR(FOCK,0,NODV,NOPV,
     &              IPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,RETUR,
     &              WORK1(KEND),LWORK1+1-KEND,MSZ,3,
     &              WORK1(KHFA),WORK1(KHFB),WORK1(KRELA),
     &              WORK1(KRELB),WORK1(KNHFA),WORK1(KNHFB),
     &              WORK1(KHFFA),WORK1(KHFFB),
     &              WORK1(KRELFA),WORK1(KRELFB),NBASIS,MAXLOP)
       ENDIF
CSSS       IF (RUNSOR) CLOSE(LUPAO,STATUS='DELETE')
      END IF
      RETURN
      END
