











      SUBROUTINE EVAL_QSTLST_PATH_RIC(QST_TANGENT, LST_TANGENT, WORK,
     &                                PMAT, LST, QST)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)

C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)


c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end



C coord.com : begin
C
      DOUBLE PRECISION Q, R, ATMASS
      INTEGER NCON, NR, ISQUASH, IATNUM, IUNIQUE, NEQ, IEQUIV,
     &        NOPTI, NATOMS
      COMMON /COORD/ Q(3*MXATMS), R(MAXREDUNCO), NCON(MAXREDUNCO),
     &     NR(MXATMS),ISQUASH(MAXREDUNCO),IATNUM(MXATMS),
     &     ATMASS(MXATMS),IUNIQUE(MAXREDUNCO),NEQ(MAXREDUNCO),
     &     IEQUIV(MAXREDUNCO,MAXREDUNCO),
     &     NOPTI(MAXREDUNCO), NATOMS

C coord.com : end



      LOGICAL XYZIN, NWFINDIF, QST, LST
      DOUBLE PRECISION LST_NORM_SQR, LST_NORM_INV, LST_TANGENT
      INTEGER TOTREDNCO

      COMMON /USINT/NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      COMMON /INPTYP/XYZIN,NWFINDIF

      DIMENSION WORK(NX*NX), LST_TANGENT(*), QST_TANGENT(*),
     &          PMAT(NX*NX)
C
C Use Eqn. 1 and 8 in Schlegel to evaluate the QST and LST normalized
C tangent vectors.
C
C Read the guess structures (user provided) for the transition
C state and the products.
C 
      CALL GETREC(20, 'JOBARC', 'REDNCORD', 1, TOTREDNCO)

      IRECTNT = 1
      IPRDUCT = IRECTNT + TOTREDNCO
      ICURENT = IPRDUCT + TOTREDNCO
      ISCR    = ICURENT + TOTREDNCO
      IEND    = ISCR    + TOTREDNCO

      IF (IEND .GT. NX*NX) CALL INSMEM("@EVAL_QSTLST_PATH_RIC",
     &                                  IEND, NX*NX)

C Cautionary note: The reactant coordinates are copied to the
C memory location desiganted as IPRDUCT. This was as a result
C of a earlier confusion in designating multiple structures. This
C has no effect on actual work done, Ajith Perera, 12/2012

      CALL GETREC(20,'JOBARC','RX_RICS_',TOTREDNCO*IINTFP,
     &            WORK(IPRDUCT))

CSSS      CALL CNVRT2RAD(WORK(IPRDUCT),3*NATOMS)
CSSS      CALL SQUSH(WORK(IPRDUCT),3*NATOMS)
C
      Write(6,"(a)") "Debug-Info: The Reactant Structure"
      Write(6,"(3(1x,F12.6))") (Work(IPRDUCT + I),I= 0,TOTREDNCO-1)

      IF (LST) THEN
        CALL GETREC(1,'JOBARC','PR_RICS_',TOTREDNCO*IINTFP,
     &              WORK(IRECTNT))

      Write(6,"(a)") "Debug-Info: Product Structure"
      Write(6,"(3(1x,F12.6))") (Work(IRECTNT + I), I= 0,TOTREDNCO-1)
        CALL DAXPY(TOTREDNCO, -1.0D0, WORK(IPRDUCT), 1,
     &             WORK(IRECTNT), 1)

      Write(6,"(a)") "Debug-Info: The Diff. Structure"
      Write(6,"(3(1x,F12.6))") (Work(IRECTNT + I), I= 0,TOTREDNCO-1)

CSSS            CALL DCOPY(NX, WORK(IPRDUCT), 1, LST_TANGENT, 1)
CSSS            CALL GETREC(20, "JOBARC", "BMATRIX ", 3*NATOMS*NOPT*
CSSS     &                  IINTFP, WORK(1))
CSSS            IOFFSET = 3*NATOMS*NOPT + 1
CSSS            CALL XGEMM("N", "N", NOPT, 1, 3*NATOMS, 1.0D0, WORK,
CSSS     &                  NOPT, LST_TANGENT, 3*NATOMS, 0.0D0, WORK
CSSS     &                  (IOFFSET), NOPT)
CSSS            CALL DCOPY(NOPT, WORK(IOFFSET), 1, LST_TANGENT, 1)
CSSS        CALL SYMUNQONLY(WORK(IPRDUCT), LST_TANGENT)

        CALL XGEMM("N", "N", TOTREDNCO, 1, TOTREDNCO, 1.0D0, 
     &              PMAT, TOTREDNCO, WORK(IRECTNT),
     &              TOTREDNCO, 0.0D0, WORK(ISCR), TOTREDNCO)
        CALL DCOPY(TOTREDNCO, WORK(ISCR), 1, WORK(IRECTNT), 1)

        CALL NORMAL(WORK(IRECTNT), TOTREDNCO)
        TNORM = DDOT(TOTREDNCO, WORK(IRECTNT),1,WORK(IRECTNT), 1)
C
        IF (TNORM .NE. 1.0D0) THEN
           WRITE(6,*)
           WRITE(6,"(a,a)") "Norm condition for the LST tangent is ",
     &                      "incorrect! Can not proceed."
            CALL ERREX
        ENDIF
C
        CALL DCOPY(TOTREDNCO, WORK(IRECTNT), 1, LST_TANGENT, 1)

            Write(6,"(a)") "Debug-Info: LST tangent"
            Write(6,"(3(1x,F12.6))") (LST_TANGENT(I), I= 1, NOPT)

      ELSE IF (QST) THEN

C Cautionary note: The reactant coordinates are copied to the
C memory location desiganted as IRECTNT. This was as a result
C of a earlier confusion in designating multiple structures. This
C has no effect on actual work done, Ajith Perera, 12/2012

        CALL GETREC(1,'JOBARC','PR_RICS_',TOTREDNCO*IINTFP,
     &              WORK(IRECTNT))

CSSS        CALL CNVRT2RAD(WORK(IRECTNT), 3*NATOMS)
CSSS        CALL SQUSH(WORK(IRECTNT), 3*NATOMS)
C
        CALL DCOPY(TOTREDNCO, R, 1, WORK(ICURENT), 1)

      Write(6,"(a)") "Debug-Info: Product Structure"
      Write(6,"(3(1x,F12.6))") (Work(IRECTNT + I), I= 0,TOTREDNCO-1)
      Write(6,"(a)") "Debug-Info: Reactant Structure"
      Write(6,"(3(1x,F12.6))") (Work(IPRDUCT + I), I= 0,TOTREDNCO-1)
      Write(6,"(a)") "Debug-Info: Current Structure"
      Write(6,"(3(1x,F12.6))") (Work(ICURENT + I), I= 0,TOTREDNCO-1)

        CALL DAXPY(TOTREDNCO, -1.0D0, WORK(ICURENT), 1,
     &             WORK(IPRDUCT), 1)
        CALL DAXPY(TOTREDNCO, -1.0D0, WORK(ICURENT), 1,
     &             WORK(IRECTNT), 1)

        QST_NORM_SQR_PX = DDOT(TOTREDNCO, WORK(IPRDUCT), 1,
     &                         WORK(IPRDUCT), 1)
        QST_NORM_SQR_RX = DDOT(TOTREDNCO, WORK(IRECTNT), 1,
     &                         WORK(IRECTNT), 1)
      Write(6,"(a)") "Debug-Info: P-X tructure"
      Write(6,"(3(1x,F12.6))") (Work(IPRDUCT + I), I= 0,TOTREDNCO-1)
      Write(6,"(a)") "Debug-Info: R-X Structure"
      Write(6,"(3(1x,F12.6))") (Work(IRECTNT + I), I= 0,TOTREDNCO-1)
        QST_NORM_SQR_PX = DDOT(TOTREDNCO, WORK(IPRDUCT), 1,
     &                         WORK(IPRDUCT), 1)
        QST_NORM_SQR_RX = DDOT(TOTREDNCO, WORK(IRECTNT), 1,
     &                         WORK(IRECTNT), 1)
        QST_NORM_PX_RX  = DDOT(TOTREDNCO, WORK(IRECTNT), 1,
     &                         WORK(IPRDUCT), 1)

        QST_NORM_PX = DSQRT(QST_NORM_SQR_PX)
        QST_NORM_RX = DSQRT(QST_NORM_SQR_RX)

        IF ((QST_NORM_SQR_PX.LE.1.0D-14).OR.
     &     (QST_NORM_SQR_RX.LE.1.0D-14)    ) THEN
          Write(6,*)
          Write(6,"(a)")"AN INPUT ERROR: INSPECT THE INPUT STRUCTURES"
          CALL ERREX
        END IF

        QST_NORM_PX_INV = 1.0D0/(QST_NORM_SQR_PX)
        QST_NORM_RX_INV = 1.0D0/(QST_NORM_SQR_RX)

        CALL DSCAL(TOTREDNCO, QST_NORM_PX_INV, WORK(IPRDUCT), 1)
        CALL DSCAL(TOTREDNCO, QST_NORM_RX_INV, WORK(IRECTNT), 1)

        CALL DAXPY(TOTREDNCO, -1.0D0, WORK(IRECTNT), 1,
     &             WORK(IPRDUCT), 1)
        DENOM = (QST_NORM_SQR_PX + QST_NORM_SQR_RX - 
     &           2.0D0*QST_NORM_PX_RX)

        IF (DENOM .LE. 1.0D-14) THEN
           Write(6,*)
           Write(6,"(a)")"AN INPUT ERROR: INSPECT THE INPUT STRUCTURES"
           CALL ERREX
        ENDIF

        A_SCALE = DSQRT((QST_NORM_SQR_PX*QST_NORM_SQR_RX)/DENOM)

      Write(6,*)
      Write(6,"(a)") "Various Norms"
      Write(6,"(a,2(1x,F12.6))") "QST_NORM_SQR_PX,RX :",QST_NORM_SQR_PX,
     &                          QST_NORM_SQR_RX
      Write(6,"(a,2(1x,F12.6))") "QST_NORM_PX,RX     :",QST_NORM_PX,
     &                          QST_NORM_RX
      Write(6,"(a,(1x,F12.6))") "QST_NORM_PX_RX     :",QST_NORM_PX_RX
CSSS            CALL DCOPY(3*NATOMS, WORK(IPRDUCT), 1, QST_TANGENT, 1)
CSSS            CALL GETREC(20, "JOBARC", "BMATRIX ", 3*NOPT*NATOMS*
CSSS     &                  IINTFP, WORK(1))
CSSS            IOFFSET = 3*NATOMS*NOPT + 1
CSSS            CALL XGEMM("N", "N", NOPT, 1, 3*NATOMS, 1.0D0, WORK,
CSSS     &                  NOPT, QST_TANGENT, 3*NATOMS, 0.0D0, WORK
CSSS     &                  (IOFFSET), NOPT)

        CALL DSCAL(TOTREDNCO, A_SCALE, WORK(IPRDUCT), 1)
C
        CALL XGEMM("N", "N", TOTREDNCO, 1, TOTREDNCO, 1.0D0,
     &              PMAT, TOTREDNCO, WORK(IPRDUCT),
     &              TOTREDNCO, 0.0D0, WORK(ISCR), TOTREDNCO)
        CALL DCOPY(TOTREDNCO, WORK(ISCR), 1, WORK(IPRDUCT), 1)

        CALL NORMAL(WORK(IPRDUCT), TOTREDNCO)

        TNORM = DDOT(TOTREDNCO, WORK(IPRDUCT),1,WORK(IPRDUCT), 1)
C
        IF (TNORM .NE. 1.0D0) THEN
           WRITE(6,*)
           WRITE(6,"(a,a)") "Norm condition for the QST tangent is ",
     &                      "incorrect! Can not proceed."
           CALL ERREX
        ENDIF
C
      Write(6,"(a)") "Checking the Norm of QST tangent"
      Write(6,"(a,1x,F12.6)")"Tangent: ", TNORM
C
CSSS        CALL SYMUNQONLY(WORK(IPRDUCT), QST_TANGENT)
          
         CALL DCOPY(TOTREDNCO, WORK(IPRDUCT), 1, QST_TANGENT, 1)
C
      Write(6,*)
      Write(6,"(a)") "Debug-Info: QST tangent"
      Write(6,"(3(1x,F12.6))") (QST_TANGENT(I), I= 1, NOPT)
         END IF

       RETURN
       END

