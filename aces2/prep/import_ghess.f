










      PROGRAM IMPORT_GHESS
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C


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



C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
C
      DIMENSION FCM(9*MXATMS*MXATMS)
C
      CALL ACES_INIT_RTE
      CALL ACES_JA_INIT

      CALL GETREC(1,'JOBARC','NREALATM', 1,  NREALS) 
C
      IR     = 33
      NC1    = 3*NREALS
      NC2    = NC1
      ISTART = 1
    
      OPEN(UNIT=IR, FILE="G_HESS",FORM="FORMATTED")
C
      ICC = 0
      DO  I=1,NC1
          IC=0
          DO MINCOL = 1, NC1, 5
             MAXCOL=MINCOL+4
             IF(MAXCOL.GT.NC1) MAXCOL=NC1
             IC=IC+1
C
             READ(IR,9020) II,ICC,
     *           (FCM(ISTART+(J-1)*NC2),J=MINCOL,MAXCOL)
             IMD100 = MOD(I,100)
C
               IF (II .NE. IMD100 .AND. ICC .NE. IC) THEN
                  WRITE(IW,9030) I,IC,II,ICC
                  CALL ERREX
               END IF
          ENDDO
            ISTART=ISTART+1
      ENDDO
C
      Write(6, "(a)") "Importing a Hessian from external sources"
      CALL PUTREC(20, 'JOBARC', 'HESSIANM', 9*NREALS*NREALS*IINTFP, 
     &            FCM)
C
      NX = 3*NREALS
      CALL OUTPUT(FCM, 1, NX, 1, NX, NX, NX, 1)
C             
 9020 FORMAT(I2,I3,5F15.8)
 9030 FORMAT(' *** ERROR READING FORCE CONSTANT MATRIX ELEMENT',2I4/
     *       '           INPUT WAS',2I4/)
      CALL ACES_JA_FIN
C
      STOP 
      END
