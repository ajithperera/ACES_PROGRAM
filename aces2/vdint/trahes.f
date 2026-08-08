










      SUBROUTINE TRAHES(SHESS,CHESS,NCOOR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION MATRIX

c IMPORTANT: The ECP library contains three routines (doitgr, ecpabl, ecpder)
c that use this parameter but do not include this file. Changing this value
c requires updating those files as well. [Yau: If we linked them together,
c then compiling ecp would require the vdint directory. Alternatively, we
c could put the mxcent.par file into the top-level include directory.]

C MXCENT : Maximum number of atoms currently allowed
C MXCOOR : The number of Cartesian that correspond to MXCENT

      INTEGER MXCENT, MXCOOR
      PARAMETER (MXCENT=200, MXCOOR=3*MXCENT)
      PARAMETER (HALF = 0.5D0)
      DIMENSION SHESS(NCOOR,NCOOR), CHESS(MXCOOR,MXCOOR),
     *          MATRIX(MXCOOR,MXCOOR)
      COMMON /TRANUC/ TRCTOS(MXCOOR,MXCOOR), TRSTOC(MXCOOR,MXCOOR)
C
      DATA AZERO,ONE /0.D0,1.D0/
C
      DO 100 I = 1, NCOOR
         DO 200 J = 1, I - 1
            SHESS(J,I) = SHESS(I,J)
  200    CONTINUE
  100 CONTINUE
      CALL XGEMM('N','N',NCOOR,NCOOR,NCOOR,ONE,
     *           TRSTOC,MXCOOR,SHESS,NCOOR,AZERO,
     *           MATRIX,MXCOOR)
      CALL XGEMM('N','T',NCOOR,NCOOR,NCOOR,ONE,
     *           MATRIX,MXCOOR,TRSTOC,MXCOOR,
     *           AZERO,CHESS,MXCOOR)
      RETURN
      END
