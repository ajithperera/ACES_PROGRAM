










      SUBROUTINE TRAGRD(SGRAD,CGRAD,NCOOR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c IMPORTANT: The ECP library contains three routines (doitgr, ecpabl, ecpder)
c that use this parameter but do not include this file. Changing this value
c requires updating those files as well. [Yau: If we linked them together,
c then compiling ecp would require the vdint directory. Alternatively, we
c could put the mxcent.par file into the top-level include directory.]

C MXCENT : Maximum number of atoms currently allowed
C MXCOOR : The number of Cartesian that correspond to MXCENT

      INTEGER MXCENT, MXCOOR
      PARAMETER (MXCENT=200, MXCOOR=3*MXCENT)
      DIMENSION SGRAD(NCOOR), CGRAD(NCOOR)
      COMMON /TRANUC/ TRCTOS(MXCOOR,MXCOOR), TRSTOC(MXCOOR,MXCOOR)
      DATA AZERO,ONE /0.D0,1.D0/
C
      CALL XGEMM('N','N',NCOOR,1,NCOOR,ONE,
     *           TRSTOC,MXCOOR,SGRAD,NCOOR,AZERO,
     *           CGRAD,NCOOR)
      RETURN
      END
