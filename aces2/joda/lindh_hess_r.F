      SUBROUTINE LINDH_HESS_R(HC, NX)
C
C Lindh model Hessian (LINDH_HESS_CART) at the Cartesian optimizer
C coordinates R of /COORD/ (GEOM_OPT=CART). 10/2026.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
#include "mxatms.par"
#include "coord.com"
C
      DIMENSION HC(NX,NX)
C
      CALL LINDH_HESS_CART(R, HC, NATOMS, IATNUM)
C
#ifdef _DEBUG_LVL0
      Write(6,"(a)") " @LINDH_HESS_R: model Cartesian Hessian"
      CALL OUTPUT(HC, 1, NX, 1, NX, NX, NX, 1)
#endif
      RETURN
      END
