










      SUBROUTINE LINDH_HESS_R(HC, NX)
C
C Lindh model Hessian (LINDH_HESS_CART) at the Cartesian optimizer
C coordinates R of /COORD/ (GEOM_OPT=CART). 10/2026.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
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


C
      DIMENSION HC(NX,NX)
C
      CALL LINDH_HESS_CART(R, HC, NATOMS, IATNUM)
C
      Write(6,"(a)") " @LINDH_HESS_R: model Cartesian Hessian"
      CALL OUTPUT(HC, 1, NX, 1, NX, NX, NX, 1)
      RETURN
      END
