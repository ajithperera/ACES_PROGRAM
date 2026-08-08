










      SUBROUTINE SAVVR (AREA,S2,CSS2,NAQRS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C     Parameters
C-----------------------------------------------------------------------

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

      parameter (nh4=4*NHT-3, nh2=NHT+NHT+1)
      parameter (mxp2=maxprim*maxprim)
C-----------------------------------------------------------------------
      COMMON /ONE/ 
     1 ALP(mxp2),BET(mxp2),S1(mxp2),CSS1(mxp2),PV(mxp2,8,3),CAB(mxp2),
     2 GAM(mxp2),DEL(mxp2),SS2(mxp2),Cjnk(mxp2),QV(mxp2,8,3),CCD(mxp2)
c      COMMON /ONE/ 
c     1   ALP(mxp2),BET(mxp2),S1(mxp2),CSS1(mxp2),PV(mxp2,8,3),CAB(mxp2)
      DIMENSION AREA(NAQRS,5),S2(2),CSS2(2)
      COMMON /CON/ RAX,RAY,RAZ,RPRX,RPRY,RPRZ,DISTPQ,DISTRS,TOLR,TP52,
     &             ZJUNK(MXP2,16)
      COMMON /CONI/ NUCT(4),NRCT(4),JSTT(4),
     1 IFD1X,IFD2X,NUCAQ,NUCRS,NXXXX(6)
      CALL MVSUM  (AREA(1,3),S1,NUCAQ,S2,NUCRS)
      CALL MVPROD (AREA(1,4),S1,NUCAQ,S2,NUCRS)
CDIR$ IVDEP
*VOCL LOOP,NOVREC
      DO 2 I = 1,NAQRS
      AREA(I,3) = 1./AREA(I,3)
    2 AREA(I,2) = AREA(I,4) * AREA(I,3)
      CALL MVPROD (AREA(1,4),CSS1,NUCAQ,CSS2,NUCRS)
CDIR$ IVDEP
*VOCL LOOP,NOVREC
      DO 3 I=1,NAQRS
      AREA(I,1)= AREA(I,4)*SQRT(AREA(I,3))*TP52
    3 CONTINUE
      RETURN
      END
