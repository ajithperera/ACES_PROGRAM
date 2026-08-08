










      SUBROUTINE CMPIN(A,AC,TOL,N,IH,T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

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

      parameter (nh4=4*nht-3, nh2=nht+nht+1)
      parameter (mxp2=maxprim*maxprim)
      parameter (mxp4=mxp2*mxp2)
      DIMENSION A(2),B(2),AC(2),BC(2),T(2)
      COMMON /ICOMP/ IBE(mxp4),K,NN
      common /icomp2/IBEM
      NN=N
      DO 1 I=1,N
    1   T(I) = ABS(A(I))
      CALL WHENFGT(N,T,1,TOL,IBE,IBEM)
      CALL GATHER(IBEM,AC,A,IBE)
      IH=IBEM
      RETURN
      ENTRY CMPR1(A,AC)
      CALL GATHER(IBEM,AC,A,IBE)
      RETURN
      ENTRY CMPR2(A,B,AC,BC)

      CALL GATHER(IBEM,AC,A,IBE)
      CALL GATHER(IBEM,BC,B,IBE)
      RETURN
      ENTRY EXPND(A,AC)
      CALL SCATTER(IBEM,A,IBE,AC)
      RETURN
      ENTRY CMPRI(A,AC,M)
      DO 44 MM=0,M-1
   44   CALL GATHER(IBEM,AC(MM*ibem+1),A(MM*NN+1),IBE)
      RETURN
      END
