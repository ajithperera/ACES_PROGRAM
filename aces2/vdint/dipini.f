










      SUBROUTINE DIPINI
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

      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2,
     *           MXAQNS=MXAQN*MXAQN*MXCONT*MXCONT)
      COMMON /SDER/ SDER0 (MXAQNS),
     *              SDERX (MXAQNS), SDERY (MXAQNS), SDERZ (MXAQNS),
     *              SDERXX(MXAQNS), SDERXY(MXAQNS), SDERXZ(MXAQNS),
     *              SDERYY(MXAQNS), SDERYZ(MXAQNS), SDERZZ(MXAQNS)
      COMMON /TDER/ TDER0 (MXAQNS),
     *              TDERX (MXAQNS), TDERY (MXAQNS), TDERZ (MXAQNS),
     *              TDERXX(MXAQNS), TDERXY(MXAQNS), TDERXZ(MXAQNS),
     *              TDERYY(MXAQNS), TDERYZ(MXAQNS), TDERZZ(MXAQNS)
      COMMON /ADER/ ADER0 (MXAQNS),
     *              IA0000, IA0X00, IA0Y00, IA0Z00,
     *              IAXX00, IAXY00, IAXZ00, IAYY00,
     *              IAYZ00, IAZZ00, IA000X, IA000Y,
     *              IA000Z, IA00XX, IA00XY, IA00XZ,
     *              IA00YY, IA00YZ, IA00ZZ, IA0X0X,
     *              IA0X0Y, IA0X0Z, IA0Y0X, IA0Y0Y,
     *              IA0Y0Z, IA0Z0X, IA0Z0Y, IA0Z0Z
      COMMON /DDER/ XINT0(MXAQNS),
     *              XINTX(MXAQNS), XINTY(MXAQNS), XINTZ(MXAQNS),
     *              YINT0(MXAQNS),
     *              YINTX(MXAQNS), YINTY(MXAQNS), YINTZ(MXAQNS),
     *              ZINT0(MXAQNS),
     *              ZINTX(MXAQNS), ZINTY(MXAQNS), ZINTZ(MXAQNS),
     *              SINT0(MXAQNS)
      COMMON /RDER/ IR00, IR0X, IR0Y, IR0Z,
     *              IRXX, IRXY, IRXZ, IRYY, IRYZ, IRZZ
      CALL ZERO(XINT0, (13*MXAQNS) )
      RETURN
      END
