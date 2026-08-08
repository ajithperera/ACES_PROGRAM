










      SUBROUTINE CLCNRXLFT(NROWXLFT)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C     Argument.
C-----------------------------------------------------------------------
      INTEGER NROWXLFT
C-----------------------------------------------------------------------
C     Common block variables.
C-----------------------------------------------------------------------
      DOUBLE PRECISION PC
      INTEGER DSTRT,NTAP,LU2,NRSS,NUCZ,ITAG,MAXLOP,MAXLOT,KMAX,NMAX,
     &        KHKT,MULT,ISYTYP,ITYPE,AND,OR,EOR,NPARSU,NPAR,MULNUC,
     &        NHKT,MUL,NUCO,NRCO,JSTRT,NSTRT,MST,JRS
C-----------------------------------------------------------------------
C     Local variables.
C-----------------------------------------------------------------------
      INTEGER IA,IB,ITOT,NC,ND
C-----------------------------------------------------------------------
C     Parameters.
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

C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
C-----------------------------------------------------------------------
      COMMON /INDX/ PC(512),DSTRT(8,MXCBF),NTAP,LU2,NRSS,NUCZ,ITAG,
     & MAXLOP,MAXLOT,KMAX,NMAX,KHKT(7),MULT(8),ISYTYP(3),ITYPE(7,28),
     & AND(8,8),OR(8,8),EOR(8,8),NPARSU(8),NPAR(8),MULNUC(Mxatms),
     &  NHKT(MXTNSH),  MUL(MXTNSH),NUCO(MXTNSH),NRCO(MXTNSH),
     & JSTRT(MXTNSH),NSTRT(MXTNSH), MST(MXTNSH), JRS(MXTNSH)
C
C-----------------------------------------------------------------------
CJDW 1/28/97.
C     Block of code to determine maximum leading dimension of XLFT.
C     Previously this was hardcoded to 4095, which is not enough for
C     the extremely common case of 14 d primitives, for example.
C
      NROWXLFT = 0
      DO 100 IA=1,KMAX
      DO  90 IB=1,IA
C
      ITOT=0
      DO  80 NC=1,NHKT(IA)
      DO  70 ND=1,NHKT(IB)
C
      ITOT = ITOT + NC + ND - 1
   70 CONTINUE
   80 CONTINUE
      ITOT = ITOT * NUCO(IA) * NUCO(IB)
C
      NROWXLFT = MAX(NROWXLFT,ITOT)
   90 CONTINUE
  100 CONTINUE
C
c     write(6,*) kmax,nhkt,nuco,itot
c     write(6,*) '  @CLCNRXLFT-I, NROWXLFT ',NROWXLFT
      IF(NROWXLFT .LT. 4095) NROWXLFT = 4095
cYAU      write(6,*) '  @CLCNRXLFT-I, NROWXLFT ',NROWXLFT
C-----------------------------------------------------------------------
      RETURN
      END
