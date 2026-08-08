










      SUBROUTINE REOREV(EVECIN,EVECOT,NBAS,ISPIN)
C
C  This routine reorders the eigenvectors from correlated order
C  to that required for the SCF.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end
C
      CHARACTER*1 ISP(2)
C
      DIMENSION EVECIN(NBAS,NBAS),EVECOT(NBAS,NBAS),ireordr(MAXBASFN)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      DATA ISP /'A','B'/
C
C ORDER EIGENVECTORS BY IRREP
C
      call getrec(-1,'JOBARC','REORDER'//isp(ispin),nbas,ireordr)


      ITOT=0
      DO 100 I=1,NBAS
        ITOT=ITOT+IREORDR(I)
  100 CONTINUE
      IF(ITOT.EQ.0) THEN
        CALL SCOPY(NBAS*NBAS,EVECIN,1,EVECOT,1)
        RETURN
      ENDIF
C
      do 10 i=1,nbas
       jtar=ireordr(i)
       call scopy(nbas,evecin(1,jtar),1,evecot(1,i),1)
10    continue
C
      call scopy(nbas*nbas,evecot,1,evecin,1)
      return
      end
