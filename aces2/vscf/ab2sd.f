










      SUBROUTINE AB2SD(DENSA,DENSB,DENSSD,LEN)
C
C  This routine constructs the singles and doubles density matrices from the
C  alpha and beta density matrices.  This is done using the formula:
C
C
C      D(doubles) = 2*D(beta)
C
C      D(singles) = D(alpha) - D(beta)
C
CEND
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION DENSA(LEN),DENSB(LEN),DENSSD(2*LEN)
      COMMON /POPUL/ NOCCA(8),NOCCB(8)
c symm2.com : begin

c This is initialized in vscf/symsiz.

c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end
      integer nirrep,      nbfirr(8),   irpsz1(36),  irpsz2(28),
     &        irpds1(36),  irpds2(56),  irpoff(9),   ireps(9),
     &        dirprd(8,8), iwoff1(37),  iwoff2(29),
     &        inewvc(maxbasfn),         idxvec(maxbasfn),
     &        itriln(9),   itriof(8),   isqrln(9),   isqrof(8),
     &        mxirr2
      common /SYMM2/ nirrep, nbfirr, irpsz1, irpsz2, irpds1, irpds2,
     &               irpoff, ireps,  dirprd, iwoff1, iwoff2, inewvc,
     &               idxvec, itriln, itriof, isqrln, isqrof, mxirr2
c symm2.com : end
C
      DATA ONEM /-1.0/
      DATA TWO /2.0/
C
      DO 100 I=1,NIRREP
        IF(NBFIRR(I).EQ.0) GOTO 100
        IF(NOCCA(I).EQ.NOCCB(I)) THEN
          CALL SCOPY(ITRILN(I),DENSB(ITRIOF(I)),1,DENSSD(ITRIOF(I)),1)
          CALL SSCAL(ITRILN(I),TWO,DENSSD(ITRIOF(I)),1)
          CALL ZERO(DENSSD(LEN+ITRIOF(I)),ITRILN(I))
        ELSE
          CALL SCOPY(ITRILN(I),DENSA(ITRIOF(I)),1,
     &               DENSSD(LEN+ITRIOF(I)),1)
          CALL SAXPY(ITRILN(I),ONEM,DENSB(ITRIOF(I)),1,
     &               DENSSD(LEN+ITRIOF(I)),1)
        ENDIF
  100 CONTINUE
C
      RETURN
      END
