










      SUBROUTINE MKFULL(WFULL,WCHUNK,SCR,LDIM1,LDIM2,NBAS,JFLAG)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION WFULL(NBAS,NBAS),WCHUNK(LDIM1)
      DIMENSION SCR(LDIM2)
C
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
      INDX2(I,J,N)=I+(J-1)*N
C
      IF(JFLAG.EQ.0) THEN
        DO 100 I=1,NIRREP
          DO 110 J=1,NBFIRR(I)
            DO 120 K=1,NBFIRR(I)
              WFULL(IRPOFF(I)+J,IRPOFF(I)+K)=WCHUNK(ISQROF(I)-1+
     &              INDX2(J,K,NBFIRR(I)))
  120       CONTINUE
  110     CONTINUE
  100   CONTINUE
      ELSEIF(JFLAG.EQ.1) THEN
        DO 101 I=1,NIRREP
          CALL EXPND2(WCHUNK(ITRIOF(I)),SCR,NBFIRR(I))
          DO 111 J=1,NBFIRR(I)
            DO 121 K=1,NBFIRR(I)
              WFULL(IRPOFF(I)+J,IRPOFF(I)+K)=SCR(INDX2(J,K,NBFIRR(I)))
  121       CONTINUE
  111     CONTINUE
  101   CONTINUE
      ENDIF
C
      RETURN
      END
