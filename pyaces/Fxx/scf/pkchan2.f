










      SUBROUTINE PKCHAN2(ICHAIN,PK,IDISSIZ,IDISNUM,BUF,IBUF,SCR,IBKSTRT,
     &                   IBKOF,MAXBUK,IBKDIS,IBKSIZ,NBUCK,NBKINT,MAXDIS,
     &                   IUHF)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER(LUSRT2=22)
C
      DIMENSION ICHAIN(3*NBUCK),IDISNUM(IBKDIS,2)
      DIMENSION PK(IBKSIZ),IDISSIZ(IBKDIS)
      DIMENSION IBKSTRT(MAXBUK),IBKOF(MAXBUK)
      DIMENSION BUF(NBKINT),IBUF(NBKINT)
      DIMENSION SCR(MAXDIS)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
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
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
C
      DATA ONE /1.0/
C
      INDX(I,J)=J+(I*(I-1))/2
C
      CALL TIMER(1)
C
      DO 50 IBK=1,NBUCK
        CALL ZERO(PK,IBKSIZ)
  198   IF(ICHAIN(NBUCK+IBK).EQ.0) GOTO 310
  199   READ(LUSRT2,REC=ICHAIN(NBUCK+IBK))BUF,IBUF,NUT,ICHAN
        ICHAIN(NBUCK+IBK)=ICHAN
        DO 200 INT=1,NUT
          IADR=IBUF(INT)-IBKOF(IBK)
          PK(IADR)=PK(IADR)+BUF(INT)
  200   CONTINUE
        IF(ICHAIN(NBUCK+IBK).NE.0) GOTO 199
C
  310   CONTINUE
C
        ISTART=IBKSTRT(IBK)
        IEND=IBKSTRT(IBK+1)-1
        IOFF=1
        DO 400 I=ISTART,IEND
          CALL GETLST(SCR,IDISNUM(I,1),1,2,1,IDISNUM(I,2))
          CALL SAXPY(IDISSIZ(I),ONE,PK(IOFF),1,SCR,1)
          CALL PUTLST(SCR,IDISNUM(I,1),1,2,1,IDISNUM(I,2))
          IOFF=IOFF+IDISSIZ(I)
  400   CONTINUE
C
   50 CONTINUE
C
      CLOSE(UNIT=LUSRT2,STATUS='DELETE')
C
      CALL TIMER(1)
      WRITE(LUOUT,5000)TIMENEW
 5000 FORMAT(T3,'@PKCHAN2-I, Chaining of PK lists required ',F10.3,
     &          ' seconds.',/)
C
      RETURN
      END
