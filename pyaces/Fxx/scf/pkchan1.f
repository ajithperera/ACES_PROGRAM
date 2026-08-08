










      SUBROUTINE PKCHAN1(ICHAIN,PK,IDISSIZ,IDISNUM,BUF,IBUF,IBKSTRT,
     &                   IBKOF,MAXBUK,IBKDIS,IBKSIZ,NBUCK,NBKINT,IUHF)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER(LUSRT1=21)
C
      DIMENSION ICHAIN(3*NBUCK),IDISNUM(IBKDIS,2)
      DIMENSION PK(IBKSIZ),IDISSIZ(IBKDIS)
      DIMENSION IBKSTRT(MAXBUK),IBKOF(MAXBUK)
      DIMENSION BUF(NBKINT),IBUF(NBKINT)
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
      INDX(I,J)=J+(I*(I-1))/2
C
      CALL TIMER(1)
C
C  Initialize all of the PK lists before we get the actual numbers.
C
c      CALL UPDMOI(0,0,0,0,1,0)
      call aces_io_remove(50,'MOINTS')
      DO 10 I=1,NIRREP
        DO 11 J=1,I
          IF(I.EQ.J) THEN
            IJUNK=INDX(I,I)
            CALL UPDMOI(IRPDS1(IJUNK),IRPDS1(IJUNK),1,IJUNK,0,0)
            IF(IUHF.NE.0) THEN
              CALL UPDMOI(IRPDS1(IJUNK),IRPDS1(IJUNK),2,IJUNK,0,0)
            ENDIF
          ELSE
            IJUNK=INDX(I-1,J)
            CALL UPDMOI(IRPDS2(2*IJUNK),IRPDS2(2*IJUNK-1),1,
     &                  36+IJUNK,0,0)
            IF(IUHF.NE.0) THEN
              CALL UPDMOI(IRPDS2(2*IJUNK),IRPDS2(2*IJUNK-1),2,
     &                    36+IJUNK,0,0)
            ENDIF
          ENDIF
   11   CONTINUE
   10 CONTINUE
C
      DO 50 IBK=1,NBUCK
        CALL ZERO(PK,IBKSIZ)
        IF(ICHAIN(IBK).EQ.0) GOTO 310
   99   READ(LUSRT1,REC=ICHAIN(IBK))BUF,IBUF,NUT,ICHAN
        ICHAIN(IBK)=ICHAN
        DO 100 INT=1,NUT
          IADR=IBUF(INT)-IBKOF(IBK)
          PK(IADR)=PK(IADR)+BUF(INT)
  100   CONTINUE
        IF(ICHAIN(IBK).NE.0) GOTO 99
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  If this is not an RHF calculation, then dump out the J integrals as
C  a separate list in addition to their contribution to the P list.
C
        IF(IUHF.NE.0) THEN
          ISTART=IBKSTRT(IBK)
          IEND=IBKSTRT(IBK+1)-1
          IOFF=1
          DO 500 I=ISTART,IEND
            CALL PUTLST(PK(IOFF),IDISNUM(I,1),1,2,2,IDISNUM(I,2))
            IOFF=IOFF+IDISSIZ(I)
  500     CONTINUE
        ENDIF
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
  310   CONTINUE
        ISTART=IBKSTRT(IBK)
        IEND=IBKSTRT(IBK+1)-1
        IOFF=1
        DO 400 I=ISTART,IEND
          CALL PUTLST(PK(IOFF),IDISNUM(I,1),1,2,1,IDISNUM(I,2))
          IOFF=IOFF+IDISSIZ(I)
  400   CONTINUE
C
   50 CONTINUE
C
      CLOSE(UNIT=LUSRT1,STATUS='DELETE')
C
      CALL TIMER(1)
      WRITE(LUOUT,5000)TIMENEW
 5000 FORMAT(T3,'@PKCHAN1-I, Chaining of PK lists required ',F10.3,
     &          ' seconds.',/)
C
      RETURN
      END
