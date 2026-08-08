










      SUBROUTINE SRTALC(IBKNUM,IDISSIZ,IDISNUM,IBKOF,IBKSTRT,IBKDIS,
     &                  MAXBUK,NBUCK,MEMSIZ,IBKSIZ)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION IBKNUM(IBKDIS),IDISSIZ(IBKDIS)
      DIMENSION IBKOF(MAXBUK),IBKSTRT(MAXBUK)
      DIMENSION IDISNUM(IBKDIS,2)
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
C
      INDX(I,J)=J+(I*(I-1))/2
C
      CALL IZERO(IBKNUM,IBKDIS)
      CALL IZERO(IDISSIZ,IBKDIS)
      CALL IZERO(IBKOF,MAXBUK)
      MEMSZ=MEMSIZ/IINTFP
      ISUM=0
      ICNT2=0
      NUMDIS=0
      NBUCK=1
      IBKSIZ=0
      IBKOF(1)=0
      IBKSTRT(1)=1
      DO 10 I=1,NIRREP
        DO 20 J=1,I
          IF(I.EQ.J) THEN
            ITOP=INDX(NBFIRR(I),NBFIRR(I))
            DO 30 K=1,ITOP
              ISUM=ISUM+IRPDS1(INDX(I,I))
              NUMDIS=NUMDIS+1
              IDISSIZ(NUMDIS)=IRPDS1(INDX(I,I))
              IF(ISUM.GT.MEMSZ) THEN
                NBUCK=NBUCK+1
C
                IF(NBUCK.LE.MAXBUK) THEN
                  IBKOF(NBUCK)=IBKOF(NBUCK-1)+ISUM-
     &                         IRPDS1(INDX(I,I))
                  IBKSIZ=MAX(IBKSIZ,ISUM-IRPDS1(INDX(I,I)))
                  IBKSTRT(NBUCK)=NUMDIS
                ENDIF
C
                ISUM=IRPDS1(INDX(I,I))
                ICNT2=ICNT2+1
                IBKNUM(ICNT2)=NBUCK
              ELSE
                ICNT2=ICNT2+1
                IBKNUM(ICNT2)=NBUCK
              ENDIF
   30       CONTINUE
          ELSE
            ITOP=INDX(NBFIRR(J),NBFIRR(J))
            DO 31 K=1,ITOP
              ISUM=ISUM+IRPDS2(2*INDX(I-1,J)-1)
              NUMDIS=NUMDIS+1
              IDISSIZ(NUMDIS)=IRPDS2(2*INDX(I-1,J)-1)
              IF(ISUM.GT.MEMSZ) THEN
                NBUCK=NBUCK+1
C
                IF(NBUCK.LE.MAXBUK) THEN
                  IBKOF(NBUCK)=IBKOF(NBUCK-1)+ISUM-
     &                         IRPDS2(2*INDX(I-1,J)-1)
                  IBKSIZ=MAX(IBKSIZ,ISUM-IRPDS2(2*INDX(I-1,J)-1))
                  IBKSTRT(NBUCK)=NUMDIS
                ENDIF
C
                ISUM=IRPDS2(2*INDX(I-1,J)-1)
                ICNT2=ICNT2+1
                IBKNUM(ICNT2)=NBUCK
              ELSE
                ICNT2=ICNT2+1
                IBKNUM(ICNT2)=NBUCK
              ENDIF
   31       CONTINUE
          ENDIF
   20   CONTINUE
   10 CONTINUE
      IF(IBKSIZ.EQ.0) THEN
        IBKSIZ=ISUM
      ENDIF
      IBKSTRT(NBUCK+1)=NUMDIS+1
C
      ICNT=0
      DO 101 I=1,NIRREP
        DO 102 J=1,I
          DO 103 K=1,NBFIRR(J)
            DO 104 L=1,K
              ICNT=ICNT+1
              IDISNUM(ICNT,1)=INDX(K,L)
              IF(I.EQ.J) THEN
                IDISNUM(ICNT,2)=INDX(I,J)
              ELSE
                IDISNUM(ICNT,2)=36+INDX(I-1,J)
              ENDIF
  104       CONTINUE
  103     CONTINUE
  102   CONTINUE
  101 CONTINUE
C
      IF(NBUCK.GT.MAXBUK) THEN
        WRITE(LUOUT,9000)NBUCK,MAXBUK
 9000   FORMAT(T3,'@SRTALC-F, Insufficient number of buckets for ',
     &            'PK sort.',/,
     &         T8,'Increase variable MAXBUK in subroutine DRVSRT.',/,
     &         T8,' Buckets required: ',I10,/,
     &         T8,'Buckets available: ',I10,/)
        CALL ERREX
      ENDIF
C
      RETURN
      END
