










      SUBROUTINE MKEVAL(XFORM,FOCK,EVAL,EVEC,SCR1,SCR2,NBAS,LDIM1,
     &                  LDIM2,LDIM3)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION XFORM(LDIM1)
      DIMENSION FOCK(2*LDIM1),EVEC(2*LDIM2),EVAL(2*NBAS)
      DIMENSION SCR1(LDIM3),SCR2(LDIM3)
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
      DATA ONE /1.0/
      DATA ZILCH /0.0/
C
      INDX2(I,J,N)=I+(J-1)*N
C
C  First we must convert the eigenvectors to the AO basis.
C
      DO 50 I=1,NIRREP
        IF(NBFIRR(I).EQ.0) GOTO 50
        CALL EXPND2(XFORM(ITRIOF(I)),SCR1,NBFIRR(I))
        CALL MXM(SCR1,NBFIRR(I),EVEC(ISQROF(I)),NBFIRR(I),
     &           SCR2,NBFIRR(I))
c YAU : old
c       CALL ICOPY(ISQRLN(I)*IINTFP,SCR2,1,EVEC(ISQROF(I)),1)
c YAU : new
        CALL DCOPY(ISQRLN(I),SCR2,1,EVEC(ISQROF(I)),1)
c YAU : end
   50 CONTINUE
C
C
      DO 100 ISPIN=1,2
        DO 110 I=1,NIRREP
          IF(NBFIRR(I).EQ.0) GOTO 110
          CALL EXPND2(FOCK((ISPIN-1)*LDIM1+ITRIOF(I)),SCR1,NBFIRR(I))
          CALL XGEMM('T','N',NBFIRR(I),NBFIRR(I),NBFIRR(I),ONE,
     &               EVEC(ISQROF(I)),NBFIRR(I),SCR1,NBFIRR(I),ZILCH,
     &               SCR2,NBFIRR(I))
          CALL XGEMM('N','N',NBFIRR(I),NBFIRR(I),NBFIRR(I),ONE,
     &               SCR2,NBFIRR(I),EVEC(ISQROF(I)),NBFIRR(I),ZILCH,
     &               SCR1,NBFIRR(I))
          DO 111 J=1,NBFIRR(I)
            EVAL((ISPIN-1)*NBAS+IRPOFF(I)+J)=SCR1(INDX2(J,J,NBFIRR(I)))
  111     CONTINUE
  110   CONTINUE
  100 CONTINUE
C
C  Copy alpha eigenvectors to beta eigenvectors.
c YAU : old
c     CALL ICOPY(ISQRLN(NIRREP+1)*IINTFP,EVEC,1,EVEC(ISQRLN(NIRREP+1)+1),1)
c YAU : new
      CALL DCOPY(ISQRLN(NIRREP+1),EVEC,1,EVEC(ISQRLN(NIRREP+1)+1),1)
c YAU : end
C
      RETURN
      END
