











      SUBROUTINE SWITCH_ORBS(IRRP,IBGN,IEND,NBAS,EVAL,EVEC,LDIM2,
     &                       ISPN,IOS)

      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION EVAL(1),EVEC(1),SCR(1),ZJUNK
      DIMENSION ILOCATE(1)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /POPUL/ NOCC(8,2)

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
      Write(6,"(a,4I6)") "ISPN,IBGN,IEND,IOS:",ISPN,IBGN,IEND,IOS
      DO ISPIN=1,2

         IF(IOS.EQ.0.AND.ISPIN.EQ.ISPN) THEN

         IPOS1=(ISPIN-1)*LDIM2+ISQROF(IRRP)-1+INDX2(1,IBGN,NBFIRR(IRRP))
         IPOS2=(ISPIN-1)*LDIM2+ISQROF(IRRP)-1+INDX2(1,IEND,NBFIRR(IRRP))

      Write(6,"(a,2I6)") "IPOS1 and IPOS2   :", IPOS1,IPOS2
         CALL DSWAP(NBFIRR(IRRP),EVEC(IPOS1),1,EVEC(IPOS2),1)
C
         IOFF1=IRPOFF(IRRP)+IBGN
         IOFF2=IRPOFF(IRRP)+IEND
         ZJUNK=EVAL((ISPIN-1)*NBAS+IOFF1)
         EVAL((ISPIN-1)*NBAS+IOFF1)=EVAL((ISPIN-1)*NBAS+IOFF2)
         EVAL((ISPIN-1)*NBAS+IOFF2)=ZJUNK

         ENDIF 
C
      ENDDO

      RETURN
      END
