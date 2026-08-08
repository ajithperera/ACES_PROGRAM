










      SUBROUTINE FAO2MOROHF(FOCK,OVLP,EVEC,SCR2,SCR3,NBAS,BACK)
C
C TRANSFORMS THE FOCK MATRIX TO THE MOLECULAR ORBITAL BASIS
C
CEND
      IMPLICIT INTEGER (A-Z)
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
      DOUBLE PRECISION FOCK,EVEC,ONE,ZILCH,SCR2,SCR3,OVLP
      LOGICAL BACK
      DIMENSION FOCK(1),EVEC(1),SCR2(1),SCR3(1),OVLP(1)
      DATA ONE /1.0/
      DATA ZILCH /0.0/
      IOFF1=1
      IOFF2=1
      IOFF3=1
      DO 10 IRREP=1,NIRREP
       NSIZ=NBFIRR(IRREP)
       if (nsiz.ne.0) then
         CALL EXPND2(FOCK(IOFF1),SCR2,NSIZ)
         IF (BACK) THEN
           CALL XGEMM('N','T',NSIZ,NSIZ,NSIZ,ONE,SCR2,NSIZ,
     &                EVEC(IOFF2),NSIZ,ZILCH,SCR3,NSIZ)
           CALL XGEMM('N','N',NSIZ,NSIZ,NSIZ,ONE,EVEC(IOFF2),
     &                NSIZ,SCR3,NSIZ,ZILCH,SCR2,NSIZ)
           call xgemm('n','n',NSIZ,NSIZ,NSIZ,ONE,SCR2,NSIZ,
     &                OVLP(IOFF3),NBAS,ZILCH,SCR3,NSIZ)
           call xgemm('n','n',NSIZ,NSIZ,NSIZ,ONE,OVLP(IOFF3),
     &                NBAS,SCR3,NSIZ,ZILCH,SCR2,NSIZ)
         ELSE
           CALL XGEMM('N','N',NSIZ,NSIZ,NSIZ,ONE,SCR2,NSIZ,
     &                EVEC(IOFF2),NSIZ,ZILCH,SCR3,NSIZ)
           CALL XGEMM('T','N',NSIZ,NSIZ,NSIZ,ONE,EVEC(IOFF2),
     &                NSIZ,SCR3,NSIZ,ZILCH,SCR2,NSIZ)
         ENDIF
         CALL SQUEZ2(SCR2,FOCK(IOFF1),NSIZ)
         IOFF1=IOFF1+(NSIZ*(NSIZ+1))/2
         IOFF2=IOFF2+NSIZ*NSIZ
         IOFF3=IOFF3+NSIZ*(NBAS+1)
        endif
10    CONTINUE
      RETURN
      END
