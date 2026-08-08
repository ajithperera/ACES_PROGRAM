











      SUBROUTINE QRHFADD(IRRP,ILOC,NBAS,EVAL,EVEC,IPOSABS,LDIM2,
     &                   ILOCATE,SCR,ISPN,IOS)
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION EVAL(1),EVEC(1),SCR(1),ZJUNK
      DIMENSION ILOCATE(1)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS/ IFLAGS(100)
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
      ISPIN=1
C
      IF(NBFIRR(IRRP)-NOCC(IRRP,ISPN).EQ.0) THEN
       WRITE(LUOUT,9980)IRRP,ISPN
9980   FORMAT(T3,'@QRHFADD-F, Irrep ',I1,' od spin ',I1,
     &           ' has no remaining virtual ',
     &           'orbitals.')
       CALL ERREX
      ENDIF
      WRITE(LUOUT,9987)ILOCATE(IPOSABS),ISPN,EVAL(IPOSABS)
9987  FORMAT(T3,'@QRHFADD-I, Electron added to orbital ',I3,
     &          ' of spin ',I2,/,
     &          T14,' with eigenvalue ',F10.5,'.') 
C
C CALCULATE ILOC ASSUMING THAT SOME EIGENVALUES HAVE BEEN REORDERED
C
      IPOSREL=INEWVC(IPOSABS)
      ILOC=IPOSREL-NOCC(IRRP,ISPN)
C
C  SWAP EIGENVECTORS IF ILOC.NE.1.  WE NEED TO SWITCH EIGENVECTORS IN
C  BOTH THE ALPHA AND BETA BLOCKS SO THAT THE TRANSFORMATION AND OTHER
C  THINGS ARE CORRECT.
C
      IF(ILOC.GT.1)THEN
       IF(ILOC.GT.NBFIRR(IRRP)-NOCC(IRRP,ISPN))THEN
        WRITE(LUOUT,9982)IRRP,NBFIRR(IRRP)-NOCC(IRRP,ISPN),ILOC
9982    FORMAT(T3,'@QRHFADD-F, Irrep ',I1,' of spin ',I1,
     &            ' has only ',I3,
     &            ' remaining virtual orbitals.',/,
     &          T3,'You asked to populate a non-existent ',
     &             'orbital (',I3,').')
        CALL ERREX
       ENDIF
       IOFF=NOCC(IRRP,ISPN)+ILOC
       INOC=NOCC(IRRP,ISPN)+1
       DO 10 ISPIN=1,2
        IF(IOS.NE.0.AND.ISPIN.NE.ISPN) GOTO 10
        IPOS1=(ISPIN-1)*LDIM2+ISQROF(IRRP)-1+INDX2(1,INOC,NBFIRR(IRRP))
        IPOS2=(ISPIN-1)*LDIM2+ISQROF(IRRP)-1+INDX2(1,IOFF,NBFIRR(IRRP))
c YAU : old
c       CALL ICOPY(NBFIRR(IRRP)*IINTFP,EVEC(IPOS1),1,SCR,1)
c       CALL ICOPY(NBFIRR(IRRP)*IINTFP,EVEC(IPOS2),1,EVEC(IPOS1),1)
c       CALL ICOPY(NBFIRR(IRRP)*IINTFP,SCR,1,EVEC(IPOS2),1)
c YAU : new
        CALL DSWAP(NBFIRR(IRRP),EVEC(IPOS1),1,EVEC(IPOS2),1)
c YAU : end
C
C   Now switch around the eigenvalues.
C
        IOFF1=IRPOFF(IRRP)+NOCC(IRRP,ISPN)+ILOC
        IOFF2=IRPOFF(IRRP)+NOCC(IRRP,ISPN)
        ZJUNK=EVAL((ISPIN-1)*NBAS+IOFF1)
        EVAL((ISPIN-1)*NBAS+IOFF1)=EVAL((ISPIN-1)*NBAS+IOFF2)
        EVAL((ISPIN-1)*NBAS+IOFF2)=ZJUNK
        ITMP=ILOCATE(IOFF1)
        ILOCATE(IOFF1)=ILOCATE(IOFF2)
        ILOCATE(IOFF2)=ITMP
C
   10  CONTINUE
      ENDIF
      NOCC(IRRP,ISPN)=NOCC(IRRP,ISPN)+1
C
      RETURN
      END
