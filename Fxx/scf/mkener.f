










      SUBROUTINE MKENER(ONEH,DENS,FOCK,SCR1,SCR2,SCR3,LDIM1,LDIM2,
     &                  REPULS,DMAX,ITER,IUHF,ETOT,IQRHF,NOCONV,scfks,
     &                  scfkslastiter,EDIFF)
C
C  The energy is caculated using the formula
C
C      E(0) = 0.5*SUM  SUM  [DENS     *(ONEH      + FOCK     )]
C                    mu   nu     mu,nu      mu,nu       mu,nu
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      parameter(luint=10)
c
      logical rohf,scfks,scfkslastiter,noconv
      DIMENSION ONEH(LDIM1),FOCK((IUHF+1)*LDIM1),DENS((IUHF+1)*LDIM1)
      DIMENSION SCR1(LDIM1),SCR2(LDIM2),SCR3(LDIM2)
C
      COMMON /FILES/ LUOUT,MOINTS
      common /flags/ iflags(100)
      common /machsp/ iintln,ifltln,iintfp,ialone,ibitwd
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
      DATA HALF /0.5/
      DATA ONE /1.0/
C
      INDX2(I,J,N)=I+(J-1)*N
C
      rohf=.false.
      if(iflags(11).eq.2) rohf=.true.
C
      if(scfks .and. .not. scfkslastiter) then
        if(scfks .and. iter.gt.1) then
          call getrec(1,'JOBARC','KSSCFENG',iintfp,etot)
         goto 7000
        end if
      end if
      CALL ZERO(SCR2,LDIM2)
      CALL ZERO(SCR3,LDIM2)
      EELEC=0.0
      DO 40 ISPIN=1,(IUHF+1)
        CALL ZERO(SCR1,LDIM1)
        CALL SAXPY(LDIM1,ONE,ONEH,1,SCR1,1)
        CALL SAXPY(LDIM1,ONE,FOCK(((ISPIN-1)*LDIM1)+1),1,SCR1,1)
        DO 50 I=1,NIRREP
          CALL EXPND2(DENS(((ISPIN-1)*LDIM1)+ITRIOF(I)),SCR2,NBFIRR(I))
          CALL EXPND2(SCR1(ITRIOF(I)),SCR3,NBFIRR(I))
C
C  Sum up the contribution from the one-electron Hamiltonian to the
C  electronic energy.
C
          DO 100 J=1,NBFIRR(I)
            DO 110 K=1,NBFIRR(I)
              EELEC=EELEC+HALF*SCR2(INDX2(J,K,NBFIRR(I)))*
     &              SCR3(INDX2(J,K,NBFIRR(I)))
  110       CONTINUE
  100     CONTINUE
C
   50   CONTINUE
   40 CONTINUE
C
      ETOT=EELEC+REPULS
 7000  continue 
      IF(IQRHF.EQ.1) GOTO 8000
c     IF(DMAX.GT.TOL.OR.ITER.EQ.0.or.rohf) THEN
      IF(NOCONV.OR.ITER.EQ.0) THEN
        WRITE(LUOUT,5000)ITER+1,ETOT,DMAX
 5000   FORMAT(T3,I6,5X,F20.10,10X,D20.10)
      ELSE
        WRITE(LUOUT,5002)ETOT
 5002   FORMAT(T3,3X,'E(SCF)= ',F20.10/)
      ENDIF

       IF (ITER .EQ .2) THEN
          CALL PUTREC(20,'JOBARC','EOLD   ',IINTFP,ETOT)
       ENDIF

       IF (ITER .GT.2) THEN
          CALL GETREC(20,'JOBARC','EOLD   ',IINTFP,EOLD)
          CALL PUTREC(20,'JOBARC','EOLD   ',IINTFP,ETOT)
          Ediff = DABS(EOLD-ETOT)
       ENDIF 

      if(iflags(16).eq.0) then
         write(luout,5010)
 5010    format(/)
         call putrec(20,'JOBARC','SCFENEG ',iintfp,etot)
cjp
cjp avoid this, will be neede when oofcor is active
cjp         close(moints,status='delete')
cjp
c         close(luint,status='keep')
c         call endscf
      endif

 8000 RETURN
      END
