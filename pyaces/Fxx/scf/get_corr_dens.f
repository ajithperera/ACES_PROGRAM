










      SUBROUTINE GET_CORR_DENS(DENS, SCR1, SCR2, SCRA, SCRB, SCRTMP,
     &                         SCRN, SCR3, SCR, LDIM1, LDIM2, NBAS, 
     &                         NBASX, IUHF)

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C


c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end



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
      DIMENSION DENS((IUHF+1)*LDIM1), SCR1(NBAS,NBAS), SCR2(LDIM2),
     &          SCR3(NBAS,NBAS), SCRA(LDIM2), SCRB(LDIM2),
     &          SCRN(NBAS,NBAS), SCRTMP(LDIM2), SCR(NBAS,NBAS)
      DATA ZILCH, ONE / 0.0D0, 1.0D0/
C
      ZTEST = 0.0D+00
      CALL ZERO(SCR1,NBAS*NBAS)
      CALL ZERO(SCRN,NBAS*NBAS)
C
      CALL GETREC(20,'JOBARC','RELDENSA',NBAS*NBAS*IINTFP,SCR1)
C
      Write(6,*) "The Alpha Density read"
      call output(SCR1, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
C
      CALL EIG(SCR1,SCRN,1,NBAS,-1)
C
C
      CALL GETREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,SCR3)
      CALL XGEMM('N', 'N', NBAS, NBAS, NBAS, ONE, SCR3, NBAS,
     &           SCRN, NBAS, ZILCH, SCR, NBAS)

      CALL PUTREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,SCR)

      Write(6,*) "The occupation number"
      call output(SCR1, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
      CALL MKDDEN_FULL(SCRN,SCR2,SCR3,SCR1,NBAS,IUHF)
C
      Write(6,*) "The Alpha Density after MKDDEN"
      call output(SCR3, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
C
      CALL MO2AO2(SCR3, SCRN, SCR1, SCR, NBAS, 1)
C
      DO I=1,NIRREP
          CALL GETBLK(SCRN,SCR2,NBFIRR(I),NBAS,IREPS(I))
          CALL SQUEZ2(SCR2,DENS(ITRIOF(I)),NBFIRR(I))
      ENDDO
C
      IF (IUHF.EQ.1) THEN
         ZTEST = 0.0D+00
         CALL ZERO(SCR1,NBAS*NBAS)
         CALL ZERO(SCRN,NBAS*NBAS)
         CALL GETREC(20,'JOBARC','RELDENSB',NBAS*NBAS*IINTFP,SCR1)
         CALL EIG(SCR1,SCRN,1,NBAS,-1)
C
C
      CALL GETREC(20,'JOBARC','SCFEVCB0',NBAS*NBAS*IINTFP,SCR3)
      CALL XGEMM('N', 'N', NBAS, NBAS, NBAS, ONE, SCR3, NBAS,
     &           SCRN, NBAS, ZILCH, SCR, NBAS)

      CALL PUTREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,SCR)


         CALL MKDDEN_FULL(SCRN,SCR2,SCR3,SCR1,NBAS,IUHF)
C
         CALL MO2AO2(SCR3, SCRN, SCR1, SCR, NBAS, 2)
C
         DO I=1,NIRREP
            CALL GETBLK(SCRN,SCR2,NBFIRR(I),NBAS,IREPS(I))
            CALL SQUEZ2(SCR2,DENS(LDIM1+ITRIOF(I)),NBFIRR(I))
         ENDDO  
C
      ENDIF

      RETURN
      END
