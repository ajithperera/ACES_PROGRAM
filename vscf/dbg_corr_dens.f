










      SUBROUTINE DBG_CORR_DENS(DENS, SCRN, SCR1, SCR2, 
     &                         LDIM1, NBAS, NBASX, IUHF)

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
      DIMENSION DENS((IUHF+1)*LDIM1), SCR1(NBAS,NBAS),
     &          SCRN(NBAS,NBAS), SCR2(NBAS,NBAS)
C
      CALL ZERO(DENS, (IUHF+1)*LDIM1)
      CALL GETREC(20,'JOBARC','RELDENSA',NBAS*NBAS*IINTFP,SCRN)
      CALL GETREC(20,'JOBARC','SCFEVECA',NBAS*NBAS*IINTFP,SCR1)
      CALL XGEMM("N", "N", NBAS, NBAS, NBAS, 1.0D0, SCR1, NBAS,
     &            SCRN, NBAS, 0.0D0, SCR2, NBAS)
      CALL XGEMM("N", "T", NBAS, NBAS, NBAS, 1.0D0, SCR2, NBAS,
     &            SCR2, NBAS, 0.0D0, SCRN, NBAS)
C
      DO I=1,NIRREP
          CALL SQUEZ2(SCRN,DENS(ITRIOF(I)),NBFIRR(I))
      ENDDO
C
      IF (IUHF.EQ.1) THEN
         CALL GETREC(20,'JOBARC','RELDENSB',NBAS*NBAS*IINTFP,SCRN)
         CALL GETREC(20,'JOBARC','SCFEVECB',NBAS*NBAS*IINTFP,SCR1)
         CALL XGEMM("N", "N", NBAS, NBAS, NBAS, 1.0D0, SCR1, NBAS,
     &            SCRN, NBAS, 0.0D0, SCR2, NBAS)
         CALL XGEMM("N", "T", NBAS, NBAS, NBAS, 1.0D0, SCR2, NBAS,
     &            SCR1, NBAS, 0.0D0, SCRN, NBAS)
C
         IOFF = 1
         DO I=1,NIRREP
            CALL SQUEZ2(SCRN,DENS(LDIM1+ITRIOF(I)),NBFIRR(I))
         ENDDO  
C
      ENDIF

      RETURN
      END
