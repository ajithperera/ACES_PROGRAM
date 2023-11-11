











      SUBROUTINE EIGVECTORS_FRM_JOBARC(SCR1, SCRN, SCRA, SCRB,
     &                                 SCRTMP, LDIM2, NBAS,
     &                                 IUHF)

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
      DIMENSION SCR1(NBAS*NBAS), SCRN(NBAS*NBAS), SCRA(LDIM2),
     &          SCRB(LDIM2),SCRTMP(LDIM2)
C
      COMMON /FLAGS/IFLAGS(100)
      DATA ZILCH /0.0/
C
C Not Bruckner and HF Stability=on
C
CSSS      IF(IFLAGS(22) .EQ. 0 .AND. IFLAGS(74) .LT. 2)THEN
C I am changing this to IFLAGS(74) .EQ. 1. The SCFEVCAS and
C SCFEVCA0 have same data (within the confines of SCF). The 
C decision has to be made  whether we are going to read 
C SCFEVCA0 or SCFEVCAS. The diffrences is that it is possible
C that SCFEVCA0 can be reorderd (Brueckner?). The
C HFSATBILTY {ON, FOLLOW} reorder the vectors and so 
C read the SCFEVCA0. Also, allow the external program to
C write SCFEVCA0 because that is the record that people are
C familiar with (Localization programs). My feeling is that
C it is better we always read SCFEVCA0. Ajith Perera, 04/2014.
C
      IF(IFLAGS(22) .EQ. 0 .AND. IFLAGS(74) .LT. 2) THEN
C
        ZTEST = 0.0D+00
        CALL ZERO(SCR1,NBAS*NBAS)
        CALL ZERO(SCRN,NBAS*NBAS)
        CALL GETREC(-1,'JOBARC','SCFEVCAS',NBAS*NBAS*IINTFP,SCRN)
        ZTEST=SNRM2(NBAS*NBAS,SCRN,1)
        IF(ZTEST.NE.ZILCH) THEN
         CALL SORTHO(SCRN,SCR1,SCRA,SCRB,SCRTMP,LDIM2,NBAS)
        ENDIF
        CALL PUTREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,SCRN)
        CALL ZERO(SCR1,NBAS*NBAS)
        CALL ZERO(SCRN,NBAS*NBAS)
C
        IF(IUHF.EQ.1) THEN
         ZTEST = 0.0D+00
         CALL ZERO(SCR1,NBAS*NBAS)
         CALL ZERO(SCRN,NBAS*NBAS)
         CALL GETREC(-1,'JOBARC','SCFEVCBS',NBAS*NBAS*IINTFP,SCRN)
         ZTEST=SNRM2(NBAS*NBAS,SCRN,1)
         IF(ZTEST.NE.ZILCH) THEN
          CALL SORTHO(SCRN,SCR1,SCRA,SCRB,SCRTMP,LDIM2,NBAS)
         ENDIF
         CALL PUTREC(20,'JOBARC','SCFEVCB0',NBAS*NBAS*IINTFP,SCRN)
         CALL ZERO(SCR1,NBAS*NBAS)
         CALL ZERO(SCRN,NBAS*NBAS)
        ENDIF
C
       ELSE

        ZTEST = 0.0D+00
        CALL ZERO(SCR1,NBAS*NBAS)
        CALL ZERO(SCRN,NBAS*NBAS)
        CALL GETREC(-1,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,SCRN)
        ZTEST=SNRM2(NBAS*NBAS,SCRN,1)

        Write(6,"(2x,a)") "The Alpha eigenvectors read from JOBARC"
        call checksum("SCF_VECA", SCRN, NBAS*NBAS)
C        call output(Scrn, 1, Nbas, 1, Nbas, Nbas, Nbas, 1)
        Write(*,"(A,F10.6)") "NORM-CHECK =", ZTEST

        IF(ZTEST.NE.ZILCH) THEN
           CALL REOREV(SCRN,SCR1,NBAS,1)
           CALL SORTHO(SCR1,SCRN,SCRA,SCRB,SCRTMP,LDIM2,NBAS)
        ENDIF
        CALL PUTREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,SCR1)
        CALL ZERO(SCR1,NBAS*NBAS)
        CALL ZERO(SCRN,NBAS*NBAS)
C
        IF(IUHF.EQ.1) THEN
         ZTEST = 0.0D+00
         CALL ZERO(SCR1,NBAS*NBAS)
         CALL ZERO(SCRN,NBAS*NBAS)
         CALL GETREC(-1,'JOBARC','SCFEVCB0',NBAS*NBAS*IINTFP,SCRN)
         ZTEST=SNRM2(NBAS*NBAS,SCRN,1)
         IF(ZTEST.NE.ZILCH) THEN
         CALL REOREV(SCRN,SCR1,NBAS,2)
         CALL SORTHO(SCR1,SCRN,SCRA,SCRB,SCRTMP,LDIM2,NBAS)
         ENDIF
         CALL PUTREC(20,'JOBARC','SCFEVCB0',NBAS*NBAS*IINTFP,SCR1)
         CALL ZERO(SCR1,NBAS*NBAS)
         CALL ZERO(SCRN,NBAS*NBAS)
        ENDIF
      ENDIF

      RETURN
      END

