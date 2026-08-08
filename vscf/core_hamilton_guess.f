











      SUBROUTINE CORE_HAMILTON_GUESS(SCR1, SCRN, ONEH, SCRA,
     &                                SCR2, SCRB, XFORM, SCRTMP,
     &                                LDIM1, LDIM2, IUHF, NBAS)

       IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
c molcas.com : begin
      logical seward, petite_list
      character*8 fnord
      integer luord
      integer ipmat(8,2), lbbt, lbbs
      common /molcas_com/ seward, petite_list, fnord, luord,
     &                    ipmat, lbbt, lbbs
c molcas.com : end
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



C
       DIMENSION SCR1(NBAS*NBAS),SCRN(NBAS*NBAS),ONEH(LDIM1),
     &           XFORM(LDIM1),SCR2(LDIM2),SCRTMP(LDIM2),
     &           SCRA(LDIM2),SCRB(LDIM2)
C
       CALL ZERO(SCR1,NBAS*NBAS)
       CALL ZERO(SCRN,NBAS*NBAS)
C
       DO 110 I=1,NIRREP
       IF(NBFIRR(I).EQ.0) GOTO 110
C
       CALL EXPND2(ONEH(ITRIOF(I)),SCR2,NBFIRR(I))
       CALL EXPND2(XFORM(ITRIOF(I)),SCRA,NBFIRR(I))
       CALL TRANSP(SCRA,SCRB,NBFIRR(I),NBFIRR(I))
       CALL MXM(SCRB,NBFIRR(I),SCR2,NBFIRR(I),SCRTMP,NBFIRR(I))
       CALL MXM(SCRTMP,NBFIRR(I),SCRA,NBFIRR(I),SCR2,NBFIRR(I))
       CALL ZERO(SCRB,NBFIRR(I)*NBFIRR(I))
       CALL EIG(SCR2,SCRB,NBFIRR(I),NBFIRR(I),0)
C
C     SCRB has core Hamiltonian eigenvectors in the "canonical" basis.
C     Transform these to the SO basis.
C
       CALL MXM(SCRA,NBFIRR(I),SCRB,NBFIRR(I),SCR2,NBFIRR(I))
C
C     Core Hamiltonian orbitals in the SO basis for this symmetry block
C     are in SCR2. Copy into SCRA and SCRB (alpha and beta sets).
C
c YAU : old
c      CALL ICOPY(NBFIRR(I)*NBFIRR(I)*IINTFP,SCR2,1,SCRA,1)
c      CALL ICOPY(NBFIRR(I)*NBFIRR(I)*IINTFP,SCR2,1,SCRB,1)
c YAU : new
       CALL DCOPY(NBFIRR(I)*NBFIRR(I),SCR2,1,SCRA,1)
       CALL DCOPY(NBFIRR(I)*NBFIRR(I),SCR2,1,SCRB,1)
c YAU : end
C
       CALL PUTBLK(SCR1,SCRA,NBFIRR(I),NBAS,IREPS(I))
       IF(IUHF.GT.0)THEN

        CALL PUTBLK(SCRN,SCRB,NBFIRR(I),NBAS,IREPS(I))
       ENDIF
C
  110  CONTINUE
C
       CALL PUTREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,SCR1)
       CALL PUTREC(20,'JOBARC','GUESSVA0',NBAS*NBAS*IINTFP,SCR1)
       IF(IUHF.GT.0)THEN
        CALL PUTREC(20,'JOBARC','SCFEVCB0',NBAS*NBAS*IINTFP,SCRN)
        CALL PUTREC(20,'JOBARC','GUESSVB0',NBAS*NBAS*IINTFP,SCRN)
       ENDIF
C
      RETURN
      END

