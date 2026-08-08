











      SUBROUTINE MAKE_GUESS_DEN(SCR1, SCR2, SCRA, SCRN, DENS, NBAS,
     &                          LDIM1, LDIM2, IUHF, ROHF, ROHFMO,
     &                          GUESS_TYPE)
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
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
      CHARACTER*11 GUESS_TYPE
      LOGICAL ROHF, ROHFMO
C
      DIMENSION SCR1(NBAS*NBAS),SCRN(NBAS*NBAS), SCR2(LDIM2),
     &          SCRA(LDIM2), DENS((IUHF+1)*LDIM1), DOCC(MAXBASFN*2),
     &          NSUM(16)
C
      COMMON /FILES/  LUOUT,MOINTS
      COMMON /FLAGS/  IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /POPUL/  NOCC(16)
C
      IF (GUESS_TYPE.NE."NDDO"        .AND.
     &    GUESS_TYPE.NE."PROJ_STO-3G" .AND.
     &    GUESS_TYPE.NE."EHT"         .AND. 
     &    GUESS_TYPE.NE. "EXPORT_DENS") THEN
C
           CALL GETREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,SCR1)
           CALL MKDP_OCCNOS(DOCC, NOCC, NSUM, NBAS, IUHF, .TRUE.)
           IROFF = 0
           DO 810 I=1,NIRREP
              CALL GETBLK(SCR1,SCR2,NBFIRR(I),NBAS,IREPS(I))
              IOFF = IROFF + 1 
              CALL MKDDEN_4IRREP(SCR2,SCRA,SCRN,NSUM(I),DOCC(IOFF),
     &                           NBFIRR(I),1,IUHF)
              IROFF = IROFF + NBFIRR(I)
              CALL SQUEZ2(SCRA,DENS(ITRIOF(I)),NBFIRR(I))
  810      CONTINUE
C
           IF (IUHF .EQ. 1) THEN
               CALL GETREC(20,'JOBARC','SCFEVCB0',NBAS*NBAS*IINTFP,
     &                     SCR1)
               IROFF = 0
               DO 820 I=1,NIRREP
                  IOFF = IROFF  + NBAS + 1
                  CALL GETBLK(SCR1,SCR2,NBFIRR(I),NBAS,IREPS(I))
                  CALL MKDDEN_4IRREP(SCR2,SCRA,SCRN,NSUM(8+I),
     &                               DOCC(IOFF),NBFIRR(I),2,IUHF)
                  IROFF = IROFF + NBFIRR(I)
                  CALL SQUEZ2(SCRA,DENS(LDIM1+ITRIOF(I)),NBFIRR(I))
  820          CONTINUE
           ENDIF
      ENDIF
C 
      IF (ROHF) THEN
         IF (IFLAGS(1) .GE. 10) WRITE(LUOUT,2050)
 2050    FORMAT(' @INITGES-I, Writing out MOVECTOR record. ')
C
         CALL GETREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,SCR1)
         CALL   ZERO(SCRN,NBAS*NBAS)
         DO 830 I=1,NIRREP
            CALL GETBLK(SCR1,SCRN(ISQROF(I)),NBFIRR(I),NBAS,IREPS(I))
  830    CONTINUE

         cALL PUTREC(20,'JOBARC','MOVECTOR',ISQRLN(NIRREP+1)*IINTFP,
     &               SCRN)
         ROHFMO = .TRUE.
      ENDIF
C
      RETURN
      END

