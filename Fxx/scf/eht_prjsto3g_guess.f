











      SUBROUTINE EHT_PRJSTO3G_GUESS(DENS, CZAO, CSO, ZAOSO,
     &                              SCR1, SCRN, DCORE, IUHF,
     &                              NBAS, NBASX, LDIM1, MAXDCOR)

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      LOGICAL MBS
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
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
C
       DIMENSION DCORE(MAXDCOR),DENS((IUHF+1)*LDIM1),SCR1(NBAS*NBAS),
     &           SCRN(NBAS*NBAS), CZAO(NBASX*NBAS,(IUHF+1)),
     &           ZAOSO(NBAS*NBASX)

C-----------------------------------------------------------------------
C     JDW, Summer 2004. This is a fudge to enable two guess options.
C
C     IFLAGS(45) = 6 The program will read a file called MOFILE to
C                    get MOs in a minimal basis set and then generate
C                    an initial guess for another basis set. This has
C                    only been coded for certain basis sets. See the
C                    routine GETDENS3 and the SPLBAS file. MOFILE must
C                    have been created previously. MOFILE is identical
C                    to AOBASMOS. IFLAGS(45) = 6 is synonymous with the
C                    2-step procedure.
C
C     IFLAGS(45) = 7 The program will generate an EHT guess and write
C                    a file called EHTMOFILE. If this is NOT a minimal
C                    basis set calculation, GETDENS3 will then try to
C                    create an initial guess for the actual basis set
C                    being used. This is the one-step procedure. Alter-
C                    natively, if this option is run with a minimal
C                    basis set, GETDENS3 will be bypassed and an STO-3G
C                    initial density generated from the EHT coeff-
C                    icients. An STO-3G scf calculation will then be
C                    performed, created AOBASMOS, which can be renamed
C                    to MOFILE for the second job of the two-step proc-
C                    edure.
C-----------------------------------------------------------------------
C     Projecting an STO-3G set of MOs into current basis set.
C     Limitations for ROHF right now --- only creates occupied MOs and
C     and alpha and beta density matrices.
C-----------------------------------------------------------------------
C
      IF(IFLAGS(45).EQ.6)THEN
       call getdens3(dens,czao,cso,zaoso,scr1,scrn,
     &               nbasx,nbas,natoms,
     &               nirrep,ldim1,itriof,ireps,nbfirr,iuhf,2)
      ENDIF
C
C-----------------------------------------------------------------------
C     Extended Huckel Guess. Only available for H-Cl and STO-3G.
C     Limitations for ROHF right now --- only creates occupied MOs and
C     and alpha and beta density matrices.
C-----------------------------------------------------------------------
C
      IF(IFLAGS(45).EQ.7)THEN
       CALL EHTDIM(NBASEHT,NBASXEHT)
       CALL CALC_S(DCORE,DCORE(1+NBASXEHT*NBASXEHT),
     &                   DCORE(1+NBASXEHT*NBASXEHT+3*NATOMS),
     &                   DCORE(1+NBASXEHT*NBASXEHT+3*NATOMS),
     &             NBASEHT,NBASXEHT,NATOMS)
C
C     CALC_S puts overlap integrals in DCORE. Length: NBASXEHT*NBASXEHT.
C
cccc This has potential      CALL SYMTRANS
       CALL EHTGSS(DENS,DCORE(1+NBASX*NBASX),DCORE,
     &             DCORE(1+2*NBASX*NBASX),DCORE(1+3*NBASX*NBASX),SCR1,
     &             EVAL,NBASEHT,NBASXEHT,LDIM1,LDIM2,IUHF)
C
       MBS = (NBASX.EQ.NBASXEHT)
C
C     Read information from EHTMOFILE to create initial guess in
C     current basis set.
C
C     EHTGSS prepares EHTMOFILE. MOFILE has to be manually prepared:
C     it is the same as AOBASMOS for a minimal basis set calculation.
C
       IF(.NOT.MBS)THEN
        call getdens3(dens,czao,cso,zaoso,scr1,scrn,
     &                nbasx,nbas,natoms,
     &                nirrep,ldim1,itriof,ireps,nbfirr,iuhf,1)
       ENDIF
      ENDIF

      Return
      End

