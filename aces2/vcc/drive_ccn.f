










C
      SUBROUTINE DRIVE_CCN(ICORE,MAXCOR,IUHF,NCYCLE,ICONTL,
     &                     DAMP_PARAMETER)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION ICORE(MAXCOR)
      DIMENSION ECORR(3)
     
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC,CC2
      LOGICAL DO_HBAR_4LCCSD
C
C A CC2 described as in CPL,243,409-418 (1995) is implemented. 
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



c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end







c This common block contains the IFLAGS and IFLAGS2 arrays for JODA ROUTINES
c ONLY! The reason is that it contains both arrays back-to-back. If the
c preprocessor define MONSTER_FLAGS is set, then the arrays are compressed
c into one large (currently) 600 element long array; otherwise, they are
c split into IFLAGS(100) and IFLAGS2(500).

c iflags(100)  ASVs reserved for Stanton, Gauss, and Co.
c              (Our code is already irrevocably split, why bother anymore?)
c iflags2(500) ASVs for everyone else

      integer        iflags(100), iflags2(500)
      common /flags/ iflags,      iflags2
      save   /flags/





      COMMON /SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC,CC2
      COMMON /ENERGY/ ENERGY(500,2),IXTRLE(500)
C
C Note at the first-iteration, T1 vector is zero and T2 is MBPT(2).
C
      CALL AMPSUM(ICORE(1),MAXCOR,IUHF,0,.TRUE.,'T')
C
C Initializations of various intermediate lists etc.
C
      INEXT = 1
      SING1 = .FALSE.
      DO_HBAR_4LCCSD = .FALSE.

      CALL RNABIJ(ICORE(INEXT),MAXCOR,IUHF,'T')
      CALL SETLST(ICORE(INEXT),MAXCOR,IUHF)

      IF (MOD(IFLAGS(21),2).EQ.0) THEN
         CALL DRRLE(ICORE(I0),MAXCOR,IUHF,RLECYC,.FALSE.)
      END IF
      IF (IFLAGS(21).EQ.1) THEN
         CALL DIISLST(1,IUHF,.TRUE.)
      END IF

      RELCYC = 0
      ICYCLE = 1

 100  CONTINUE
     
      RELCYC = RLECYC + 1

      CALL INITIN(ICORE(INEXT),MAXCOR,IUHF)
      CALL INITSN(ICORE(INEXT),MAXCOR,IUHF)
C
C Generate W(mbej), W(mnij) W(mbej, F(ea), F(ij) and F(me)
C intermediates. CC2 like methods requires that the next to last 
C argument be 2. 
C
      CALL GENINT(ICORE(INEXT),MAXCOR,IUHF,2,DO_HBAR_4LCCSD)
C 
C Evalaute Sum_e T(i,e)F(a,e) contribution to T1.
C
      CALL FEACONT_MODF(ICORE(INEXT),MAXCOR,IUHF)
C
C Evaluate Sum_e T(m,a)F(m,i) contribution to T1.
C
      CALL FMICONT_MODF(ICORE(INEXT),MAXCOR,IUHF)
C
C Evalaute Sum_e T(m,a)F(m,i) contribution to T1.
C  
      IF (SING1) THEN

         CALL FMECONT(ICORE(INEXT),MAXCOR,IUHF,1)
         IF(IUHF.NE.0) CALL FMECONT(ICORE(INEXT),MAXCOR,IUHF,2)
C
C Evaluate -P(ab)Sum_m T(m,a)<mb||ij> contribution to T2.
C
         CALL T1INT2A(ICORE(INEXT),MAXCOR,IUHF)
C
C Evaluate +P(ij)Sum_m T(i,e)<ab||ej> contribution to T2.
C
         CALL T1INT2B(ICORE(INEXT),MAXCOR,IUHF)
C
C Evaluate -Sum_nf T(n,f)<na||if> contribution to T1.
C
         CALL T1INT1(ICORE(INEXT),MAXCOR,IUHF,1)
         IF(IUHF.NE.0) CALL T1INT1(ICORE(INEXT),MAXCOR,IUHF,2)

      ENDIF
C
C
C DO W INTERMEDIATE CONTRIBUTION TO T2 EQUATION
C
C Evaluate the following three contributions to T2:
C
C   -P(ij)P(ab) Sum_me T(i,e)T(m,a) <mb||ej> (T12INT2)
C   +P(ij)P(ab) Sum_me T(i,e)T(m,a) W(mb,ej) (DRRNG)
C   +1/2 Sum_mn Tau(mn,ab)W(mn,ij) +
C    1/2 Sum_ef Tau(ij,ef)W(ab,ef)           (DRLAD)
C
      CALL DRE3EN(ICORE(INEXT),MAXCOR,IUHF,0)
C
C Evaluate the following three contributions to T1: Also,
C denominator weigh the T1 increments
C
C   -1/2 Sum_mef T(im,ef)<ma||ef>  (T2T1AA1, T2T1AB1)
C   -1/2 Sum_men T(nm,ei)<nm||ei>  (T2T1AA2, T2T1AB2)
C
      CALL E4S(ICORE(INEXT),MAXCOR,IUHF,DUMMY)
C 
C Denominator weigh the T2 increments. 
C
      CALL NEWT2(ICORE(INEXT),MAXCOR,IUHF)
C
C At this point, we should have all the t1 contributions 
C computed. So, we can look at the difference between
C the current and previouse T1 vector. 
C 
      CALL AMPSUM(ICORE(INEXT),MAXCOR,IUHF,0,.TRUE.,'T')
      CALL DRTSTS(ICORE(INEXT),MAXCOR,ICYCLE,IUHF,ICONVG,ICONTL,
     &            SING1,0,'T')
C
C Compute correlation energies for each cycle.
C
      CALL CMPENG(ICORE(INEXT),MAXCOR,60,2,ECORR,ENERGY(ICYCLE+1,1),
     &            ENERGY(ICYCLE+1,2),IUHF,1)
C

      IF (IFLAGS(21).EQ.1) THEN
         CALL DODIIS0(ICORE(INEXT),MAXCOR/IINTFP,IUHF,1,ICYCLE,
     &                ICONVG,ICONTL,SING1,44,61,90,0,90,2,70,
     &                '     ',DAMP_PARAMETER)
      END IF

      IF (ICONVG.NE.0) THEN

         CALL DRMOVE(ICORE(I0),MAXCOR,IUHF,0,SING1)
         CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')

      END IF

      IF (ICONVG.EQ.0) THEN

         CALL CMPENG(ICORE(INEXT),MAXCOR,43,0,ECORR,ENERGY(ICYCLE+1,1),
     &               ENERGY(ICYCLE+1,2),IUHF,1)
     
         CALL AMPSUM(ICORE(INEXT),MAXCOR,IUHF,0,SING1,'T')
         IF (IUHF.NE.0) CALL S2PROJ(ICORE(INEXT),MAXCOR,IUHF,SING1)
C
         CALL FINISH(ICYCLE+1)
         CALL DDMPTGSS(ICORE(INEXT), MAXCOR/IINTFP, IUHF, 0, 'TGUESS  ')

         CALL INIT2(IUHF)
         CALL DRMOVE(ICORE(INEXT),MAXCOR,IUHF,100,SING1)
         IF (IUHF.EQ.0.AND..NOT.UCC) CALL RESET(ICORE(INEXT),
     &                                          MAXCOR,IUHF)

         CALL ACES_FIN 
         WRITE(6,1020)
 1020 FORMAT(/,77('-'),/,32X,'Exiting xvcc',/,77('-'),/)

         CALL ACES_EXIT(0)

      ENDIF 

      IF (MOD(IFLAGS(21),2).EQ.0) THEN
         CALL DRRLE(ICORE(INEXT),MAXCOR,IUHF,RLECYC,.FALSE.)
         CALL DRMOVE(ICORE(INEXT),MAXCOR,IUHF,0,SING1)
         CALL RNABIJ(ICORE(INEXT),MAXCOR,IUHF,'T')
      END IF
C
      ICYCLE = ICYCLE +  1
      IF (ICYCLE .LE. NCYCLE) GOTO 100
       
      IF (ICONVG .NE. 0) THEN
         CALL AMPSUM(ICORE(INEXT),MAXCOR,IUHF,0,SING1,'T')
         Write(6,"(a,a)") "The CC2 equations did not converge",
     &                     " in alloted number of cyclces."
         call aces_exit(1)
      ELSE
 
         Write(6,"(A,1x,i3,1x,A)")"T1 iterations converged in",
     &                             ICYCLE, "cyclces."
      ENDIF 
       
C
      RETURN
      END
