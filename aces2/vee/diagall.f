      SUBROUTINE DIAGALL(SCR,MAXCOR,IUHF)
C
C THIS ROUTINE (WHICH SHOULD NOT BE USED ROUTINELY) EXPLICITLY
C CALCULATES THE EFFECTIVE HAMILTONIAN MATRIX ELEMENTS AND THEN
C PROCEEDS TO DIAGONALIZE THE WHOLE MONSTER.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ZILCH,FSTR
      DOUBLE PRECISION SCR,EVAL,SDOT,X,Z0
      DOUBLE PRECISION POLTOT, R
      DOUBLE PRECISION EVALIMX,RESMAX,RESID,DNRM2
      LOGICAL ESPROP,CORE_SEARCH,DOUBLE_EES
      DIMENSION SCR(MAXCOR),IOFFX(2)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/EXTINF/NDIMR,IOLDEST
      COMMON/EXTINF3/IROOT,LOCROOT,ITROOT
      COMMON/RMAT/ R(10000)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREP0(255,2),DIRPRD(8,8)
      COMMON/CALCINFO/NROOT(8)
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
      COMMON/PROPGRAD/ESPROP,IOPTROOT,IOPTSYM
      COMMON/POLAR/POLTOT(3,3)
      COMMON/FLAGS/IFLAGS(100)
C
      DATA ONE,ZILCH /1.0D0,0.0D0/
C
      CALL ZERO(POLTOT,9)
      FSTR        = ZILCH
      ESPROP      = .TRUE.
C
C CORE_SEARCH/DOUBLE_EES only switch on TDENS's beyond-dipole moments;
C not wanted here. (The original CORE_SEACRH typo left CORE_SEARCH an
C undefined implicit INTEGER, and the TDENS calls below were missing
C the DOUBLE_EES argument added to TDENS after this routine was last
C used.)
C
      CORE_SEARCH = .FALSE.
      DOUBLE_EES  = .FALSE.
C 
      WRITE(6,6000)
C
C Full left/right MO transition densities of every root are written
C by TDENS to eom_tdm.dat (EE_ROOT_FIND=FULL_DIAG) for the Lehmann
C (sum-over-states) screened interaction in scrnc. Header: NMO, IUHF.
C
      OPEN(UNIT=83,FILE='eom_tdm.dat',FORM='FORMATTED',
     &     STATUS='UNKNOWN')
      WRITE(83,*) NOCCO(1)+NVRTO(1),IUHF
C
C READ IN THE NUMBER OF ROOTS REQUESTED IN EACH SYMMETRY BLOCK
C
      CALL GETREC(20,'JOBARC','EESYMINF',NIRREP,NROOT)
      ITROOT=0
      DO 10 IRREPX=1,NIRREP
       ISIDE=1
       IF(NROOT(IRREPX).EQ.0)GOTO 10
       CALL ZERO(R,1000)
       CALL NEWLST(IRREPX,SCR,MAXCOR*IINTFP,IUHF)
       WRITE(6,10000)NSIZEC
       IF(NSIZEC.EQ.0)GOTO 10
       WRITE(6,10001)IRREPX
C
C SET UP ADDRESSES FOR HBAR CONSTRUCTION
C
       I0A=1+NSIZEC
       I0=I0A+NSIZEC
       IF(ESPROP)THEN
        I0L=I0+NSIZEC*NSIZEC
        IA=I0L+NSIZEC*NSIZEC
       ELSE
        IA=I0+NSIZEC*NSIZEC
       ENDIF 
       I000=IA+NSIZEC
       I010=I000+NSIZEC
       I020=I010+NSIZEC
       CALL ZERO(SCR(I0),NSIZEC*NSIZEC)
       IOFFX(1)=I0
       IOFFX(2)=I0L
C
C BUILD HBAR MATRIX
C
        DO 50 I=1,NSIZEC
        CALL ZERO(SCR(IA),NSIZEC)
        SCR(IA-1+I)=ONE
        IF(IUHF.EQ.0)THEN
c 
c for spin adapted basis
c          CALL SPNTSNG2(IRREPX,SCR(IA),SCR(I000),MAXCOR-I000,.TRUE.)
c
c for spin orbital basis
         CALL SYMT2AB(IRREPX,SCR(IA),SCR(I000),MAXCOR-I000)
c 
        ENDIF
C
        CALL UPDATES(IRREPX,SCR(IA),444,0,490,IUHF)
        ITOP=1
        IF(ESPROP)ITOP=2
        DO 55 ISIDE=1,ITOP
         CALL HBARXC(SCR(I000),(MAXCOR-I000+1)*IINTFP,IUHF,ISIDE,
     &               IRREPX)
         CALL LOADVEC1(IRREPX,SCR(IOFFX(ISIDE)),1,IUHF,490,2,460,
     &                 NSIZEC,.FALSE.)
         IF(IUHF.EQ.0)THEN
c
c spin orbital basis
c
          CALL SYMT2AB(IRREPX,SCR(IOFFX(ISIDE)),SCR(I000),MAXCOR-I000)
c
c spin-adapted basis
c
c           CALL SPNTSNG2(IRREPX,SCR(IOFFX(ISIDE)),SCR(I000),
         ENDIF
         IOFFX(ISIDE)=IOFFX(ISIDE)+NSIZEC
55      CONTINUE
50     CONTINUE
C
C HBAR CONSTRUCTION COMPLETE... NOW PROCEED TO DIAGONALIZATION
C
C MATRIX HELD IN SCR(I0) AND SCR(I0L)
C
       NDIMR=NSIZEC
       I000=i000
       I010=I000+NDIMR
       I020=I010+NDIMR
       I030=I020+NDIMR*NDIMR
       I040=I030+NDIMR
       IF(ESPROP)THEN
        I050=I040+NDIMR
        I060=I050+NDIMR
        I070=I060+NDIMR*NDIMR
       ENDIF 
c
CSSS       write(6,*)' denominator '
CSSS       write(6,'((4f20.10))')(scr(j),j=i0,i0+ndimr*ndimr-1,ndimr+1)
       IF (IFLAGS(1) .GE. 40) THEN
        write(6,*)' right hbar matrix '
        call output(scr(i0),1,ndimr,1,ndimr,ndimr,ndimr,1)
       ENDIF
CSSS       Write(6,*)' left hbar matrix '
CSSS      call output(scr(i0l),1,ndimr,1,ndimr,ndimr,ndimr,1)
c
c normalization for spin adapted basis
c
c      if(iuhf.eq.0)call sscal(ndimr*ndimr,0.5d0,scr(i0),1)
cmn       CALL RG(NDIMR,NDIMR,SCR(I0),SCR(I000),SCR(I010),1,
cmn     &          SCR(I020),SCR(I030),SCR(I040),IERR)
c
       CALL MN_GEEV(NDIMR,NDIMR,SCR(I0),SCR(I000),SCR(I010),
     &          SCR(I020),SCR(I030),MAXCOR-I030+1,IERR)
       CALL AFTERRG(NDIMR,SCR(I000),SCR(I010),SCR(I020),SCR(I0),
     &              SCR(I040),11)
       EVALIMX=ZILCH
       DO I=1,NDIMR
        EVALIMX=MAX(EVALIMX,ABS(SCR(I010-1+I)))
       ENDDO
       WRITE(6,10002) EVALIMX
       IF(ESPROP)THEN
C
C Left eigenvectors as the rows of R**-1, not from a second GEEV on the
C left matrix: AFTERRG sorts the two eigenvector sets independently by
C eigenvalue, so a separate left solve pairs L_i with R_i only if both
C sorts agree. R**-1 makes <L_i|R_j> = delta_ij exact by construction.
C The independently built left matrix H_L (HBARXC, ISIDE=2, still
C intact at SCR(I0L)) then checks it: H_L*l_i = w_i*l_i must hold, which
C also verifies H_L = transpose(H_R) in this (orthonormal UHF) basis.
C
        CALL SCOPY(NDIMR*NDIMR,SCR(I020),1,SCR(I060),1)
        IPIV=I070
        IWRK=IPIV+NDIMR
        LWRK=64*NDIMR
        CALL DGETRF(NDIMR,NDIMR,SCR(I060),NDIMR,SCR(IPIV),IERR)
        IF(IERR.NE.0)THEN
         WRITE(6,*)' @DIAGALL-F, DGETRF failed, info = ',IERR
         CALL ERREX
        ENDIF
        CALL DGETRI(NDIMR,SCR(I060),NDIMR,SCR(IPIV),SCR(IWRK),LWRK,
     &              IERR)
        IF(IERR.NE.0)THEN
         WRITE(6,*)' @DIAGALL-F, DGETRI failed, info = ',IERR
         CALL ERREX
        ENDIF
C       rows of R**-1 -> columns of SCR(I060)
        CALL TRANSP(SCR(I060),SCR(I0),NDIMR,NDIMR)
        CALL SCOPY(NDIMR*NDIMR,SCR(I0),1,SCR(I060),1)
        RESMAX=ZILCH
        DO I=1,NDIMR
         IOFFL=I060+(I-1)*NDIMR
         CALL DGEMV('N',NDIMR,NDIMR,ONE,SCR(I0L),NDIMR,SCR(IOFFL),1,
     &              ZILCH,SCR(I040),1)
         CALL SAXPY(NDIMR,-SCR(I000-1+I),SCR(IOFFL),1,SCR(I040),1)
         RESID=DNRM2(NDIMR,SCR(I040),1)/DNRM2(NDIMR,SCR(IOFFL),1)
         RESMAX=MAX(RESMAX,RESID)
        ENDDO
        WRITE(6,10003) RESMAX
       ENDIF
C
C SUMMARIZE
C
C TDENS gets its own scratch past every array still needed below (the
C eigenvalues at I000, right/left eigenvectors at I020/I060): handing
C it SCR(I0) with the full MAXCOR, as before, let its density-build
C scratch overwrite the eigenvectors of roots not yet processed.
C
       ITDS=I070+NDIMR+64*NDIMR
       IF(ITDS.GT.MAXCOR)THEN
        WRITE(6,*)' @DIAGALL-F, insufficient memory for TDENS'
        CALL ERREX
       ENDIF
       IOFFEVCR=I020
       IOFFEVCL=I060
       DO 100 I=1,NDIMR
        EVAL=SCR(I000-1+I)
        IF(ABS(EVAL).GT.1.D-6.AND.ESPROP)THEN
c
c also need reference state contribution if this is the righthand
c solution and irrepx=1
c
         IF(IRREPX.EQ.1)THEN
          CALL LOADVEC1(IRREPX,SCR(I040),MAXCOR,IUHF,93,0,13,NSIZEC,
     &                  IUHF.EQ.0)
          Z0=SDOT(NSIZEC,SCR(I040),1,SCR(IOFFEVCR),1)/EVAL
         ELSE
          Z0=ZILCH
         ENDIF
C
C NORMALIZE RIGHT HAND VECTOR
C
         CALL SCOPY(NSIZEC,SCR(IOFFEVCR),1,SCR(I040),1)
         IF(IUHF.EQ.0)THEN
          CALL SPNTSING(IRREPX,SCR(I040),SCR(I050),1)
         ENDIF   
         X=DSQRT(SDOT(NSIZEC,SCR(IOFFEVCR),1,SCR(I040),1)+Z0*Z0)
         CALL SSCAL(NSIZEC,ONE/X,SCR(IOFFEVCR),1)
         Z0=Z0/X
C
C NORMALIZE LEFT HAND VECTOR SUCH THAT <L|R>=1
C
         CALL SCOPY(NSIZEC,SCR(IOFFEVCL),1,SCR(I040),1)
         IF(IUHF.EQ.0)THEN
          CALL SPNTSING(IRREPX,SCR(I040),SCR(I050),1)
         ENDIF   
         IF(ESPROP)THEN
          X=SDOT(NSIZEC,SCR(IOFFEVCR),1,SCR(I040),1)
          CALL SSCAL(NSIZEC,ONE/X,SCR(IOFFEVCL),1)
         ENDIF
C
         ISIDE=1
         CALL UPDATES(IRREPX,SCR(IOFFEVCR),444,0,490,IUHF)
         CALL PUTLST (SCR(IOFFEVCR),1+ISIDE,1,1,1,472)
         CALL TDENS  (IRREPX,SCR(ITDS),MAXCOR-ITDS+1,IUHF,ISIDE,
     &                (ONE+DFLOAT(1-IUHF)),Z0,EVAL,CORE_SEARCH,
     &                DOUBLE_EES)
         ISIDE=2
         CALL UPDATES(IRREPX,SCR(IOFFEVCL),444,0,490,IUHF)
         CALL PUTLST (SCR(IOFFEVCL),1+ISIDE,1,1,1,472)
         CALL TDENS  (IRREPX,SCR(ITDS),MAXCOR-ITDS+1,IUHF,ISIDE,
     &                (ONE+DFLOAT(1-IUHF)),Z0,EVAL,CORE_SEARCH,
     &                DOUBLE_EES)
C
         CALL PRINTSUM(EVAL,FSTR,.FALSE.,.FALSE.)
c
        ENDIF
        IOFFEVCR=IOFFEVCR+NDIMR
        IOFFEVCL=IOFFEVCL+NDIMR
100    CONTINUE
10    CONTINUE
      CLOSE(UNIT=83)
      RETURN

6000  FORMAT(T3,'@DIAGALL-I, Excitation energies computed by the ',
     &          'EOM-CCSD method.')
10000 FORMAT(T3,'@DIAGALL-I, Physical dimension of H-bar is ',I12,'.')
10001 FORMAT(T3,'@DIAGALL-I, Forming matrix for symmetry block ',I1,'.')
10002 FORMAT(T3,'@DIAGALL-I, Max |Im(eigenvalue)| = ',1PE10.3)
10003 FORMAT(T3,'@DIAGALL-I, Max left-eigenvector residual ',
     &          '|H_L*l - w*l|/|l| = ',1PE10.3)
      end
