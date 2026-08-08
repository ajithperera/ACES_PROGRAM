













































































































































































































      SUBROUTINE QRHFIT(PK,ONEH,DENS,FOCK,EVAL,EVEC,ICORE,MAXMEM,LDIM1,
     &                  LDIM2,NBAS,IPKSIZ,REPULS,ETOT,IUHF,IOS,
     &                  naobasfn,scfks,scfksexact,scfkslastiter,
     &                  V,z1,ksa,ksb,
     &                  screxc,scr,scr2,
     &                  valao,valgradao,totwt,
     &                  max_angpts,natoms,
     &                  intnumradpts,ncount,kshf)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end
c Nevin 6/5/95 modified so QRHF will work with FOCK=AO
      parameter(luint=10)
      logical aofil, closed, noconv, kshf
C
      DIMENSION ONEH(LDIM1),DENS((IUHF+1)*LDIM1),FOCK((IUHF+1)*LDIM1)
      DIMENSION EVEC((IUHF+1)*LDIM2),EVAL((IUHF+1)*NBAS)
      DIMENSION ICORE(1),PK(IPKSIZ)
      DIMENSION IDUMMY(MAXBASFN),IDUMMY2(MAXBASFN),ILOCATE(MAXBASFN)
      DIMENSION IDUMMY3(MAXBASFN), NSUM(8, 2)
      DIMENSION IADD(8,2),IREM(8,2),DOCC(MAXBASFN*2)
c---------------------------------------------------------------------







c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler























cYAU - ACES3 stuff . . . we hope - #include <aces.par>






c This contains flags that are set in the INTGRT namelist.  See the
c file initintgrt.F for a description of each of them.
c
c The following are exceptions:
c    int_ks           : .true. if we are doing Kohn-Sham
c    int_ks_finaliter : .true. if this is the final iteration of KS
c    int_ks_exch      : which potential to use to calculate exchange
c    int_ks_corr      : which potential to use to calculate correlation
c    int_kspot        : which hybrid functional to use to calculate potential
c                       (if equal to fun_special, use int_ks_exch and
c                       int_ks_corr)
c    int_dft_fun      : which functional to use with any SCF density
c                       Added and modified by Stan Ivanov
c                       (if fun_special the functional is user-defined
c                        if fun_hyb_name then use hybrid functional) 
c    int_printlev     : 0 if we are doing a dft calculation, 1 if we
c                       are doing the final iteration of a KS calculation,
c                       2 if we are doing a KS iteration.
c These are set in the calling routines, NOT in the namelist.
c                     : Additions by S. Ivanov
c     num_acc_ks      : .true.  if numerical accelerator is used for KS 
c                        Default is .true.
c     ks_exact_ex     : .true. if exact LOCAL exchange is used for KS
c                        Deafult is .false.
c     int_tdks        : .true. if time-dependent KS calculation is
c                        requested
c                        Default is .false.
c     int_ks_scf      : .true. if the actual KS SCF energy is being
c                        calculated and printed out. Default is .false.
c
      integer int_numradpts,int_radtyp,int_partpoly,int_radscal,
     &    int_parttyp,int_fuzzyiter,int_defenegrid,int_defenetype,
     &    int_defpotgrid,int_defpottype,int_kspot,
     &    int_ks_exch,int_ks_corr,int_dft_fun,

     &    int_printlev,
     &    int_printscf,int_printint,int_printsize,int_printatom,
     &    int_printmos,int_printocc,
     &    potradpts, numauxbas,int_ksmem,int_overlp

      logical int_ks,num_acc_ks,ks_exact_ex,int_tdks,int_ks_scf,
     &        int_ks_finaliter

      double precision
     &    int_radlimit,coef_pot_nonlocal

      common /intgrtflags/  int_numradpts,int_radtyp,int_partpoly,
     &    int_radscal,int_parttyp,int_fuzzyiter,int_defenegrid,
     &    int_defenetype,int_defpotgrid,int_defpottype,int_kspot,
     &    int_ks_exch,int_ks_corr,int_dft_fun,

     &    int_printlev,
     &    int_printscf,int_printint,int_printsize,int_printatom,
     &    int_ks_finaliter,int_printmos,int_printocc,
     &    potradpts, numauxbas,int_ksmem,int_overlp
c

c  prakash added int_ksmem to the common block
      common /intgrtflagsd/ int_radlimit,coef_pot_nonlocal
      common /intgrtflagsl/ int_ks,num_acc_ks,ks_exact_ex,int_tdks,
     &                      int_ks_scf

      save /intgrtflags/
      save /intgrtflagsl/
      save /intgrtflagsd/

c The following are parameters used in the namelist

      integer int_prt_never,int_prt_dft,int_prt_ks,int_prt_always
      parameter (int_prt_never   =1)
      parameter (int_prt_dft     =2)
      parameter (int_prt_ks      =3)
      parameter (int_prt_always  =4)

      integer int_radtyp_handy,int_radtyp_gl
      parameter (int_radtyp_handy=1)
      parameter (int_radtyp_gl   =2)

      integer int_partpoly_equal,int_partpoly_bsrad,
     &    int_partpoly_dynamic
      parameter (int_partpoly_equal  =1)
      parameter (int_partpoly_bsrad  =2)
      parameter (int_partpoly_dynamic=3)

      integer int_radscal_none,int_radscal_slater
      parameter (int_radscal_none  =1)
      parameter (int_radscal_slater=2)

      integer int_parttyp_rigid,int_parttyp_fuzzy
      parameter (int_parttyp_rigid=1)
      parameter (int_parttyp_fuzzy=2)

      integer int_gridtype_leb
      parameter (int_gridtype_leb=1)


c
       integer iuhf,naobasfn,ncount,natoms,max_angpts,
     &         intnumradpts
       logical scfks,scfksexact,scfkslastiter
       integer z1(naobasfn,2)
       double precision V(naobasfn,naobasfn,iuhf+1),
     &                  ksa(ldim2),ksb(ldim2),
     &                  screxc(naobasfn,naobasfn),
     &                  scr(naobasfn,naobasfn,iuhf+1),
     &                  scr2(nbas,nbas),coef_nonloc,
     &                  valao(naobasfn,intnumradpts,
     &                        max_angpts,ncount),
     &                  valgradao(naobasfn,intnumradpts,
     &                            max_angpts,ncount,3),
     &                  totwt(ncount,intnumradpts,max_angpts)     
c----------------------------------------------------------------------------
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /POPUL/ NOCC(8,2)
      LOGICAL BRUEK
c Nevin 6/5/95 modified so QRHF will work with FOCK=AO
      common /fock/ aofil
c symm2.com : begin

c This is initialized in vscf/symsiz.

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
      DATA ONE /1.0/
      DATA TWO /2.0/
C
      INDX2(I,J,N)=I+(J-1)*N
C
      IONE=1
      CALL IZERO(IREM,16)
      CALL IZERO(IADD,16)
C
C SET UP ARRAY WHICH MATCHES EIGENVALUE POSITION WITH ORIGINAL
C POSITION IN ARRAY.
C
      DO 1000 I=1,NBAS
       ILOCATE(I)=I
1000  CONTINUE
C
C  Make the adjustment in the occupation vector based on the QRHF
C  flags.
C
      CALL GETREC(20,'JOBARC','QRHFTOT ',IONE,NMODIFY)
cYAU - This is probably unnecessary since one cannot drop more orbitals
c      than actually exist.
      if (nmodify.gt.maxbasfn) then
         print *, '@QRHFIT: Assertion failed.'
         print *, '         nmodify = ',nmodify
         print *, '         allowed = ',maxbasfn
         call errex
      end if
C
C If Brueckner orbitals are requested with QRHF we need to save
C the occupation vector that correspond to the pre-qrhf (see also
C initges.F in vscf). A. Perera, 10/2013.
C
      BRUEK = .FALSE.
      BRUEK = (IFLAGS(22) .GT. 0)
      IF (BRUEK) THEN
          CALL PUTREC(20, "JOBARC", "ORGOCCA", 8, NOCC(1,1))
          CALL PUTREC(20, "JOBARC", "ORGOCCB", 8, NOCC(1,2))
      ENDIF
C
      CALL GETREC(20,'JOBARC','QRHFIRR ',NMODIFY,IDUMMY)
      CALL GETREC(-1,'JOBARC','QRHFLOC ',NMODIFY,IDUMMY2)
      CALL GETREC(-1,'JOBARC','QRHFSPN ',NMODIFY,IDUMMY3)
      Write(6,"(a)") " QRHF parameters"
      write(6,"(6(1x,4i1))")(Idummy(j), j=1,nmodify)
      write(6,"(6(1x,4i1))")(Idummy2(j),j=1,nmodify)
      write(6,"(6(1x,4i1))")(Idummy3(j),j=1,nmodify)
      IF(IFLAGS(77).NE.0) THEN
       DO 401 I=1,NMODIFY 
        IRRP=IDUMMY(I)
        ILOC=IDUMMY2(I)
        ISPN=IDUMMY3(I)

        IF(IRRP.LT.0)THEN
         IF(ISPN.NE.1) ISPN=2
         IRRP=-IRRP
         IPOSORG=IRPOFF(IRRP)+NOCC(IRRP,ISPN)+1-MAX(ILOC,1)+
     &                IREM(IRRP,ISPN)-IADD(IRRP,ISPN)
         IPOSABS=ILOCATE(IPOSORG)
         IREM(IRRP,ISPN)=IREM(IRRP,ISPN)+1
         CALL QRHFPOP(IRRP,ILOC,NBAS,EVAL,EVEC,IPOSABS,LDIM2,
     &                ILOCATE,ICORE,ISPN,IOS)
        ELSE
         IF(ISPN.NE.2) ISPN=1
         IPOSORG=IRPOFF(IRRP)+NOCC(IRRP,ISPN)+MAX(ILOC,1)-
     &                IADD(IRRP,ISPN)+IREM(IRRP,ISPN)
         IPOSABS=ILOCATE(IPOSORG)
         IADD(IRRP,ISPN)=IADD(IRRP,ISPN)+1
         CALL QRHFADD(IRRP,ILOC,NBAS,EVAL,EVEC,IPOSABS,LDIM2,
     &                ILOCATE,ICORE,ISPN,IOS)
        ENDIF
401    CONTINUE
C
      ENDIF
C
C  Write out eigenvalues
C
      CALL PUTREC(20,'JOBARC','QRHFEVAL',NBAS*IINTFP,EVAL)
C
C  Now construct a new density matrix.
C
      I000=1
      I010=I000+MXIRR2*IINTFP
      I015=I010+MXIRR2*IINTFP
      IF(I015-I000.GT.MAXMEM) THEN
        CALL NOMEM('Construct density','{MKDEN} <-- QRHFIT <-- VSCF',
     &             I010-I000,MAXMEM)
      ENDIF
C
      DO 100 ISPIN=1,(IUHF+1)
        IROFF = (ISPIN - 1)*NBAS
        DO 101 I=1,NIRREP
          IF(NBFIRR(I).EQ.0) GOTO 101
          CALL MKDEN(EVEC((ISPIN-1)*LDIM2+ISQROF(I)),ICORE(I000),
     &               NOCC(I,ISPIN),NBFIRR(I),IUHF)
          CALL SQUEZ2(ICORE(I000),DENS((ISPIN-1)*LDIM1+ITRIOF(I)),
     &                NBFIRR(I))
  101   CONTINUE
  100 CONTINUE

C
C  Now construct the new Fock matrix.
C
c Nevin 6/5/95 modified so QRHF will work with FOCK=AO
      if(aofil) then
        ilnbuf=600
        i010=i000+itriln(nirrep+1)*iintfp
        i020=i010+iintfp*ilnbuf
        i030=i020+ilnbuf
        i040=i030+nbas*nbas
        IF(I040-I000.GT.MAXMEM) THEN
          CALL NOMEM('Make Fock matrix','{MKUHFF} <-- QRHFIT <-- VSCF',
     &             I040-I000,MAXMEM)
        ENDIF
        call mkuhff(fock,fock(1+itriln(nirrep+1)),dens,
     &              dens(1+itriln(nirrep+1)),icore(i000),oneh,
     &              icore(i010),icore(i020),itriln(nirrep+1),
     &              nbas,nbfirr,icore(i030),ilnbuf,luint,.true.,
     &              naobasfn,iuhf,scfks,scfksexact,
     &              scfkslastiter,
     &              V,z1,ksa,ksb,
     &              screxc,scr,scr2,
     &              valao,valgradao,totwt,
     &              max_angpts,natoms,
     &              intnumradpts,ncount,kshf)

      ELSE
        I010=I000+(IUHF+1)*ITRILN(NIRREP+1)*IINTFP
        I020=I010+MXIRR2*IINTFP
        i030=i020+nbas*nbas*iintfp
        IF(I030-I000.GT.MAXMEM) THEN
          CALL NOMEM('Make Fock matrix','{MKFOCK} <-- QRHFIT <-- VSCF',
     &             I030-I000,MAXMEM)
        ENDIF
C
        CALL MKFOCK(PK,ONEH,FOCK,DENS,ICORE(I000),ICORE(I010),
     &            ICORE(I020),IPKSIZ,ITRILN(NIRREP+1),
     &            MXIRR2,nbas,IUHF)
      ENDIF
c
C  Now determine the QRHF energy.  The value for TOL is irrelevant, since
C  it is used to determine whether to use a different format statement for
C  the energy.  This section in MKENER is skipped for a QRHF case.
C
      TOL=0.0
      DMAX=0.0
      ITER=0
      noconv=.false.
      I010=I000+ITRILN(NIRREP+1)*IINTFP
      I020=I010+MXIRR2*IINTFP
      I030=I020+MXIRR2*IINTFP
      IF(I030-I000.GT.MAXMEM) THEN
        CALL NOMEM('QRHF energy','{MKENER} <-- QRHFIT <-- VSCF',
     &             I030-I000,MAXMEM)
      ENDIF
C
      CALL MKENER(ONEH,DENS,FOCK,ICORE(I000),ICORE(I010),ICORE(I020),
     &            ITRILN(NIRREP+1),MXIRR2,REPULS,DMAX,ITER,IUHF,ETOT,
     &            1,NOCONV,scfksiter,scfkslastiter)
C
      WRITE(LUOUT,5000)ETOT
 5000 FORMAT(T3,'The QRHF reference energy is ',F20.10,/)
C
C  At this point, the values needed for the QRHF calculation have
C  been placed in the appropriate arrays and will be dumped out
C  with everything else in DMPJOB.
C
C  We also need to calculate the new S**2 value, based on the new
C  occupancies.  This reference function is an eigenfunction of
C  spin, so we only determine the uncontaminated value.
C
      NALPHA=0
      NBETA=0
      DO 200 I=1,NIRREP
        NALPHA=NALPHA+NOCC(I,1)
        NBETA=NBETA+NOCC(I,2)
  200 CONTINUE
      SAVG=ABS(NALPHA-NBETA)/TWO
      S2=SAVG*(SAVG+ONE)
      AMULT=SQRT(1.0+4.0*S2)
      WRITE(LUOUT,5100)AMULT,S2
 5100 FORMAT(/,T3,'     The QRHF average multiplicity is ',F12.7,/,
     &         T3,'The QRHF expectation value of S**2 is ',F12.7,/)
      CALL PUTREC(20,'JOBARC','S2SCF   ',IINTFP,S2)
C
Cmn print some further information
c
      write(6,*)
      write(6,*) ' Final QRHF population'
      WRITE(LUOUT,5302)(NOCC(I, 1),I=1,NIRREP)
      WRITE(LUOUT,5303)(NOCC(I, 2),I=1,NIRREP)
 5302 FORMAT(T8,'   Alpha population by irrep: ',8(I3,2X))
 5303 FORMAT(T8,'    Beta population by irrep: ',8(I3,2X))
      write(6,*)
C
C CHANGE TO CLOSED SHELL CALCULATION IF THE STATE IS REOCCUPIED TO
C CLOSED SHELL, AND SCF-REF WAS RHF OR ROHF  AND MAKERHF=ON
C I Think  better logic is to check whether it is actually a
C closed shell and MAKERHF=ON, 03/2011, Ajith Perera.

CSSS      IF ((IFLAGS(11) .EQ. 0 .OR. IFLAGS(11) .EQ. 2)
CSSC     $     .AND. IFLAGS2(140) .eq. 1) THEN
         Nalpha = 0
         Nbeta  = 0
         CLOSED = .FALSE.
         DO IRREP = 1, NIRREP
            Nalpha = Nalpha + NOCC(IRREP,1)
            Nbeta  = Nbeta  + NOCC(IRREP,2)
         ENDDO
         If (Nalpha .EQ. Nbeta) CLOSED=.TRUE.

         IF (CLOSED .AND. IFLAGS2(140) .eq. 1) THEN
            WRITE(6,*) ' The QRHF reference state is closed shell'
            WRITE(6,*)
            IUHF = 0
            iflags(11) = iuhf
            nflags=100
            call putrec(20, 'JOBARC','IFLAGS  ',nflags,iflags)
         endif
CSSS      endif
C
      RETURN
      END
