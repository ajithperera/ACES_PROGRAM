











C THE PRIMARY PURPOSE OF THE INTProc MODULE IS TO PRODUCE THE
C  ORDERED MO INTEGRAL LISTS WHICH ARE EVENTUALLY WRITTEN TO
C  THE MOINTS FILE.  ALSO, A NUMBER OF TWO INDEX QUANTITIES
C  SUCH AS FOCK MATRIX LISTS ARE FORMED AS WELL.
C
C CODED BY J.F. STANTON AND J. GAUSS
cjp MR-BW-CC extension added by Jiri Pittner (1998-2000)





































































































































































































      SUBROUTINE INTPROC(ICORE,MAXCOR,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER ICORE(MAXCOR),MAXCOR,IUHF,I0
      INTEGER DIRPRD,USEC,ABCDTYPE
      LOGICAL DOPHPH,ABIJ,DOALL,NOABCD
      LOGICAL COMPRESS,SYMMETRIC
      LOGICAL MACRO_ITER
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /INCORE/ ICREC(2),USEC,IXT(2),IMOD(2)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /ABCD/ ABCDTYPE
      COMMON /DOINTS/ DOALL,ABIJ,NOABCD
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMLOC/ ISYMOFF(8,8,25)
C
CJDW  KKB stuff.
C
      COMMON /SYMPOP2/ IRP_PD(8,22)
      COMMON /SHIFT/   ISHIFT,NDRGEO
CJDW END
cjp
cjp switches to be later implemented as aces namelist options
cjp


cjp
cjp data for multireference state specific Brillouin-Wigner CC method
cjp coded by Jiri Pittner 1998-2000
cjp
      logical isbwcc,masik,isactive,bwgossip,useeq429,scfrefread
      logical bwwarning
      character*256 bwwarntext
      real*8 ecorrbw0,ecorrbw,epsilon0,fockcontr,denomblow,fockcd
      real*8 heff,heffevalr,heffevali,heffevecl,heffevecr,hdiagcontr
      real*8 fock2elcontr,enerscf,hcore,lambdahomotop,hfakt,diishonset
      real*8 fock2elcontr0,enerscf0,cbwstate,totmaxdenom,heffevecrold
      real*8 intruder,hfaktmax
      integer maxorb,maxref,nref, iref, iocc,iocc0
      integer iphnum,invpnum,invhnum,nbwstates,ibwstate
      integer ibwconvg,internfrom,internto,internnum,internindex
      integer ibwpass
      integer nactive,numactive,ihubaccorr,ihomotop
      integer internfrom1,internto1,internindex1,internnum1
      integer ihefferank, iheffefrom, iheffeto, iheffespin, maxexcit
      integer correctiontype
      integer maxbwwarnings, nproc, myproc
c
      parameter(maxorb=512,maxref=32,maxexcit=9,maxbwwarnings=10)
c      NOTE!!! change of maxorb parameter requires format change
c              and character* change in bwread routine!!!
      parameter(denomblow=1d250)
c

cjp common has been splitted in order to avoid problems
cjp with padding on different 32 and 64 bit architectures

      common/bwccint/isbwcc, masik, nref, iref,iocc(maxorb,maxref,2),
     +     iphnum(maxorb,maxref,2),invpnum(maxorb,maxref,2),
     +     invhnum(maxorb,maxref,2),
     +     isactive(maxorb,2),nbwstates,ibwstate(maxref+1),
     +     internfrom(maxref*(maxref-1)/2,maxref,3),
     +     internto(maxref*(maxref-1)/2,maxref,3),
     +     internindex(maxref*(maxref-1)/2,maxref,3),
     +     internnum(maxref,3),
     +     internfrom1(maxref*(maxref-1)/2,maxref,2),
     +     internto1(maxref*(maxref-1)/2,maxref,2),
     +     internindex1(maxref*(maxref-1)/2,maxref,2),
     +     internnum1(maxref,2),
     +     ibwpass,ibwconvg(maxref),bwgossip,useeq429,
     +     nactive(2),numactive(maxorb,2),ihubaccorr, ihomotop,
     +     iocc0(maxorb,2),scfrefread,
     +     ihefferank(maxref,maxref),iheffefrom(maxexcit,maxref,maxref),
     +     iheffeto(maxexcit,maxref,maxref),
     +     iheffespin(maxexcit,maxref,maxref),
     +     correctiontype,bwwarning(maxbwwarnings),
     +     bwwarntext(maxbwwarnings),nproc,myproc

      common/bwccreal/ecorrbw,epsilon0,cbwstate(maxref+1),
     +     fockcontr(maxorb*(maxorb+1)/2,2),fockcd(maxorb,maxref,2),
     +     heff(maxref,maxref),heffevalr(maxref),heffevali(maxref),
     +     heffevecl(maxref,maxref),heffevecr(maxref,maxref),
     +     hdiagcontr(maxref), fock2elcontr(maxorb,2),
     +     enerscf(maxref), hcore(maxorb,2),
     +     lambdahomotop,hfakt,diishonset,enerscf0,
     +     fock2elcontr0(maxorb,2),ecorrbw0,totmaxdenom,
     +     heffevecrold(maxref,maxref),intruder,hfaktmax

c
c
cjp BRIEF DESCRIPTION OF VARIABLES INTRODUCED FOR THE MR-BWCC ROUTINES
cjp IN FACT, A LOT OF THAT COULD BE USEFUL FOR ANY HILBERT-SPACE MR-CC
c
c
c isbwcc ... flag for doing bwcc calculation
c maxbwwarnings, bwwarning, bwwarntext ... serious warnings will be
c    summarized at the end of xvcc output for the user's' convenience
c ihefferank(jref,iref) ... degree of excitation between jref and iref
c iheffefrom(maxexcit,jref,iref) , iheffeto, iheffespin ... list of
c    indices of that excitation, sorted according to spin and then the
c    indices, numbers stored are defined as effective particle-hole
c    indices of reference iref
c ihomotop ... whether to use homotopic transition to the
c    size-extensivity correction, after which iteration (if .ne.0)
c lambdahomotop scaling factor of the geometrical series of
c    lambda 1->0 transition
c hfakt ... current value of the homotopy parameter
c hfaktmax ... maximal homotopy parameter allowed to consider cc
c    equations converged
c diishonset ... at which value of hfact restart diis convergence acceleration
c masik ... prepare sorted integral file for the program by Masik and stop
c nref ... number of reference configurations
c iref ... current reference configuration and fermi vacuum
c bwgossip ... switch on debugging output
c ibwpass ... routines like newt2 have to be splitted in two passes -
c    construction of Heff and amplitude update after heff is diagonalized
c    for backw. compatibility, instead of introducing a new routine
c    the same routine does different things being called twice with
c    different ibwpass value
c ecorrbw ... correlation energy from BWCC - Heff(iref,iref) ...
c    denominator correction
c ecorrbw0 .... ecorrbw, but not scaled by the homotopic factor hfact
c denomblow ... huge number to cause division underflow - used for
c    zeroing out the internal amplitudes automatically
c nactive(spin): total count of active spinorbitals
c numactive(i=1..nactive,spin): number of i-th active spinorbital
c    in sequential numbering
c isactive(maxorb,spin): belongs given orbital to the active space?
c    for RHF, the beta ones must be initialized to be identical with alpha ones
c iocc(maxorb,1..nref,spin): defines the nref reference configurations
c    for both RHF and UHF:  iocc(i,iref,spin)=0 or 1
c iphnum(orbital no, iref, spin): gives the effective number of orbital
c    (both particle and hole ones are counted starting from 1)
c invpnum(eff.p.orb.,iref,spin): gives true orbital no. from the
c    effective particle one
c invhnum(eff.h.orb.,iref,spin): gives true orbital no. from the
c    effective hole one
c    all these three ones must be in RHF case initialised to be equal in the
c    alpha and beta parts to keep the code unique and simple
c internfrom(sequence counter n,iref,ispin) is for the ab spin case, the
c    other ones have to be iuhf-indexed
c    internfrom, internto: they are first and second index of n-th internal
c    excitation when processing reference given by second index to the array
c
c internindex(sequence counter,iref,ispin) is the position of
c    corresponding denominator in the denominator list
c internnum(iref,ispin) - number of internal excitation in that
c    category = max sequence counter here ispin=1,2,3 for AA,BB,AB
c internfrom1 etc. are analogous quantities for monoexcitations, here
c    ispin =1,2 note for later: all intern.... quantities are irrep-specific!
c fockcontr(findex(i,j),ispin) ... addition to the fock matrix of
c    reference no.1 to obtain the fock matrix of current reference
c    (fermi vacuum)
c fockcd(i,iref,ispin) ... diagonal part of that correction for ref. no. iref.
c hcore(i,ispin) ... one electron diagonal hamiltonian elements
c fock2elcontr(i,ispin) ... 2el contribution to the diagonal fock element
c    used temporarily
c hdiagcontr(iref) ... contribution of differences of HF energies of
c    different Fermi vacua to diagonal Heff elements
c enerscf(iref) ... HF energy of iref-th fermi vacuum
c ihubaccorr ... =1 ... calculate the size extenzivity correction for BWCC
c                =2,3 ... second and third pass of that calculation
c iocc0(maxorb,2) ... like iocc, but for dummy reference configuration
c    corresponding to SCF WF
c fock2elcontr0(maxorb,2) ... like fock2elcontr0 but for dummy reference
c    configuration
c enerscf0 ... like enerscf, but for dummy reference config
c scfrefread ... tells to bwprep that SCF reference has been read from input
c    and should not be generated automatically from nocc(ispin)
c nbwstates ... how many states to average
c ibwstates(1..nbwstates),cbwstates() ... their numbers and coefficients
c correctiontype ... 0=DC,L T2 term is removed/scaled, 1=DC/L term is not
c    removed/scaled
c totmaxdenom ... max 1/denom found for given reference's' fermi vacuum -
c    as indication of possible intruder problem
c intruder ... limit of 1/denom to be considered intruder and its
c    amplitude zeroed
c for parallelization
c nproc ... number of processors (counted from 1)
c myproc ... number of the processor currently executing the code
c    (counted from 1)


cjp
C Parametrized CC parameters. 

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 

c
      call bwread(.false.,iuhf)
      call parread(iuhf)
C
C INITIALIZE MACHINE DEPENDENT PARAMETERS
C
c   o the pyaces conversion left this subroutine's own top-level
c     `CALL CRAPSI(ICORE(1),IUHF,0)` in place even though ICORE/MAXCOR/
c     IUHF are now caller-supplied dummy args, not this routine's own
c     allocation -- calling it again here corrupted the caller-supplied
c     MAXCOR (observed live: MAXCOR jumped from the correct caller value
c     to IFLAGS(36), CRAPSI's own internally-recomputed request size,
c     immediately after this call returned), causing an out-of-bounds
c     write a few lines below and a SIGSEGV in every standalone run.
c     The caller (see the PROGRAM wrapper at the bottom of this file)
c     already does the equivalent CRAPSI/aces_init call once per
c     process, which is all that's needed -- removed here.
      if (isbwcc) iref=1
      if (isbwcc) call bwprep(nocco,nvrto,iuhf)
C
      CALL SETMET_PROC
      I0     = 1
      ITWO   = 2
cjp prepare for reopenmo
      do iref=1,nref
         call storemoio
      end do


      do 1234 iref = 1, nref
      if (iref.gt.nproc) call reopenmo

      totmaxdenom=0.

c----------------------------------------------------------------
c----    To handle the drop-mo in energy gradient  --  KB -------
c----------------------------------------------------------------
      call getrec(20,'JOBARC','NDROPGEO',1,NDRGEO)
      if (ndrgeo.eq.2) then
         ISHIFT = 300
      else
         ISHIFT = 0
      end if
c----------------------------------------------------------------

      NMO = NOCCO(1) + NVRTO(1)
cYAU      I0IRRBET=I0+MAXCOR-NMO
cYAU      I0IRRALP=I0IRRBET-NMO
cYAU      CALL GETREC(20,'JOBARC','IRREPALP',NMO,ICORE(I0IRRALP))
cYAU      CALL GETREC(20,'JOBARC','IRREPBET',NMO,ICORE(I0IRRBET))
cYAU      IPASS=I0IRRALP
cYAU      MAXCOR = MAXCOR + ( -2 * NMO )
      MAXCOR = MAXCOR - NMO - NMO
      IPASS = I0 + MAXCOR
      CALL GETREC(20,'JOBARC','IRREPALP',NMO,ICORE(IPASS      ))
      CALL GETREC(20,'JOBARC','IRREPBET',NMO,ICORE(IPASS + NMO))

      ABIJ   = .FALSE.
      NOABCD = .FALSE.
      DOALL  = .TRUE.
      IF (IFLAGS(83).EQ.0) THEN
         IF (IFLAGS(2).LE.1) DOALL=.FALSE.
         IF (IFLAGS(2).EQ.0) NOABCD=.TRUE.
         IF (IFLAGS(2).EQ.1.AND.
     &       IFLAGS(3).LE.1) NOABCD=.TRUE.
         IF (IFLAGS(2).EQ.1.AND.
     &       IFLAGS(18).GE.3.AND.IFLAGS(18).LT.9) THEN
            NOABCD=.TRUE.
         END IF
CMN
         IF (IFLAGS(2).EQ.1.AND.
     &       IFLAGS(3).EQ.0.AND.IFLAGS(74).EQ.0.AND.IFLAGS(22).EQ.0.AND.
     &       IFLAGS(87).EQ.0) THEN
            ABIJ=.TRUE.
         END IF
CMN END
cSG
c 4/25/96 For P-EOM-MBPT(2) we only need <AB|CD> for gradients
         IF (IFLAGS(87).GT.0) THEN
            IF ((IFLAGS(87).EQ.7.OR.IFLAGS(87).EQ.8).AND.
     &          IFLAGS(3).EQ.0.AND.IFLAGS(2).EQ.1) THEN
               NOABCD = .TRUE.
               DOALL  = .FALSE.
            ELSE
               NOABCD = .FALSE.
               DOALL  = .TRUE.
            END IF
         END IF
C
      ELSE IF (IFLAGS(83).EQ.1) THEN
         DOALL  = .TRUE.
         NOABCD = .FALSE.
      ELSE IF(IFLAGS(83).EQ.2) THEN
         DOALL  = .FALSE.
         NOABCD = .TRUE.
      END IF
      IF (IFLAGS(93).EQ.2) NOABCD=.TRUE.
CJDW
c---------------------------------------------------------------------
c---  To handle the drop-mo in energy gradient    --------  KB -------
c---  For Drop-Mo cases : without drop any virtual mo   --------------
c---  MOABCD is not recalculated in the second run of intprc   -------
c---------------------------------------------------------------------
      CALL GETREC(20,'JOBARC','NDROTVRT',1,NDROTVRT)
c     CALL GETREC(20,'JOBARC','IFLDMOGD',1,IFLDMGD)
      IF (.NOT.NOABCD .AND. NDRGEO.EQ.2) THEN
         IF (NDROTVRT.EQ.0) THEN
c        IF (NDROTVRT.EQ.0 .AND. IFLDMGD.NE.0 ) THEN
c           IFLDMGD = 2
            NOABCD  = .TRUE.
c           CALL PUTREC(20,'JOBARC','IFLDMOGD',1,IFLDMGD)
         END IF
      END IF
c---------------------------------------------------------------------
CJDW END
C
      CALL GETREC(20,'JOBARC','NOCCORB ',ITWO,NOCCO)
      CALL GETREC(20,'JOBARC','NVRTORB ',ITWO,NVRTO)
      IF (IUHF.EQ.0) THEN
         NOCC = NOCCO(1)
         NVRT = NVRTO(1)
      END IF
C
C  DO INITIAL SORTING OF INTEGRALS.  FORM PPPP, PPPH, ETC. FILES.
C
      call aces_io_reset
      call aces_cache_reset

      NSTO=NOCCO(1)+NVRTO(1)
C
CJDW  KKB stuff.
C     CALL DGMOI(ICORE(I0),MAXCOR,NSTO,NIRREP,IUHF)
C
      IF (NDRGEO.EQ.1) THEN
         CALL DGMOID(ICORE(I0),MAXCOR,NSTO,NIRREP,IUHF)
      ELSE
         CALL DGMOI(ICORE(I0),MAXCOR,NSTO,NIRREP,IUHF)
      END IF
CJDW END
C
      CALL CLMOIO(ICORE(I0),MAXCOR,NOCCO,NVRTO,IUHF)
C
C PROCESS HHHH INTEGRALS.  THIS IS DONE IN CORE.
C
      IF (.NOT.ABIJ) THEN
         IF (IUHF.NE.0) THEN
            CALL DS16AA(ICORE(I0),MAXCOR,1,IUHF,1,
     &                  ICANT,NMO,ICORE(IPASS))
            CALL DS16AA(ICORE(I0),MAXCOR,1,IUHF,2,
     &                  ICANT,NMO,ICORE(IPASS))
         END IF
         CALL DS16AB(ICORE(I0),MAXCOR,1,IUHF,ICANT,NMO,ICORE(IPASS))
C
C PROCESS PHHH INTEGRALS.  ALSO DONE IN CORE.  DONE FOR
C  E(2) GRADIENT AND HIGHER ENERGY CALCULATIONS.
C
            CALL DS25AA(ICORE(I0),MAXCOR,2,IUHF,1,IOUT,NMO,ICORE(IPASS))
         IF (IUHF.EQ.1) THEN
            CALL DS25AA(ICORE(I0),MAXCOR,2,IUHF,2,IOUT,NMO,ICORE(IPASS))
            CALL DS25AB(ICORE(I0),MAXCOR,2,IUHF,       NMO,ICORE(IPASS))
         ELSE
         IF (IOUT.EQ.0) THEN
            CALL DS25AB(ICORE(I0),MAXCOR,2,IUHF,       NMO,ICORE(IPASS))
         END IF
         END IF
      END IF
C
C PROCESS PPHH AND PHPH INTEGRALS.  THIS IS DONE IN CORE.
C
      DOPHPH = .NOT.ABIJ
      DO I = 1+IUHF, 1, -1
         CALL DS3AA(ICORE(I0),MAXCOR,3,IUHF,I,NMO,ICORE(IPASS))
         IF (DOPHPH) THEN
            CALL DS4AA(ICORE(I0),MAXCOR,4,IUHF,I,NMO,ICORE(IPASS))
         END IF
      END DO
      IF (IUHF.NE.0) THEN
         CALL DS3AB(ICORE(I0),MAXCOR,3,IUHF,3,NMO,ICORE(IPASS))
         IF (DOPHPH) THEN
            CALL DS4AB(ICORE(I0),MAXCOR,4,IUHF,3,NMO,ICORE(IPASS))
         END IF
      END IF
C
C Giving the users to choose whether they want no REDUNDANT lists
C option. This may save some  (i doubt that it is siginificant)
C disk space with a modest (i hope) increase in CPU time. This
C option was included to address some of Anthony's concerns about
C REDUNDANT storage despite these are not the disk intensive
C terms. Ajith Perera 07/2002

      IF (iFlags2(155).eq.0) THEN
         CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'W')
      END IF
C
C  DO PPPH INTEGRALS.
C
      IF (.NOT.ABIJ) THEN
         IF (IUHF.NE.0) THEN
            CALL DS25AA(ICORE(I0),MAXCOR,5,IUHF,1,IOUT,NMO,ICORE(IPASS))
            CALL DS25AA(ICORE(I0),MAXCOR,5,IUHF,2,IOUT,NMO,ICORE(IPASS))
         END IF
         CALL DS25AB(ICORE(I0),MAXCOR,5,IUHF,NMO,ICORE(IPASS))
      END IF
C
C NOW DO THE PPPP INTEGRALS.
C
      IF (.NOT.NOABCD) THEN
c         IMODE3 = 1
         call aces_io_remove(52,'MOABCD')
         IMODE3 = 0
         IF (IUHF.NE.0) THEN
            CALL INIPCK(1,1,1,231,IMODE3,0,1)
            CALL INIPCK(1,2,2,232,IMODE3,0,1)
            CALL INIPCK(1,13,13,233,IMODE3,0,1)
         ELSE
C
CMN  DETERMINE IF ABCD INTEGRALS ARE TO BE COMPRESSED OR NOT
C
            COMPRESS = IFLAGS2(107) .EQ. 2
            IF (.NOT. COMPRESS) THEN
               CALL INIPCK(1,13,13,233,IMODE3,0,1)
            ELSE
               CALL INIPCK(1,5,13,233,IMODE3,0,1)
            END IF
         END IF
         IF (IUHF.NE.0) THEN
            CALL DS16AA(ICORE(I0),MAXCOR,6,IUHF,1,
     &                  ICANT,NMO,ICORE(IPASS))
            CALL DS16AA(ICORE(I0),MAXCOR,6,IUHF,2,
     &                  ICANT,NMO,ICORE(IPASS))
         END IF
         CALL DS16AB(ICORE(I0),MAXCOR,6,IUHF,ICANT,NMO,ICORE(IPASS))
      END IF
C
C PROCESS FOCK MATRIX
C
      CALL DFOCK(ICORE(I0),MAXCOR,IUHF)
CJDW
c-------------------------------------------------------------------
c-----    end of w2 processing   --------- May 94,  KB  ------------
c-------------------------------------------------------------------
      if (ishift.ne.0) then
         CALL PUTSTF(ICORE(I0),MAXCOR,IUHF)
         RETURN
      end if
c-------------------------------------------------------------------
CJDW END
C
C WRITE DENOMINATOR ARRAYS AND FORM INITIAL T2 VECTOR.
C
      SYMMETRIC = (IFLAGS(2) .EQ. 48 .OR.
     +             IFLAGS(2) .EQ. 49)

      CALL DDIJAB(ICORE(I0),MAXCOR,IUHF)
       CALL DFRMT2(ICORE(I0),MAXCOR,IUHF)
      IF (IUHF.NE.0) THEN
         CALL DMOSMT(ICORE(I0),MAXCOR,IUHF)
         CALL S2PROJ(ICORE(I0),MAXCOR,IUHF,.FALSE.)
      END IF
C
C PRINT OUT DOMINANT AMPLITUDES FROM FIRST-ORDER WAVEFUNCTION.
C
c     CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,.FALSE.,'T')
      IF (.NOT. SYMMETRIC) THEN
         CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,isbwcc,'T')
      ELSE
         WRITE(6,*)
         WRITE(6,"(16x,a)") " ----------WARNING--------------"
         WRITE(6,"(a,a)") " For rCCD or drCCD the second-order",
     &                    " energy printed here is incorrect!" 
      ENDIF 
      CALL PUTSTF(ICORE(I0),MAXCOR,IUHF)
C
C PERFORM HARTREE-FOCK STABILITY ANALYSIS
C
      IF (IFLAGS(74).NE.0) CALL STABLE(ICORE(I0),MAXCOR,IUHF)
      IF (IFLAGS(1).GE.5) CALL aces_io_summary
      IF (IFLAGS(1).GE.1) THEN
         CALL MEMUSE
         CALL FLOPUSE(IUHF)
      END IF

      if (isbwcc) then
         write(*,*) '@INTPRC: iref=',iref,
     &              ' TOTAL MAX 1/denomin: ',totmaxdenom
      end if

      if (nref.gt.1) then
         call aces_cache_flush
         call aces_cache_reset
         call storemoio
      end if

CSSS      call check_ints(icore(i0),Maxcor,iuhf,.false.)

 1234 continue

      RETURN
      END


C ----------------------------------------------------------------------
C Thin PROGRAM wrapper restoring the standalone xINTPRC entry point
C (MAIN__), lost when this module was converted to a SUBROUTINE for
C pyaces (see project memory: ACES_PROGRAM unification, 2026-08). This
C calls the shared aces_init/aces_fin pair so the standalone binary
C manages its own memory exactly as the classic driver's own separate
C process used to, before the pyaces conversion removed that from
C inside the subroutine itself.
C ----------------------------------------------------------------------
      PROGRAM XINTPRC
      IMPLICIT NONE


c icore.com : begin

c icore(1) is an anchor in memory that allows subroutines to address memory
c allocated with malloc. This system will fail if the main memory is segmented
c or parallel processes are not careful in how they allocate memory.

      integer icore(1)
      common / / icore

c icore.com : end





c istart.com : begin
      integer         i0, icrsiz
      common /istart/ i0, icrsiz
      save   /istart/
c istart.com : end
      INTEGER IUHF
      CALL CRAPSI(ICORE, IUHF, 0)
      CALL INTPROC(ICORE(I0), ICRSIZ, IUHF)
      CALL ACES_FIN
      STOP
      END
