










C
C VECTORIZED CC/MBPT PROGRAM
C
C WRITTEN BY J.F. STANTON, J. GAUSS, J.D. WATTS AND W.J. LAUDERDALE
C TD-CC EXTENSION ADDED BY P.G. SZALAY (ROUTINES WITH 'OS')
C




































































































































































































      SUBROUTINE VCC(ICORE,ICRSIZ,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISTSZ,RLECYC,T1_CYCLE
      INTEGER POP(8,2),VRT(8,2)
      LOGICAL BPRINT, NO_REDUNDANT
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,QCISD,LINCC
      LOGICAL CICALC,NONHF,TRIPIT,TRIPNI,TRIPNI1,T3STOR,PRESNT,INIT,
     &        DORESET,UCC,BRUECK,BRKCNV,READT,CCN,CC2
      LOGICAL ROHF4,ITRFLG,OSFLAG,INITL1,SVD
      LOGICAL CIS,EOM,ACT_SPC_CC, LEOM_CC, PEOM_X
      LOGICAL DO_HBAR_4LCCSD,SING_4TCC
      LOGICAL EXTRNL_CCSD
      LOGICAL pCCD,pCCDS,pCCDTS,pCCDTSD,OO_CC
      CHARACTER*4 ACT
      DIMENSION ECORR(3)
      INTEGER ICORE(ICRSIZ),ICRSIZ,IUHF,I0
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /FLAGS/  IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC,CC2
      COMMON /TRIPLES/TRIPNI,TRIPNI1,TRIPIT,T3STOR
cAP - 500 is the maximum number of CC iterations (not basis functions)
      COMMON /ENERGY/ ENERGY(500,2),IXTRLE(500)
      COMMON /LINEAR / LINCC,CICALC
      COMMON /NHFREF/ NONHF
      COMMON /ROHF/ ROHF4,ITRFLG
      COMMON /CORENG/ ELAST
      COMMON /BRUECKNER/ BRUECK
      COMMON /SYMLOC/ ISYMOFF(8,8,25)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP,VRT,NT(2),NFMI(2),NFEA(2)
      COMMON /AUXIO / DISTSZ(8,100),NDISTS(8,100),INIWRD(8,100),LNPHYR,
     &                NRECS,LUAUX
      COMMON /T3IOOF/ IJKPOS(8,8,8,2),IJKLEN(36,8,4),IJKOFF(36,8,4),
     &                NCOMB(4)
      COMMON /EXCITE/ CIS,EOM


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


      integer ncycle
      common /ncycle/ncycle

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 

C
      EQUIVALENCE (IFLAGS(2),METHOD)
 1000 FORMAT(/,77('-'),/,32X,'Entering xvcc',/,77('-'),/)
      WRITE(6,1000)
C
      call bwread(.true.,iuhf)
      if (isbwcc) call bwcc(.false.)
C
      iref = 1
      ibwpass = 0
      I0 = 1

      bPrint = iflags(1).ge.10
      NO_REDUNDANT = iflags2(155).ne.0

C Set the Damp parameter. The damping of CC amplitudes during the
C iterations is usefull to improve the convergence in difficult cases.
C The default of this is set to 1.0

      DAMP_PARAMETER = IFLAGS2(170)*0.010D0
      EXTRNL_CCSD    = IFLAGS(2) .EQ. 50 .OR. 
     &                 IFLAGS(2) .EQ. 51
      pCCD    = IFLAGS(2)  .EQ. 52
      pCCDS   = IFLAGS(2)  .EQ. 55
      pCCDTS  = IFLAGS(2)  .EQ. 56
      pCCDTSD = IFLAGS(2)  .EQ. 58
      OO_CC   = (IFLAGS(2) .EQ. 52 .OR.
     &           IFLAGS(2) .EQ. 53 .OR.
     &           IFLAGS(2) .EQ. 54)
      IF (PCCD) THEN
         METHOD=8
      ELSEIF (pCCDS) THEN
         METHOD=10
      ELSEIF (pCCDTS) THEN
         METHOD=10
      ELSEIF (pCCDTSD) THEN
         METHOD=10
      ENDIF 
      CALL GETREC(0,"JOBARC","EXTRNLCC",LENGTH_EXTREC,1)
      WRITE(6,*) '@VCC-DEBUG: LENGTH_EXTREC=',LENGTH_EXTREC,
     &           ' METHOD=',METHOD
C
C This is an additional argument for GENINT. This is false except
c when Hbar(mb,ej) is constructed for linear CCSD. 05//2014, Ajith
C perera. 
C
      DO_HBAR_4LCCSD = .FALSE.

      Write(*,*)
      write(*,"(3x,a,F6.4)") "The CC damp parameter = ", DAMP_PARAMETER
      Write(*,*)

      CALL ZERO(ENERGY,1000)
      CALL IZERO(IXTRLE,500)
      CALL SETMET_VCC(pCCD,pCCDS,pCCDTS,pCCDTSD)
      CALL SETLOG

      ICONVG=1

      if (method .eq. 10 .or. method .eq. 8
     $     .or. method .eq. 22 .or. method  .eq. 13
     $     .or. method .eq. 14 .or. method  .eq. 16
     $     .or. method .eq. 18 .or. method  .eq. 50
     $     .or. method .eq. 51) then 
         coulomb = .False. 
         call parread(iuhf)
         if (ispar) then
           write(6,*) '  Parameterized CC calculation is performed'
           Write(6,*)
           write(6,2010) paralpha
           write(6,2011) parbeta
           write(6,2012) pargamma
           write(6,2013) pardelta 
           write(6,2014) parepsilon
 2010      format(' PCCSD   alpha parameter : ', F14.6)
 2011      format(' PCCSD    beta parameter : ', F14.6)
 2012      format(' PCCSD   gamma parameter : ', F14.6)
 2013      format(' PCCSD   delta parameter : ', F14.6)
 2014      format(' PCCSD epsilon parameter : ', F14.6)
           if (coulomb) Write(6,"(a,a)") " The Coulomb integrals are ",
     $                    "used in W(mbej) intermediate."
           write(6,*)
           Fae_scale    = (Paralpha - 1.0D0)
           Fmi_scale    = (Parbeta  - 1.0D0)
           Wmnij_scale  = Pargamma
           Wmbej_scale  = Pardelta
           Gae_scale    = Paralpha 
           Gmi_scale    = Parbeta 
         else
           write(6,*) '  Regular CC calculation is performed'
           write(6,*)
           Fae_scale    = 0.0D0
           Fmi_scale    = 0.0D0
           Wmnij_scale  = 1.0D0
           Wmbej_scale  = 1.0D0
           Gae_scale    = 1.0D0
           Gmi_scale    = 1.0D0
         endif
      endif

c
      if (ispar .and. no_redundant .and.
     $     abs(parbeta-1.0d0) .gt. 1.0d-6) then
         write(6,*)' PCCSD is not supported with NO_REDUNDANT option'
         write(6,*)' on Need to fix subroutine WMBEJNR for this to work'
         write(6,*)' PCCSD Calculation will still work if parbeta=1.0'
         call aces_exit(1)
      endif
c
      IF(NO_REDUNDANT) THEN
         CALL NONRED_INIT(IUHF)
      ENDIF
C
      CALL GETREC(-1,'JOBARC','OSCALC',1,IOS)
      OSFLAG = (IOS.NE.0)
      BRUECK = (IFLAGS(22).EQ.1 .OR. IFLAGS2(179) .EQ. 2) 
      ITRFLG = .TRUE.
      ROHF4  = .FALSE.

c     INITIALIZATION FOR OPENSHELL SINGLET CALCULATION
      IF (OSFLAG) CALL OSINIC(ICORE,ICRSIZ,IUHF,I0)
C
c Nevin added to insure maxcor is aligned from the top as well as the bottom
c      MAXCOR=ICRSIZ
      MAXCOR=ICRSIZ-iand(icrsiz,1)

      IF (LENGTH_EXTREC .GT.0) THEN
         call check_ints(icore(i0),maxcor,iuhf,.false.)
         write(6,"(a)") " The T1/T2 old"
         call check_t2in(icore(i0),Maxcor,Iuhf) 
         Write(6,"(3x,2a)")"Starting a second-phase of the external CC",
     +                     " calculation"
         IF (.NOT. NONHF) NONHF = (IFLAGS(38) .GT.0)
         CALL ADD_EXTERNAL_T3_TO_H12(ICORE(I0),MAXCOR,IUHF,NONHF)
         CALL SET_INTMS_2ZERO_4EXTRNL_CC(ICORE(I0),MAXCOR,IUHF)
         call check_ints(icore(i0),maxcor,iuhf,.false.)
         write(6,*)
      ENDIF

C INITIALIZE DOUBLES LISTS
      IF (METHOD.GT.1 .AND. LENGTH_EXTREC .LT. 0) 
     &    CALL INMBPT(ICORE(I0),MAXCOR,IUHF)

C  INITIALIZE LISTS FOR UCC CALCULATION
      IF (UCC) CALL UCCLST(ICORE(I0),MAXCOR,IUHF)
C
C INITIALIZE SINGLE LISTS
      IF (METHOD.GE.3.OR.NONHF) CALL INSING(METHOD,IUHF)
C ALSO DO SINGLES IF pCCD and BRUECKNER, Ajith Perera,10/2021
      IF (pCCD .AND. BRUECK) CALL INSING(10,IUHF)
C
C-----------------------------------------------------------------------
C     Create and initialize triples lists if these are needed.
C-----------------------------------------------------------------------
      T3STOR = .FALSE.
      METACT = 0
      IF (TRIPNI.OR.TRIPNI1.OR.TRIPIT) THEN
         T3STOR = METHOD.EQ.12.OR.
     &            METHOD.EQ.17.OR.
     &            METHOD.EQ.18.OR.
     &           (METHOD.GE.26.AND.METHOD.LE.31).OR.
     &           (TRIPIT.AND.IFLAGS(38).GT.0.AND.IFLAGS(39).EQ.0)
c 3/14/97     &   .OR.
c 3/14/97     &        (METHOD.EQ.14.AND.IFLAGS(87).EQ.3).OR.
c 3/14/97     &        (METHOD.EQ.16.AND.IFLAGS(87).EQ.3).OR.
c 3/14/97     &        (METHOD.EQ.33.AND.IFLAGS(87).EQ.3).OR.
c 3/14/97     &        (METHOD.EQ.34.AND.IFLAGS(87).EQ.3)
c11/ 3/96     &        (METHOD.EQ.34.AND.IFLAGS(87).EQ.3).OR.
c11/ 3/96     &        (METHOD.EQ.22.AND.IFLAGS(87).EQ.3 .AND.
c11/ 3/96     &           (IFLAGS2(124).EQ.4 .OR. IFLAGS2(124).EQ. 6 .OR.
c11/ 3/96     &            IFLAGS2(124).EQ.8 .OR. IFLAGS2(124).EQ.10) )
         IF (METHOD.EQ.29.OR.METHOD.EQ.30) THEN
            METACT = METHOD
            METHOD = 18
         END IF
C         T3STOR = .TRUE.
         IF (T3STOR) THEN
 1010       FORMAT('  @VCC-I, Triples lists will be created. ')
            WRITE(6,1010)
C           PRESNT indicates whether file already exists or not.
C           At this stage in a calculation it should not exist.
            PRESNT = .FALSE.
            INIT   = .TRUE.
            IRREPX = 1
            NLIST  = 8
            IF (IFLAGS(87).EQ.3.AND.IFLAGS(91).GT.0) THEN
               IF (METHOD.EQ.13) NLIST = 12
               IF (METHOD.EQ.14) NLIST = 12
               IF (METHOD.EQ.16) NLIST = 12
               IF (METHOD.EQ.33) NLIST = 12
               IF (METHOD.EQ.34) NLIST = 12
               IF (METHOD.EQ.18) NLIST = 16
               IF (METHOD.EQ.22) NLIST =  4
C               IF (METHOD.EQ.22) NLIST =  8
            END IF
            CALL AUXIOI(PRESNT,INIT,IRREPX,NLIST,IUHF)
c        END IF (T3STOR)
         END IF
c     END IF (TRIPNI.OR.TRIPNI1.OR.TRIPIT)
      END IF
C-----------------------------------------------------------------------

C SET UP SOME INFORMATION BEFORE CALCULATION.
C
CJDW  3/20/96 IEVERY is set equal to IFLAGS(13), so that TAMP_SUM
C             keyword works.
C      IEVERY = 5
      IEVERY = IFLAGS(13)
      SING1  = .FALSE.
      READT  = .FALSE.
C
      IF (CC) THEN
         NCYCLE = IFLAGS(7)
         NKEEP  = IFLAGS(12)
         ICONTL = IFLAGS(4)
      END IF
      IF (bPrint) THEN
         WRITE(6,801)
  801    FORMAT(T3,' Initial T amplitudes: ')
         CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,SING1,'T')
         WRITE(6,800)
  800    FORMAT(T3,' Correlation energies computed from initial ',
     &          'T amplitudes: ')
      END IF

      CALL CMPENG(ICORE(I0),MAXCOR,43,0,ECORR,
     &            ENERGY(1,1),ENERGY(1,2),IUHF,1)
      ELAST = ENERGY(1,2)

c ----------------------------------------------------------------------

      IF (.NOT.CC) THEN
CMN
CJDW 10/24/95. Block of MN's code.
C
C GENERATE INTERMEDIATES FOR MBPT(2) BASED EOM CALCULATION; RPA and
C DRPA extensions are added. The T2 for RPA is generated in RPA code
C and written to lists 44, 45, 46. There is no change to those 
C here. Only create all the things that need to set for Hbar
C construction. Ajith Perera, (iflags2(117)=eomref=mbpt(2),
C rpa and drpa) 07/2012.
C
         IF ((IFLAGS2(117).EQ. 2  .OR. IFLAGS2(117).EQ. 4 .OR.
     &        IFLAGS2(117).EQ. 5  .OR. IFLAGS2(117).EQ. 10)
     &       .AND. .NOT. NONHF) THEN
            CALL INMBPT(ICORE(I0),MAXCOR,IUHF)
            CALL INSING(10,IUHF)
            CALL SETLST(ICORE(I0),MAXCOR,IUHF)
            CALL INITIN(ICORE(I0),MAXCOR,IUHF)
            CALL INITSN(ICORE(I0),MAXCOR,IUHF)
            CALL INCOR(I0,ICRSIZ,IUHF)
            MAXCOR = ICRSIZ
            CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')
            CALL GENINT(ICORE(I0),MAXCOR,IUHF,2,DO_HBAR_4LCCSD)
            CALL INIT2(IUHF,PCCD.OR.PCCDS.OR.PCCDTS.OR.PCCDTSD)
C
C  PUT CONVERGED AMPLITUDES AS INITIAL GUESS FOR SOLVING THE
C  LAMBDA EQUATIONS ON THE GAMLAM FILE
C
            DO ISPIN = 3, 3-2*IUHF, -1
               CALL MVSYM(ICORE(I0),MAXCOR,IUHF,1,43+ISPIN,143+ISPIN)
            END DO
            initl1 = .false.
            if (initl1) then
               DO ISPIN=1,IUHF+1
                  CALL UPDMOI(1,NT(ISPIN),ISPIN,190,0,0)
                  CALL UPDMOI(1,NT(ISPIN),ISPIN+2,190,0,0)
               END DO
C              ZERO THE L1 LIST (ALL CYCLES)
               CALL ZERLST(ICORE,NT(1),1,1,1,190)
               IF (IUHF.NE.0) CALL ZERLST(ICORE,NT(2),1,1,2,190)
            end if
C
            IF (IUHF.EQ.0) CALL RESET(ICORE(I0),MAXCOR,IUHF)
            IF (IFLAGS(35).NE.0) THEN
               CALL ACES_AUXCACHE_FLUSH
               CALL ACES_AUXCACHE_RESET
            END IF
            WRITE(6,1020)
            RETURN
c        END IF (IFLAGS2(117).EQ.2)
         ENDIF
C
C Standard EOM-MBPT(2) ends here. Start the PEOM_X. What is it. The
C first-step of it to solve the Brueckner MBPT(2). That is equivalent
C to solving Q2H(bar,T2(1)) = 0.0 Q1Hbar(T2(1),T1)=0.
C 

         PEOM_X = (Iflags2(117)  .EQ. 8 .AND.
     &             Iflags(87)    .NE. 0 .AND.
     &             Iflags(2) .EQ. 1) 

         IF (PEOM_X) THEN

C
C Initialize the singles lists (list 90; 1,2,3,4; 3 and 4 are
C T1 increments). The SETLST form all Wmbej intermediate lists.
C SING1 does not need to be .TRUE. unless T1 guess is read 
C from somwhere else. The INMBPT at this point is not necessary
C (perhaps) but sure does do no harm. 
C
            SING1 = .FALSE.
C
C First Initialize some list that do not get initialized by
C intilization calls (these are normally done INSING)
C      
            DO ISPIN = 1, IUHF+1
               CALL UPDMOI(1,NFMI(ISPIN),ISPIN,  91,0,0)
               CALL UPDMOI(1,NFEA(ISPIN),ISPIN,  92,0,0)
               CALL UPDMOI(1,NT(ISPIN),  ISPIN,  93,0,0)
            END DO

            CALL UPDMOI(1,NT(1),3,90,0,0)
            CALL UPDMOI(1,NT(1),3,93,0,0)
            IF (IUHF.EQ.1) THEN
               CALL UPDMOI(1,NT(2),4,90,0,0)
               CALL UPDMOI(1,NT(2),4,93,0,0)
            ENDIF

            CALL INMBPT(ICORE(I0),MAXCOR,IUHF)
            CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')
            CALL SETLST(ICORE(I0),MAXCOR,IUHF)
C           
            DIFF     = 0.1D0
            T1_CYCLE = 0
C
C Iterate the T1 equation for a fixed T2 untill the desired 
C convergence is achived.
C
            TOL = 10.0**(-IFLAGS(4))

            Write(6,*) TOL, IFLAGS(7)
            DO WHILE (DIFF .GT. TOL .AND.
     &                T1_CYCLE .LE. IFLAGS(7))
                
               call aces_list_memset(3,90,0)
               call aces_list_memset(1,91,0)
               call aces_list_memset(1,92,0)
               call aces_list_memset(1,93,0)
               call aces_list_memset(3,93,0)
      
               call zersym(icore(i0),51)
               call zersym(icore(i0),52)
               call zersym(icore(i0),54)
               call zersym(icore(i0),56)
               call zersym(icore(i0),58)

               if (iuhf .ne. 0) then
                   call aces_list_memset(4,90,0)
                   call aces_list_memset(2,91,0)
                   call aces_list_memset(2,92,0)
                   call aces_list_memset(2,93,0)
                   call aces_list_memset(4,93,0)
                   call zersym(icore(i0),53)
                   call zersym(icore(i0),55)
                   call zersym(icore(i0),57)
                   call zersym(icore(i0),59)
               endif

               CALL COMPUTE_T1_4T2(ICORE(I0),MAXCOR,IUHF,DIFF)
               T1_CYCLE = T1_CYCLE + 1
     
            ENDDO 

            IF  (T1_CYCLE .GT. IFLAGS(7) .OR. 
     &           DIFF .GT. IFLAGS(4)) THEN
                Write(6,"(A,A)") "T1 iterations did not converged",
     &                           " in alloted number of cyclces."
                call ERREX
            ELSE

                Write(6,"(A,1x,i3,1x,A)")"T1 iterations converged in",
     &                           T1_CYCLE, "cyclces."
            ENDIF
C
C It is important to have the Brueckner key-word set. So, the test
C Breuckner is not absolutely necessary. 
C
            IF (BRUECK) THEN
               CALL BRUECKIT(ICORE(I0),MAXCOR,IUHF,BRKCNV,pCCD)
               RETURN
            END IF
            
            IF (BRKCNV) THEN
            CALL INMBPT(ICORE(I0),MAXCOR,IUHF)
            CALL INSING(10,IUHF)
            CALL SETLST(ICORE(I0),MAXCOR,IUHF)
            CALL INITIN(ICORE(I0),MAXCOR,IUHF)
            CALL INITSN(ICORE(I0),MAXCOR,IUHF)
            CALL INCOR(I0,ICRSIZ,IUHF)
            MAXCOR = ICRSIZ
            CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')
            CALL GENINT(ICORE(I0),MAXCOR,IUHF,2,DO_HBAR_4LCCSD)
            CALL INIT2(IUHF,PCCD.OR.PCCDS.OR.PCCDTS.OR.PCCDTSD)

C The FEACONT_MODF, FMICONT_MODF and FMECONT_MODF are simply
C modifications to the orgianal routines (as the name indicates).
C They compute only FT1 + FT2 to T1 equation. The original 
C also compute 

             CALL FEACONT_MODF(ICORE(I0),MAXCOR,IUHF)
             CALL FMICONT_MODF(ICORE(I0),MAXCOR,IUHF)
             CALL FMECONT(ICORE(I0),MAXCOR,IUHF,1)
C 
C This will compute the WT2 to the T1 equation. 


C Built Hbar(p,q) and Hbar(mb,ej) and Hbar(mn,ij) with the 
C MBPT(2) T2. The Hbar(ij,ka), Hbar(ab,ci) and Hbar(ab,cd) is
C handlled in lambda code.
C Hbar(p,q)  :  f(p,q) + T1*W + Tau2*W2 terms.
C Hbar(pq,rs):  W2 + T1*W + Tau*W2 terms (note that during LCCD
C iteration only W2 terms are invloved).
C
C This may not be necessary; but safer than being sorry. All
C the (mb,ej) intermediates are set to zero inside genint. So,
C it is not done here.

               IF (IUHF .NE. 0) THEN
                   CALL ZERSYM(ICORE,51)
                   CALL ZERSYM(ICORE,52)
               ELSE
                   CALL ZERSYM(ICORE,53)
               ENDIF
C
               CALL GENINT(ICORE(I0),MAXCOR,IUHF,2,DO_HBAR_4LCCSD)
               IF (IUHF .EQ. 0) CALL RESET(ICORE(I0),MAXCOR,IUHF)
C
C I do not know whether these are necessry. There is no L1 and L2
C for LCCSD. I believe our EOM codes looks for these lists for
C the left hand vectors. So, this is simply a trick to get the
C left side (may be gradients to work).

               DO ISPIN = 3, 3-2*IUHF, -1
                  CALL MVSYM(ICORE(I0),MAXCOR,IUHF,1,43+ISPIN,60+ISPIN)
                  CALL MVSYM(ICORE(I0),MAXCOR,IUHF,1,43+ISPIN,143+ISPIN)
               ENDDO
               DO ISPIN=1,IUHF+1
                   CALL UPDMOI(1,NT(ISPIN),ISPIN,190,0,0)
                   CALL UPDMOI(1,NT(ISPIN),ISPIN+2,190,0,0)
               END DO

                CALL GETLST(ICORE(I0), 1, 1, 1, 1, 90)
                CALL PUTLST(ICORE(I0), 1, 1, 1, 1, 190)

                IF (IUHF .NE. 0) THEN
                   CALL GETLST(ICORE(I0), 1, 1, 1, 2, 90)
                   CALL PUTLST(ICORE(I0), 1, 1, 1, 2, 190)
               ENDIF
C
         ENDIF 
               
         ENDIF
CMN END

C CALL SPECIAL ROUTINE FOR NONHF PERTURBATION THEORY

         IF (NONHF) THEN

            CALL INMBPT(ICORE(I0),MAXCOR,IUHF)
            CALL ROHFPT(ICORE(I0),MAXCOR,IUHF)
C
C This block handle NONHF EOM-MBPT(2),RPA,DRAP calcualttions
C Ajith Perera, 02/2017.  

         IF ((IFLAGS2(117).EQ. 2  .OR. IFLAGS2(117).EQ. 4 .OR.
     &        IFLAGS2(117).EQ. 5  .OR. IFLAGS2(117).EQ. 10))THEN

            CALL INMBPT(ICORE(I0),MAXCOR,IUHF)
            CALL INSING(10,IUHF)
            CALL SETLST(ICORE(I0),MAXCOR,IUHF)
            CALL INITIN(ICORE(I0),MAXCOR,IUHF)
            CALL INITSN(ICORE(I0),MAXCOR,IUHF)
            CALL INCOR(I0,ICRSIZ,IUHF)
            MAXCOR = ICRSIZ
            CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')
            CALL GENINT(ICORE(I0),MAXCOR,IUHF,2,DO_HBAR_4LCCSD)
            CALL INIT2(IUHF,PCCD.OR.PCCDS.OR.PCCDTS.OR.PCCDTSD)
C
C  PUT CONVERGED AMPLITUDES AS INITIAL GUESS FOR SOLVING THE
C  LAMBDA EQUATIONS ON THE GAMLAM FILE
C
            DO ISPIN = 3, 3-2*IUHF, -1
               CALL MVSYM(ICORE(I0),MAXCOR,IUHF,1,43+ISPIN,143+ISPIN)
            END DO
            initl1 = .false.
            if (initl1) then
               DO ISPIN=1,IUHF+1
                  CALL UPDMOI(1,NT(ISPIN),ISPIN,190,0,0)
                  CALL UPDMOI(1,NT(ISPIN),ISPIN+2,190,0,0)
               END DO
C              ZERO THE L1 LIST (ALL CYCLES)
               CALL ZERLST(ICORE,NT(1),1,1,1,190)
               IF (IUHF.NE.0) CALL ZERLST(ICORE,NT(2),1,1,2,190)
            end if
C
            IF (IUHF.EQ.0) CALL RESET(ICORE(I0),MAXCOR,IUHF)
            IF (IFLAGS(35).NE.0) THEN
               CALL ACES_AUXCACHE_FLUSH
               CALL ACES_AUXCACHE_RESET
            END IF
            WRITE(6,1020)
            RETURN
c        END IF (IFLAGS2(117).EQ.2)
         ENDIF
C
C IN CASE OF MBPT BASED ON BRUECKNER ORBITALS, GET NEW ORBITALS
C BY CALLING BRUECKIT

            IF (BRUECK) THEN
               write(*,5100)
5100           FORMAT('   Residual T1 amplitudes ')
               CALL AMPT1(ICORE(I0),MAXCOR,IUHF)
               CALL BRUECKIT(ICORE(I0),MAXCOR,IUHF,BRKCNV,pCCD)

               IF (BRKCNV) THEN

                  IF (Iflags2(117)  .EQ. 2 .AND.
     &                Iflags(87)    .NE. 0) THEN
                      CALL INMBPT(ICORE(I0),MAXCOR,IUHF)
                      CALL INSING(10,IUHF)
                      CALL SETLST(ICORE(I0),MAXCOR,IUHF)
                      CALL INITIN(ICORE(I0),MAXCOR,IUHF)
                      CALL INITSN(ICORE(I0),MAXCOR,IUHF)
                      CALL INCOR(I0,ICRSIZ,IUHF)
                      MAXCOR = ICRSIZ
CSSS                      CALL DRMOVE(ICORE(I0),MAXCOR,IUHF,0,.FALSE.)
                      CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')
                      CALL GENINT(ICORE(I0),MAXCOR,IUHF,2,
     &                            DO_HBAR_4LCCSD)
                      CALL INIT2(IUHF,PCCD.OR.PCCDS.OR.PCCDTS
     &                           .OR.PCCDTSD)
                      IF (IUHF.EQ.0) CALL RESET(ICORE(I0),MAXCOR,
     &                                          IUHF)
                  ENDIF 
               ENDIF 

            END IF
C
C FOR GRADIENT CALCULATION, INITIALIZE ``GAMLAM-LISTS'' HERE
C
            CALL INCOR(I0,ICRSIZ,IUHF)
            MAXCOR = ICRSIZ
C
            IF (GRAD.AND.METHOD.GT.1) THEN
C
C FOR ROHF-MBPT(3) GRADIENRS, COPY T[2] TO L[2] LISTS
C
C  T[1] LISTS 44,45,46   T[2] LISTS 61,62,63
C  L[1] IS  T[1]         L[2] LISTS 143,144,145
C
C   SINGLES : T[1]: 1,2;90   T[2]:  3,4;90
C             L[1] IS T[1]   L[2]:  1,2;190
C
               CALL INIT2(IUHF,PCCD.OR.PCCDS.OR.PCCDTS.OR.PCCDTSD)
               CALL DRMOVE(ICORE(I0),MAXCOR,IUHF,100,.TRUE.)
            END IF
C
C HF CASES
C
         ELSE IF (MBPT3) THEN
            CALL SETLST(ICORE(I0),MAXCOR,IUHF)
            CALL INCOR(I0,ICRSIZ,IUHF)
            MAXCOR = ICRSIZ
            CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')
            CALL INMBPT(ICORE(I0),MAXCOR,IUHF)
            CALL GENINT(ICORE(I0),MAXCOR,IUHF,0,DO_HBAR_4LCCSD)
            CALL DRE3EN(ICORE(I0),MAXCOR,IUHF,0)
            CALL NEWT2(ICORE(I0),MAXCOR,IUHF)
            CALL CMPENG(ICORE(I0),MAXCOR,60,0,ECORR,ENERGY(2,1),
     &                  ENERGY(2,2), IUHF,0)
            IF (GRAD) CALL FORMT2(ICORE(I0),MAXCOR,IUHF,43,60)
         ELSE IF (MBPT4) THEN
C
C FOURTH-ORDER LOGIC.  COMPUTE D3 AND D4 ON FIRST PASS, THEN GET QUADS.
C
            CALL SETLST(ICORE(I0),MAXCOR,IUHF)
            CALL INCOR(I0,ICRSIZ,IUHF)
            MAXCOR=ICRSIZ
            CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')
            CALL INMBPT(ICORE(I0),MAXCOR,IUHF)
            CALL GENINT(ICORE(I0),MAXCOR,IUHF,0,DO_HBAR_4LCCSD)
            CALL DRE3EN(ICORE(I0),MAXCOR,IUHF,0)
            CALL NEWT2(ICORE(I0),MAXCOR,IUHF)
            CALL GETD4(ICORE(I0),MAXCOR,ENERGY(3,1),IUHF)
            CALL CMPENG(ICORE(I0),MAXCOR,60,0,ECORR,ENERGY(2,1),
     &                  ENERGY(2,2),IUHF,0)
            WRITE(6,4001)ENERGY(2,1)
4001        FORMAT(T3,'E3D = ',F14.10,' a.u.')
4002        FORMAT(T3,'E4D = ',F14.10,' a.u.')
4003        FORMAT(T3,'E4Q = ',F14.10,' a.u.')
C
C  FOR GRADIENT CALCULATION MOVE T2 TO THE LAMBDA LISTS
C
            IF (GRAD) THEN
               CALL INIT2(IUHF,PCCD.OR.PCCDS.OR.PCCDTS.OR.PCCDTSD)
               CALL DRMOVE(ICORE(I0),MAXCOR,IUHF,100,.FALSE.)
            END IF
            WRITE(6,4002)ENERGY(3,1)
            CALL INMBPT(ICORE(I0),MAXCOR,IUHF)
            CALL GENINT(ICORE(I0),MAXCOR,IUHF,1,DO_HBAR_4LCCSD)
            CALL FEACONT(ICORE(I0),MAXCOR,IUHF)
            CALL FMICONT(ICORE(I0),MAXCOR,IUHF)
            CALL DRE3EN(ICORE(I0),MAXCOR,IUHF,1)
            CALL NEWT2(ICORE(I0),MAXCOR,IUHF)
            CALL CMPENG(ICORE(I0),MAXCOR,60,0,ECORR,ENERGY(4,1),
     &                  ENERGY(4,2),IUHF,0)
            WRITE(6,4003)ENERGY(4,1)
            CALL E4S(ICORE(I0),MAXCOR,IUHF,ENERGY(5,1))
            IF (TRIPNI) CALL TRPS(ICORE(I0),MAXCOR,IUHF,ENERGY(6,1),
     &                            .FALSE.,EXTRNL_CCSD)
            IF (IUHF.EQ.0) CALL RESET(ICORE(I0),MAXCOR,IUHF)

c        END IF (NONHF)
         END IF

         CALL MBPTOUT

c     ELSE IF (CC) THEN
      ELSE

      CALL INMBPT(ICORE(I0),MAXCOR,IUHF)
C
C CC CALCULATION LOGIC BEGINS HERE.
C
C SG 7/24/98
C Try to read in a T guess.  If we can, set SING1 to be true.
C
      CALL DRDTGSS(ICORE(I0),MAXCOR/IINTFP,IUHF,'TGUESS  ',0,READT)
C
C The CC2, (eventually CC3..) starts here. The CC2 shares most
C of the CCSD calls. When the ode enters to DRIVE_CCN, it does not
C return back to the main program. Both converged and else is 
C handled within the DRIVE_CCN. The CC2 changes are in, setmet
C setlst, feacont, fmicont, rngprd, t1w1aa, t1w1ab, ladaa1, ladab1,
C q1all routines. Ajith Perera, 04/2015. 
C
      CCN = .FALSE.
      CCN = CC2

      IF (CCN) THEN
         CALL DRIVE_CCN(ICORE(I0),MAXCOR,IUHF,NCYCLE,ICONTL,
     &                  DAMP_PARAMETER,PCCD.OR.PCCDS.OR.PCCDTS
     &                  .OR.PCCDTSD) 
      END IF 
C
C Note that we need to test to see whether we are doing a CCD
C calculation before we turn the SING1 on. This was not tested
C originaly and CCD optimizations were failing. 01/2001, Ajith Perera
C 
      IF ((IFLAGS(2).NE.5.OR.IFLAGS(2).NE.8).AND.READT) SING1 = .TRUE.
      IF (.NOT.SING1) THEN
         CALL GETLST(ICORE(I0),1,1,1,1,90)
         DTMP  = DNRM2(NT(1),ICORE(I0),1)
         SING1 = (DTMP.GT.1.D-7)
C
C If NON-HF flag is set and calc=CCD the above incorrectly turns on the single
C flag.
        
         IF (NONHF .AND. SING1 .AND. 
     &      (IFLAGS(2).EQ.5.OR.IFLAGS(2).EQ.8).OR.OO_CC) THEN
              SING1 = .FALSE.
              CALL IZERO(ICORE(I0),NT(1)) 
              CALL PUTLST(ICORE(I0),1,1,1,1,90)
              IF (IUHF .NE. 0) THEN
                 CALL IZERO(ICORE(I0),NT(2)) 
                 CALL PUTLST(ICORE(I0),1,1,1,2,90)
              ENDIF
         ENDIF  

C The finite field CCD (or no singles methods) response property
C calcs should not turn on the singles even when there are non-diagonal
c elements in fock. A check is not performed on the calc level since
c this block will not get executed for methods that do not include
c singles, Ajith Perera, 01/2016. 
C
         IF (iflags(19) .eq. 1 .and.
     &      (iflags(23)  .ne. 0 .or.
     &       iflags(24)  .ne. 0 .or.
     &       iflags(25)  .ne. 0)) 
     &       SING1 = .FALSE.
         
      END IF

C
C Start the T2 lists by adding <ab||ij>, This is actually a 
C contribution to T2.
C
      WRITE(6,*) '@VCC-TRACE: pre-RNABIJ isytyp(61)=',ISYTYP(1,61),
     &   ISYTYP(2,61),' NO_REDUNDANT=',NO_REDUNDANT
      IF (.NOT. NO_REDUNDANT) CALL RNABIJ(ICORE(I0),MAXCOR,
     &                                    IUHF,'T')
      WRITE(6,*) '@VCC-TRACE: post-RNABIJ isytyp(61)=',ISYTYP(1,61),
     &   ISYTYP(2,61)
C
C GENERATE W LISTS FOR INTERMEDIATES IF THIS IS NOT MBPT(3) OR LCCD.
C
C      IF(METHOD.NE.2.AND.METHOD.NE.5)THEN
         CALL SETLST(ICORE(I0),MAXCOR,IUHF)
C      ENDIF
      WRITE(6,*) '@VCC-TRACE: post-SETLST isytyp(61)=',ISYTYP(1,61),
     &   ISYTYP(2,61)
C**********************************************************************
C
C START OF CC LOOP
C
      RLECYC = 0
      IF (MOD(IFLAGS(21),2).EQ.0) THEN
         CALL DRRLE(ICORE(I0),MAXCOR,IUHF,RLECYC,.FALSE.)
      END IF
      WRITE(6,*) '@VCC-TRACE: post-DRRLE isytyp(61)=',ISYTYP(1,61),
     &   ISYTYP(2,61)
      IF (IFLAGS(21).EQ.1) THEN
         CALL DIISLST(1,IUHF,METHOD.GE.6.AND.METHOD.NE.8)
      END IF
      WRITE(6,*) '@VCC-TRACE: post-DIISLST isytyp(61)=',ISYTYP(1,61),
     &   ISYTYP(2,61)
      CALL INCOR(I0,ICRSIZ,IUHF)
      WRITE(6,*) '@VCC-TRACE: post-INCOR isytyp(61)=',ISYTYP(1,61),
     &   ISYTYP(2,61)
      MAXCOR = ICRSIZ

      icycle = 1
555   continue
C
      RLECYC = RLECYC + 1
C
C IF (ICYCLE.EQ.1) CALL UPDVC(ICORE(I0),ECORR,1,NCYCLE,IUHF)
C FILL T2 INCREMENTS WITH <IJ||AB> INTEGRALS
C AND ZERO THE T1 INCREMENTS (ONLY CCSD AND QCISD)
C JDW 10/24/95. Last argument was removed from INITIN, so remove it here.
C CALL INITIN(ICORE(I0),MAXCOR,IUHF,IFLAGS(93).EQ.2)
C
C Start the T2 lists by adding <ab||ij>, This is actually a
C contribution to T2.
C
      WRITE(6,*) '@VCC-TRACE: pre-INITIN(main) isytyp(61)=',
     &   ISYTYP(1,61),ISYTYP(2,61)
      CALL INITIN(ICORE(I0),MAXCOR,IUHF)
      WRITE(6,*) '@VCC-TRACE: post-INITIN(main) isytyp(61)=',
     &   ISYTYP(1,61),ISYTYP(2,61)

       call Check_T2_VCC(icore(i0),Maxcor,Iuhf)
      IF (METHOD.GE.6.AND.METHOD.NE.8) THEN
C
C Zero out the T1 lists for all cycles.
C
         CALL INITSN(ICORE(I0),MAXCOR,IUHF)
      END IF

c
c calculate the contributions for pccsd
c
CSSS#ifdef _DCC_FLAG

      If (Ispar) Then 
         Write(6,*)
         Write(6,"(a,a)") " ***pCC/DCC FAE and FMI to T2 equation is",
     &                        " constructed***"
          call quad2(icore(i0),maxcor,iuhf, 5,Fmi_scale)
          call fmicont2(icore(i0),maxcor,iuhf)
          call quad3(icore(i0),maxcor,iuhf, 5,Fae_scale)
          call feacont2(icore(i0),maxcor,iuhf)
          call zerlst(icore(i0),nfmi(1),1,1,1,91)
          call zerlst(icore(i0),nfea(1),1,1,1,92)
          if (iuhf .ne. 0) then
             call zerlst(icore(i0),nfmi(2),1,1,2,91)
             call zerlst(icore(i0),nfea(2),1,1,2,92)
          endif 
       Endif 

CSSS#else
CSSS      if (ispar) call pccsd_modification(icore(i0), maxcor, iuhf)
CSSS#endif 
C
C GENERATE W AND F INTERMEDIATES.
C
      IF (LINCC) THEN
         INTTYP = 3
      ELSE
         INTTYP = 2
      END IF
C 
C Generate W(mbej), W(mnij) W(mbej, F(ea), F(ij) and F(me)
C intermediates.  
C
      Write(6,"(a)") " The checksums of all the intms before genint"
      Write(6,"(a,2l)") " Flags singles and non-hf:",SING1,NONHF
      Call checkintms(icore(i0),maxcor,Iuhf,1)
      write(6,"(a)") " The T1/T2 residuals" 
      call Check_T2_VCC(icore(i0),Maxcor,Iuhf) 
C
      CALL GENINT(ICORE(I0),MAXCOR,IUHF,INTTYP,DO_HBAR_4LCCSD)
C
C COMPUTE F INTERMEDIATE CONTRIBUTION TO T2 INCREMENT.
C
C Evaluate P(ab)Sum_e T(ij,ae){F(b,e)-1/2 sum_m T(m,b)F(m,e)}
C contribution to T2 and Sum_e T(i,e)F(a,e) contribution to T1.
C
       call PROBE_T1_90('PRE-FEA-A90/3','PRE-FEA-B90/4',Iuhf)
      CALL FEACONT(ICORE(I0),MAXCOR,IUHF)
       call Check_T2_VCC(icore(i0),Maxcor,Iuhf)
       call PROBE_T1_90('PST-FEA-A90/3','PST-FEA-B90/4',Iuhf)
C
C Evaluate -P(ij)Sum_e T(im,ab){F(mj)-1/2 sum_m T(j,e)F(m,e)}
C contribution to T2 and Sum_e T(m,a)F(m,i) contribution to T1.
C
      CALL FMICONT(ICORE(I0),MAXCOR,IUHF)
       call Check_T2_VCC(icore(i0),Maxcor,Iuhf)
       call PROBE_T1_90('PST-FMI-A90/3','PST-FMI-B90/4',Iuhf)
C
      IF ((METHOD.GT.9.AND.SING1).OR.(METHOD.EQ.6.AND.SING1)) THEN

C Evaluate Sum_e T(im,ae)F(m,e) contribution to T1.
C
       CALL FMECONT(ICORE(I0),MAXCOR,IUHF,1)
       IF (IUHF.NE.0) CALL FMECONT(ICORE(I0),MAXCOR,
     &                                  IUHF,2)
C
       call Check_T2_VCC(icore(i0),Maxcor,Iuhf)
       call PROBE_T1_90('PST-FME-A90/3','PST-FME-B90/4',Iuhf)
C
C Evaluate -P(ab)Sum_m T(m,a)<mb||ij> contribution to T2.
C
         CALL T1INT2A(ICORE(I0),MAXCOR,IUHF)
C
C Evaluate +P(ij)Sum_m T(i,e)<ab||ej> contribution to T2.
C
         CALL T1INT2B(ICORE(I0),MAXCOR,IUHF)
C
C Evaluate -Sum_nf T(n,f)<na||if> contribution to T1.
C
         CALL T1INT1(ICORE(I0),MAXCOR,IUHF,1)
         IF (IUHF.NE.0) CALL T1INT1(ICORE(I0),MAXCOR,
     &                             IUHF,2)
C
       call Check_T2_VCC(icore(i0),Maxcor,Iuhf)
       call PROBE_T1_90('PST-T1I-A90/3','PST-T1I-B90/4',Iuhf)
      END IF
C 
      IF (UCC) THEN
         CALL T1INT2A(ICORE(I0),MAXCOR,IUHF)
         CALL T1INT2B(ICORE(I0),MAXCOR,IUHF)
      END IF
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
       call Check_T2_VCC(icore(i0),Maxcor,Iuhf)
      CALL DRE3EN(ICORE(I0),MAXCOR,IUHF,0)
       call PROBE_T1_90('PST-DR3-A90/3','PST-DR3-B90/4',Iuhf)
C
C DO TdegerWT CONTRIBUTION FOR UCC AND XCC
C 
      IF (UCC) CALL TTWT(ICORE(I0),MAXCOR,IUHF)
C
C NOW DO T3->T2 AND T3->T1 FOR ITERATIVE TRIPLE METHODS
C
      IF (TRIPIT) THEN
         CALL TRPINT(ICORE(I0),MAXCOR,IUHF)
         CALL TRPS(ICORE(I0),MAXCOR,IUHF,ZJUNK,.FALSE.,EXTRNL_CCSD)
      END IF
C 
      IF (OSFLAG) CALL OS(ICORE(I0),MAXCOR,IUHF,W)

      IF (LENGTH_EXTREC .GT. 0) THEN
          CALL ADD_T3IN2T2_ABIJ(ICORE(I0),MAXCOR,IUHF)
      ENDIF 
C
C COMPUTE THE CONTRIBUTION OF DOUBLES TO T1 AND
C DENOMINATOR WEIGHT THE NEW INCREMENTS
C (ONLY CCSD AND QCISD METHODS)
C 
      IF ((METHOD.GE.9.OR.METHOD.EQ.6.OR.METHOD.EQ.7)
     &    .AND. .NOT. OO_CC) THEN
C
C Evaluate the following three contributions to T1:
C
C   -1/2 Sum_mef T(im,ef)<ma||ef>  (T2T1AA1, T2T1AB1)
C   -1/2 Sum_men T(nm,ei)<nm||ei>  (T2T1AA2, T2T1AB2)
C
       call Check_T2_VCC(icore(i0),Maxcor,Iuhf)
C Also, Do the T1 = T1/{f(i,i) - f(a,a)} to get a new T1
C
         CALL E4S(ICORE(I0),MAXCOR,IUHF,EDUMMY)
       call PROBE_T1_90('PST-E4S-A90/3','PST-E4S-B90/4',Iuhf)
      END IF
      CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,SING1,'T')

C Also the call E4S for pCCD (OO-CCD) and Brueckner. This will do T2 into T1
C terms. Ajith Perera,10/2021. 

       IF (OO_CC .AND. BRUECK) THEN
          CALL E4S(ICORE(I0),MAXCOR,IUHF,EDUMMY)
       ENDIF 
C
C If calculation is set to extrnl_ccd (Iflags(2)=50) Then do not
C include singles

      IF (IFLAGS(2) .EQ. 50) THEN
         CALL ZERLST(ICORE(I0),NT(1),1,1,3,90)
         IF (IUHF.NE.0) CALL ZERLST(ICORE(I0),NT(2),1,1,4,90)
      ENDIF 
C
C DENOMINATOR WEIGHT T2 INCREMENTS TO FORM NEW T2.
C 
C Do the T2 = T2/{f(i,i) + f(j,j) - f(a,a) - f(b,b)}
C to get a new T2.
C
      IF (pCCD .OR. pCCDS .OR. pCCDTS) THEN
C Ajith Perera, 10/2021.
         Write(6,"(20x,a)") " ----------Warning------------"
         Write(6,"(2a)")  " A pCCD or pCCD like calculation is being",
     +                    " performed and the off-diagonal"
         Write(6,"(a)")   " blocks of T2 is set to zero"
         Write(6,"(20x,a)") " -----------------------------"
         call pccd_reset_vcc(icore(i0),Maxcor,Iuhf,63)
      ENDIF 
       call Check_T2_VCC(icore(i0),Maxcor,Iuhf)
      CALL NEWT2(ICORE(I0),MAXCOR,IUHF)

       call Check_T2_VCC(icore(i0),Maxcor,Iuhf)
       CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,SING1,'T')
       WRITE(6,*) '@VCC-GAP-DEBUG: after AMPSUM'
       call Check_T2_VCC(icore(i0),Maxcor,Iuhf)

      IF (DAMP_PARAMETER .NE. 1.0D0) THEN
         CALL DAMP_CC_RESIDUAL(ICORE(I0),MAXCOR/IINTFP,IUHF,SING1,
     &                         44,61,90,90,DAMP_PARAMETER)
      ENDIF
       WRITE(6,*) '@VCC-GAP-DEBUG: after DAMP_CC_RESIDUAL block'
       call Check_T2_VCC(icore(i0),Maxcor,Iuhf)
C
C WRITE OUT RMS AND MAX DIFFERENCES FOR NEW T2 AND PUT
C NEW ERROR VECTOR
C
      CALL DRTSTS(ICORE(I0),MAXCOR,ICYCLE,IUHF,ICONVG,ICONTL,
     &            SING1,0,'T')
       WRITE(6,*) '@VCC-GAP-DEBUG: after DRTSTS'
       call Check_T2_VCC(icore(i0),Maxcor,Iuhf)
       WRITE(6,*) '@VCC-SING1-DEBUG: ICYCLE=',ICYCLE,' SING1=',SING1,
     &            ' READT=',READT,' NONHF=',NONHF
      CALL CMPENG(ICORE(I0),MAXCOR,60,2,ECORR,ENERGY(ICYCLE+1,1),
     &            ENERGY(ICYCLE+1,2),IUHF,1)
      WRITE(6,*) '@VCC-ENERGY-DEBUG: ICYCLE=',ICYCLE,
     &           ' E1=',ENERGY(ICYCLE+1,1),' E2=',ENERGY(ICYCLE+1,2)
      CALL FLUSH(6)
C
      IF (OSFLAG) THEN
         WRITE(6,81)ENERGY(ICYCLE+1,1)+W
 81      FORMAT(T3,' The total correlation energy for the singlet ',
     &          'is ',F15.12,' a.u.')
         WRITE(6,82)ENERGY(ICYCLE+1,1)-W
 82      FORMAT(T3,' The total correlation energy for the triplet ',
     &          'is ',F15.12,' a.u.')
         ENERGY(ICYCLE+1,1)=ENERGY(ICYCLE+1,1)+W
         ENERGY(ICYCLE+1,2)=ENERGY(ICYCLE+1,2)+W
      END IF
      ELAST = ENERGY(ICYCLE+1,2)

C IF ICONVG IS EQUAL ZERO, CONVERGENCE HAS BEEN ACHIEVED, EXIT
C VIA FINISH
C
      IF (IFLAGS(21).EQ.1 .OR. IFLAGS2(179) .EQ.2) THEN
         CALL DODIIS0(ICORE(I0),MAXCOR/IINTFP,IUHF,1,ICYCLE,
     &                ICONVG,ICONTL,SING1,44,61,90,0,90,2,70,
     &                '     ',DAMP_PARAMETER)
      END IF

CJDW 10/24/95.
C    Do not update if we have converged. This allows us to maintain a
C    correspondence between the amplitudes and existing Hbar elements.
C    This is important for MN's finite-order EOM stuff.
C
      IF (ICONVG.NE.0) THEN
         CALL DRMOVE(ICORE(I0),MAXCOR,IUHF,0,SING1)
         IF (.NOT. NO_REDUNDANT) CALL RNABIJ(ICORE(I0),MAXCOR,
     &                                       IUHF,'T')
      END IF
CJDW END

      IF (ICONVG.EQ.0) THEN
         CALL CMPENG(ICORE(I0),MAXCOR,43,0,ECORR,ENERGY(ICYCLE+1,1),
     &               ENERGY(ICYCLE+1,2),IUHF,1)
         IF (OSFLAG) THEN
            WRITE(6,81)ENERGY(ICYCLE+1,1)+W
            WRITE(6,82)ENERGY(ICYCLE+1,1)-W
            ENERGY(ICYCLE+1,1)=ENERGY(ICYCLE+1,1)+W
            ENERGY(ICYCLE+1,2)=ENERGY(ICYCLE+1,2)+W
         END IF
      END IF

C SG 7/21/98 Check if we are ready for a Brueckner step
      IF (BRUECK .AND. (ICYCLE.GT.1)) THEN
         CALL CHKBKSTP(ICORE(I0), MAXCOR/IINTFP, IUHF, ICONVG)
      END IF

      IF (ICONVG.EQ.0) THEN

C  CALCULATE THE UCC AND XCC ENERGY CONTRIBUTION
C   !!!! NOTE THAT THIS ROUTINE REWRITES WMBEJ LISTS 54-59 !!!!
C   !!!! NOTE THAT THIS ROUTINE REWRITES T2INC LISTS 41-46 !!!!
         IF (UCC) CALL EUCC(ICORE(I0),MAXCOR,IUHF,ENERGY(ICYCLE+1,1))

CJDW 3/20/96. Print final amplitudes.
   
         CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,SING1,'T')
C
C  PUT CONVERGED AMPLITUDES AS INITIAL GUESS FOR SOLVING THE
C  LAMBDA EQUATIONS ON THE GAMLAM FILE
         WRITE(6,4000)
 4000    FORMAT(T3,' The CC iterations have converged.',/)

         IF ((GRAD.AND..NOT.UCC)) THEN
C
CSSS     &       (IFLAGS(87).EQ.3 
CSSS     &       .AND. (TRIPIT.OR.TRIPNI1)) THEN
C Modified for EOM-CCSD(T) variants, 03/2014, Ajith Perera

            CALL INIT2(IUHF,PCCD.OR.PCCDS.OR.PCCDTS.OR.PCCDTSD)
            CALL DRMOVE(ICORE(I0),MAXCOR,IUHF,100,SING1)
         END IF
C
         IF (IFLAGS(87).EQ.3 .OR. IFLAGS(87).EQ.11) THEN
            CALL INIT2(IUHF,PCCD.OR.PCCDS.OR.PCCDTS.OR.PCCDTSD)
            CALL DRMOVE(ICORE(I0),MAXCOR,IUHF,100,SING1)
         ENDIF

         IF (IUHF.NE.0) CALL S2PROJ(ICORE(I0),MAXCOR,IUHF,SING1)

C-----------------------------------------------------------------------
C     In EOM triples calculations, RESET is sometimes called by RESRNG.
C     Determine if we need to call RESET here or not.
C-----------------------------------------------------------------------
         DORESET = IFLAGS(87).NE.3 .OR.
     &            (IFLAGS(87).EQ.3 .AND. IFLAGS(2).NE.16 .AND.
     &                                   IFLAGS(2).NE.18 .AND.
     &                                   IFLAGS(2).NE.33 .AND.
     &                                   IFLAGS(2).NE.34 .AND.
     &           .NOT. (IFLAGS(2).EQ.22 .AND. IFLAGS2(124).GE.5))
         IF (DORESET) THEN
            IF(IUHF.EQ.0.AND..NOT.UCC)CALL RESET(ICORE(I0),MAXCOR,IUHF)
         END IF
C
C Construct the F(ae),F(mi),W(mnij),W(mbej) lists for pCC

         If (Ispar) Then
            Write(6,"(a)") " Constructing pCC lists post-CC steps"
            CALL POST_DCC_MODS(ICORE(I0), MAXCOR, IUHF) 
         Endif 
CSSS         IF (IFLAGS(87).EQ.3.AND.(TRIPIT.OR.TRIPNI)) THEN
C Modified for EOM-CCSD(T) variants, 03/2014, Ajith Perera
C
         IF (IFLAGS(87).EQ.11) THEN

            CALL RESRNG(ICORE(I0),MAXCOR,IUHF)
         END IF
C-----------------------------------------------------------------------

         CALL FINISH(ICYCLE+1,pCCD,PCCDS,PCCDTS,PCCDTSD)

C SG 7/21/98 
         CALL DDMPTGSS(ICORE(I0), MAXCOR/IINTFP, IUHF, 0, 'TGUESS  ')
         IF (BRUECK) CALL BRUECKIT(ICORE(I0),MAXCOR,IUHF,BRKCNV,pCCD)
         IF (TRIPNI .AND. (.NOT.BRUECK .OR. BRKCNV)) THEN
C
C Lets check whether this is a NON-HF CCSD(T) calculation with
C semi-canonical orbitals. If that is the case, then take the
C time-honored path. Otherwise, iterate the T3 equation (non-HF CCSDT-
C like for converged T2). 

           IF (.NOT. TRIPIT .AND. 
     &         (IFLAGS(39) .EQ. 1  .AND.
     &          IFLAGS(38)    .EQ. 1)  .OR .
     &         (IFLAGS(39) .EQ. 0  .AND.
     &          IFLAGS(38)    .EQ. 0)) THEN

CSSS               Write(*,*) "Enter T3ITER"
CSSS               CALL T3_ITER(ICORE(I0),MAXCOR,IUHF)
CSSS               CALL AUXIOO(.TRUE.)
               CALL TRPS(ICORE(I0),MAXCOR,IUHF,EDUMMY,.FALSE.,
     &                   EXTRNL_CCSD)
           ELSE
               CALL T3_ITER(ICORE(I0),MAXCOR,IUHF, EXTRNL_CCSD)
               CALL AUXIOO(.TRUE.)
           END IF

         ENDIF 
C-----------------------------------------------------------------------
C     Close the T3 file. Delete it except in certain noniterative meth-
C     ods and in EOMEE triples methods.
C-----------------------------------------------------------------------
         IF (T3STOR) THEN
            IF ((METACT.EQ.29).OR.(METACT.EQ.30)) METHOD = METACT
            IF ((METHOD.EQ.12).OR.(METHOD.GE.26.AND.METHOD.LE.31).OR.
     &          (IFLAGS(87).EQ.3.AND.(TRIPIT.OR.TRIPNI))) THEN
               CALL AUXIOO(.FALSE.)
            ELSE
               CALL AUXIOO(.TRUE.)
            END IF
         END IF

C-----------------------------------------------------------------------
CJDW 10/24/95. Write energy and coordinates in POLYRATE jobs.
C-----------------------------------------------------------------------
         IF (IFLAGS2(113).GT.0) THEN
            IONE=1
            CALL GETREC(20,'JOBARC','NATOMS  ',IONE,NATOMS)
            CALL POLYPRT0(NATOMS,IINTFP,IFLAGS(1),
     &                    .TRUE.,.TRUE.,.FALSE.,.FALSE.)
         END IF
C-----------------------------------------------------------------------
         IF (IFLAGS(35).NE.0) THEN
            CALL ACES_AUXCACHE_FLUSH
            CALL ACES_AUXCACHE_RESET
         END IF

C
C For linearized coupled cluster methods, the lambda is T^(t).
C We can get the rest of the modules (EOM) and gradients to
C work by replacing L1 and L2 with converged T1 and T2.
C The density code also needs to be modified.
C  Ajith Perera, 04/2014. 

         LEOM_CC = ((Iflags(2) .EQ. 5 .OR.
     &              Iflags(2) .EQ. 6) .AND.
     &              Iflags(87)    .NE. 0)

         IF (LEOM_CC) THEN

            IF (Iflags(2) .EQ. 5) THEN
C
C This is for LCCD. There are no singles (T1 and L1). The
C Hbar(p,q)  :  f(p,q) + T2*W2 terms.
C Hbar(pq,rs):  W2 + T2*W2 terms (note that during LCCD 
C iteration only W2 terms are invloved).
C
               DO ISPIN = 1, IUHF+1
                  CALL UPDMOI(1,NFMI(ISPIN),ISPIN, 91,0,0)
                  CALL UPDMOI(1,NFEA(ISPIN),ISPIN, 92,0,0)
                  CALL UPDMOI(1,NT(ISPIN),ISPIN,90,0,0)
                  CALL UPDMOI(1,NT(ISPIN),ISPIN+2,90,0,0)
                  CALL UPDMOI(1,NT(ISPIN),ISPIN,190,0,0)
                  CALL UPDMOI(1,NT(ISPIN),ISPIN+2,190,0,0)
                  CALL UPDMOI(1,NT(ISPIN), ISPIN,  93,0,0)

                  call aces_list_memset(ispin,91,0)
                  call aces_list_memset(ispin,92,0)
                  call aces_list_memset(ispin,93,0)
                  call aces_list_memset(ispin,90,0)
                  call aces_list_memset(ispin+2,90,0)
                  call aces_list_memset(ispin,190,0)
                  call aces_list_memset(ispin+2,190,0)
               END DO
C
C This may not be necessary; but safer than being sorry.All
C the (mb,ej) intermediates are set to zero inside genint, so
C it is not done here.
C
               IF (IUHF .NE. 0) THEN
                   CALL ZERSYM(ICORE,51)
                   CALL ZERSYM(ICORE,52)
               ELSE
                   CALL ZERSYM(ICORE,53)
               ENDIF
C
C Built Hbar(p,q) and Hbar(mb,ej) and Hbar(mn,ij) with converged 
C LCCD T2. The Hbar(ij,ka), Hbar(ab,ci) and Hbar(ab,cd) is 
C handlled in lambda code. The T2 is copied to L2 ( this is necessary 
C since L2 =T2 for LCCD). 
C
               CALL GENINT(ICORE(I0),MAXCOR,IUHF,2,DO_HBAR_4LCCSD)
               IF (IUHF .EQ. 0) CALL RESET(ICORE(I0),MAXCOR,IUHF)
C
C Move T2 list to L2 lists for LCCD. I believe our EOM codes 
C looks for these lists for the left hand vectors. 
C So, this is simply a trick to get the left side (may be gradients to work).
C

               DO ISPIN = 3, 3-2*IUHF, -1
                  CALL MVSYM(ICORE(I0),MAXCOR,IUHF,1,43+ISPIN,60+ISPIN)
                  CALL MVSYM(ICORE(I0),MAXCOR,IUHF,1,43+ISPIN,143+ISPIN)
               ENDDO
 
             ELSE IF (Iflags(2) .EQ. 6) THEN

C Built Hbar(p,q) and Hbar(mb,ej) and Hbar(mn,ij) with converged
C LCCD T2. The Hbar(ij,ka), Hbar(ab,ci) and Hbar(ab,cd) is
C handlled in lambda code. 
C Hbar(p,q)  :  f(p,q) + T1*W + Tau2*W2 terms.
C Hbar(pq,rs):  W2 + T1*W + Tau*W2 terms (note that during LCCD
C iteration only W2 terms are invloved).
C
C This may not be necessary; but safer than being sorry. All
C the (mb,ej) intermediates are set to zero inside genint. So, 
C it is not done here. 

               IF (IUHF .NE. 0) THEN
                   CALL ZERSYM(ICORE,51)
                   CALL ZERSYM(ICORE,52)
                   call aces_list_memset(2,91,0)
                   call aces_list_memset(2,92,0)
                   call aces_list_memset(2,93,0)
               ELSE
                   CALL ZERSYM(ICORE,53)
                   call aces_list_memset(1,91,0)
                   call aces_list_memset(1,92,0)
                   call aces_list_memset(1,93,0)
               ENDIF
C
C Before calling genint, we need to save <mn||ij> integrals to built
C Hbarijka for LinCCSD. Also, Hbarijka needs f(a,i) instead of F(i,e)
C (for NONHF) but f(a,i) lists are  untouched in genint.
C
               IRREPX = 1
               IMODE  = 0
               IF (IUHF .EQ. 1) THEN

                  LENGTH_MNIJAA = IDSYMSZ(IRREPX,ISYTYP(1,11),
     &                                    ISYTYP(2,11))
                  LENGTH_MNIJBB = IDSYMSZ(IRREPX,ISYTYP(1,12),
     &                                    ISYTYP(2,12))
                  LENGTH_MNIJAB = IDSYMSZ(IRREPX,ISYTYP(1,13),
     &                                    ISYTYP(2,13))
                  CALL INIPCK(IRREPX,3,3,111,IMODE,0,1)
                  CALL GETALL(ICORE(I0),LENGTH_MNIJAA,IRREPX,11)
                  CALL PUTALL(ICORE(I0),LENGTH_MNIJAA,IRREPX,111)
                  CALL INIPCK(IRREPX,4,4,112,IMODE,0,1)
                  CALL GETALL(ICORE(I0),LENGTH_MNIJBB,IRREPX,12)
                  CALL PUTALL(ICORE(I0),LENGTH_MNIJBB,IRREPX,112)
                  CALL INIPCK(IRREPX,14,14,113,IMODE,0,1)
                  CALL GETALL(ICORE(I0),LENGTH_MNIJAB,IRREPX,13)
                  CALL PUTALL(ICORE(I0),LENGTH_MNIJAB,IRREPX,113)
               ELSE
                  LENGTH_MNIJAB = IDSYMSZ(IRREPX,ISYTYP(1,13),
     &                                           ISYTYP(2,13))
                  CALL INIPCK(IRREPX,14,14,113,IMODE,0,1)
                  CALL GETALL(ICORE(I0),LENGTH_MNIJAB,IRREPX,13)
                  CALL PUTALL(ICORE(I0),LENGTH_MNIJAB,IRREPX,113)
               ENDIF 

               DO_HBAR_4LCCSD = .TRUE.
C
               CALL GENINT(ICORE(I0),MAXCOR,IUHF,2,DO_HBAR_4LCCSD)
               IF (IUHF .EQ. 0) CALL RESET(ICORE(I0),MAXCOR,IUHF)

C Move T2 list to L2 lists for LCCD. I believe our EOM codes 
C looks for these lists for the left hand vectors. 
C So, this is simply a trick to get the left side (may be gradients to work).
C
               DO ISPIN = 3, 3-2*IUHF, -1
                  CALL MVSYM(ICORE(I0),MAXCOR,IUHF,1,43+ISPIN,143+ISPIN)
               ENDDO
               DO ISPIN=1,IUHF+1
                   CALL UPDMOI(1,NT(ISPIN),ISPIN,190,0,0)
                   CALL UPDMOI(1,NT(ISPIN),ISPIN+2,190,0,0)
               END DO

                CALL GETLST(ICORE(I0), 1, 1, 1, 1, 90)
                CALL PUTLST(ICORE(I0), 1, 1, 1, 1, 190)

                IF (IUHF .NE. 0) THEN
                   CALL GETLST(ICORE(I0), 1, 1, 1, 2, 90)
                   CALL PUTLST(ICORE(I0), 1, 1, 1, 2, 190)
                ENDIF

C
            ENDIF
         ENDIF
C
         WRITE(6,1020)
         RETURN

c     END IF (ICONVG.EQ.0)
      END IF

C DO THE RLE EXTRAPOLATION AND GET THE UPDATED RLE ENERGY.
C
      IF (MOD(IFLAGS(21),2).EQ.0) THEN
         CALL DRRLE(ICORE(I0),MAXCOR,IUHF,RLECYC,.FALSE.)
         CALL DRMOVE(ICORE(I0),MAXCOR,IUHF,0,SING1)
         CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')
      END IF
C
C PRINT AMPLITUDES EVERY IEVERY CYCLES.
C
      IF (bPrint .and. iconvg.ne.0) THEN
         IF (MOD(ICYCLE,IEVERY).EQ.0) THEN
            CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,SING1,'T')
         END IF
      END IF

      icycle=icycle+1
      if (icycle.le.ncycle) goto 555

      WRITE(*,*) '@VCC: The Coupled-Cluster equations did not converge!'
      CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,SING1,'T')
      call aces_exit(1)

c     END IF (.NOT.CC)
      END IF
C
CJDW 10/24/95. Write energy and coordinates in POLYRATE jobs.
      IF (IFLAGS2(113).GT.0) THEN
         IONE=1
         CALL GETREC(20,'JOBARC','NATOMS  ',IONE,NATOMS)
         CALL POLYPRT0(NATOMS,IINTFP,IFLAGS(1),
     &                 .TRUE.,.TRUE.,.FALSE.,.FALSE.)
      END IF
C
C IN THE CASE OF OUT OF CORE ALGORITHM SAVE AMPLITUDE AND INTERMEDIATE
C LISTS IF NECESSARY
      IF (IFLAGS(35).NE.0) THEN
         CALL ACES_AUXCACHE_FLUSH
         CALL ACES_AUXCACHE_RESET
      END IF

 1020 FORMAT(/,77('-'),/,32X,'Exiting xvcc',/,77('-'),/)
      WRITE(6,1020)

      RETURN
      END

