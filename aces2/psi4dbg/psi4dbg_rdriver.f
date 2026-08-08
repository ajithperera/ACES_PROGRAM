










      SUBROUTINE PSI4DBG_RDRIVER(IOCRE,MAXCOR,IUHF,MICRO_ITER)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION ICORE

      DIMENSION ICORE(MAXCOR)
      DIMENSION ECORR(3)
      INTEGER RELCYC
     
      LOGICAL PCCD,CCD,LCCD,MBPT2
      LOGICAL PCCDS,PCCDTS
      LOGICAL NONHF
      LOGICAL CIS,EOM
      LOGICAL BRUECK,BRKCNV
      LOGICAL MICRO_ITER
      LOGICAL READT
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
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end


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


c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end

      COMMON/METH/PCCD,CCD,LCCD,MBPT2
      COMMON /ENERGY/ ENERGY(500,2)
      COMMON/EXCITE/EOM,CIS
      COMMON/NHFREF/NONHF
      COMMON/EXTRAPO/RLE

C Note at the first-iteration, T1 vector is zero and T2 is MBPT(2).

      CALL PSI4DBG_SETMET()
C
C Initializations of various intermediate lists etc.
C
      INEXT   = 1
      RLECYC  = 0
      IREF    = 1
      IBWPASS = 0
      I0      = 1
      BRUECK  = .FALSE.
      ICONTL  = IFLAGS(4)
      NCYCLE  = IFLAGS(7)
      NKEEP   = IFLAGS(12)
      IEVERY  = IFLAGS(13)
      BRUECK  = (IFLAGS(22) .EQ. 1)
      MBPT    = IFLAGS(2) .EQ. 1
      pCCD    = IFLAGS(2) .EQ. 52
      pCCDS   = IFLAGS(2) .EQ. 55
      pCCDTS  = IFLAGS(2) .EQ. 56
      If (pCCD) THEN
         METHOD=8
      ELSEIF (pCCDS) THEN
         METHOD=10
      ELSEIF (pCCDTS) THEN
         METHOD=10
      ELSEIF (MBPT2) THEN
         METHOD = 1
      ENDIF 

      IF (METHOD.GE.1.OR.NONHF) CALL INSING(METHOD,IUHF)

      CALL CMPENG(ICORE(I0),MAXCOR,43,0,ECORR,
     &            ENERGY(1,1),ENERGY(1,2),IUHF,1)
      ELAST = ENERGY(1,2)

      IF ((IFLAGS2(117).EQ. 2  .OR. IFLAGS2(117).EQ. 4 .OR.
     &        IFLAGS2(117).EQ. 5  .OR. IFLAGS2(117).EQ. 10)
     &       .AND. .NOT. NONHF) THEN

      CALL INMBPT(ICORE(I0),MAXCOR,IUHF)
      CALL INSING(10,IUHF)
      CALL SETLST(ICORE(I0),MAXCOR,IUHF)
      CALL INITIN(ICORE(I0),MAXCOR,IUHF)
      CALL INITSN(ICORE(I0),MAXCOR,IUHF)
      CALL INCOR(I0,MAXCOR,IUHF)
      CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')
      CALL GENINT(ICORE(I0),MAXCOR,IUHF,2,.FALSE.)
      CALL INIT2(IUHF,PCCD.OR.PCCDS.OR.PCCDTS)
C
C  PUT CONVERGED AMPLITUDES AS INITIAL GUESS FOR SOLVING THE
C  LAMBDA EQUATIONS ON THE GAMLAM FILE
C
      DO ISPIN = 3, 3-2*IUHF, -1
         CALL MVSYM(ICORE(I0),MAXCOR,IUHF,1,43+ISPIN,143+ISPIN)
      END DO
      INTL1 = .FALSE.
      IF (INTL1) THEN
         DO ISPIN=1,IUHF+1
            CALL UPDMOI(1,NT(ISPIN),ISPIN,190,0,0)
            CALL UPDMOI(1,NT(ISPIN),ISPIN+2,190,0,0)
         END DO
C ZERO THE L1 LIST (ALL CYCLES)
         CALL ZERLST(ICORE,NT(1),1,1,1,190)
         IF (IUHF.NE.0) CALL ZERLST(ICORE,NT(2),1,1,2,190)
      ENDIF 
C
      IF (IUHF.EQ.0) CALL RESET(ICORE(I0),MAXCOR,IUHF)
      IF (IFLAGS(35).NE.0) THEN
          CALL ACES_AUXCACHE_FLUSH
          CALL ACES_AUXCACHE_RESET
      END IF
      WRITE(6,1020)
      ENDIF
C
      CALL PCCD_INITL2D2(IUHF)
      IF (MICRO_ITER) THEN 
         Write(6,*)
         CALL INMBPT(ICORE(I0),MAXCOR,IUHF)
         CALL PSI4DBG_DRDTGSS(ICORE(I0),MAXCOR,IUHF,'TGUESS  ',0,READT)
CSSS         CALL PSI4DBG_RETRIVE(ICORE(I0),MAXCOR,IUHF)
         RETURN
      ENDIF 
      CALL PCCD_RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')
      IF (NONHF) THEN
         CALL PSI4DBG_ROHFPT(ICORE(I0),MAXCOR,IUHF)
         CALL INMBPT(ICORE(I0),MAXCOR,IUHF)
         CALL INSING(10,IUHF)
         CALL SETLST(ICORE(I0),MAXCOR,IUHF)
         CALL INITIN(ICORE(I0),MAXCOR,IUHF)
         CALL INITSN(ICORE(I0),MAXCOR,IUHF)
         CALL INCOR(I0,MAXCOR,IUHF)
CSSS         CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')
CSSS         CALL GENINT(ICORE(I0),MAXCOR,IUHF,2,.FALSE.)
         CALL INIT2(IUHF,PCCD.OR.PCCDS.OR.PCCDTS)
C
C PUT CONVERGED AMPLITUDES AS INITIAL GUESS FOR SOLVING THE
C LAMBDA EQUATIONS ON THE GAMLAM FILE
C
         DO ISPIN = 3, 3-2*IUHF, -1
            CALL MVSYM(ICORE(I0),MAXCOR,IUHF,1,43+ISPIN,143+ISPIN)
         END DO
         INTL1 = .FALSE.
         IF (INTL1) THEN
            DO ISPIN=1,IUHF+1
               CALL UPDMOI(1,NT(ISPIN),ISPIN,190,0,0)
               CALL UPDMOI(1,NT(ISPIN),ISPIN+2,190,0,0)
            END DO
C ZERO THE L1 LIST (ALL CYCLES)
            CALL ZERLST(ICORE,NT(1),1,1,1,190)
            IF (IUHF.NE.0) CALL ZERLST(ICORE,NT(2),1,1,2,190)
         ENDIF 
C
         IF (IUHF.EQ.0) CALL RESET(ICORE(I0),MAXCOR,IUHF)
         IF (IFLAGS(35).NE.0) THEN
             CALL ACES_AUXCACHE_FLUSH
             CALL ACES_AUXCACHE_RESET
         END IF
CSSS         CALL PSI4DBG_STORE(ICORE(I0),MAXCOR,IUHF)
         CALL PSI4DBG_DDMPTGSS(ICORE(I0),MAXCOR,IUHF,0,'TGUESS  ')
         CALL PSI4DBG_MBPTOUT
         WRITE(6,1020)
      ELSE
         CALL DDMPTGSS(ICORE(I0),MAXCOR,IUHF,0,'TGUESS  ')
         CALL PSI4DBG_MBPTOUT
      ENDIF 

 1020 FORMAT(/,77('-'),/,32X,'Exiting xvcc',/,77('-'),/)
C
      RETURN
      END
