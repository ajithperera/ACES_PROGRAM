



























































































































































































      SUBROUTINE PCCD_ULDRIVER(ICORE,MAXCOR,IUHF)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION ICORE(MAXCOR)
      DIMENSION ECORR(3)
      INTEGER RELCYC

      LOGICAL NONHF
      LOGICAL CIS,EOM
      LOGICAL DO_TAU
      LOGICAL PCCD,CCD,LCCD
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

      COMMON/CALC/PCCD,CCD,LCCD
      COMMON /ENERGY/ ENERGY(500,2)
      COMMON/EXCITE/EOM,CIS
      COMMON/NHFREF/NONHF
      COMMON/EXTRAPO/RLE

C Initializations of various intermediate lists etc.

      I0      = 1 
      INEXT   = 1
      RLECYC  = 0
      IREF    = 1
      IBWPASS = 0
      IRREPX  = 1
      ICONTL  = IFLAGS(4)
      NCYCLE  = IFLAGS(7)
      NKEEP   = IFLAGS(12)
      IEVERY  = IFLAGS(13)
      IPRTENG = 1

C-----------------------------------------------------------------------

      CALL PCCD_INILAM(IUHF)
      CALL ZERO(ENERGY,1000)

      IF (LCCD) THEN
         Write(6,"(a)") "      Lambda is constructed for LCCD."
         CALL MVSYM(ICORE(I0),MAXCOR,IUHF,IRREPX,46,146)
         CALL GETLST(ICORE(I0),1,1,1,1,90)
         CALL PUTLST(ICORE(I0),1,1,1,1,190)
         CALL LGENINT(ICORE(I0),MAXCOR,IUHF)
         WRITE(6,1020)
         RETURN
      ENDIF

      WRITE(6,801)
      CALL AMPSUM(ICORE(i0),MAXCOR,IUHF,100,SING,'L')
801   FORMAT('  Initial lambda amplitudes: ')
      WRITE(6,800)
  800 FORMAT(T3,' Correlation energies computed from initial ',
     &          'L amplitudes: ')
      CALL CMPENG2(ICORE(i0),MAXCOR,143,0,ECORR,ENERGY(1,1),
     &              ENERGY(1,2),IUHF,IPRTENG)
C
      IF (IFLAGS(21) .EQ.1) THEN
         CALL PCCD_DIISLST(IRREPX,IUHF,.FALSE.,"L")
      ENDIF 
C
      CALL RNABIJ(ICORE(i0),MAXCOR,IUHF,'L')
      CALL STLLST(ICORE(i0),MAXCOR,IUHF)
C
C CALCULATE FOR CCSD THE W(AB,EF) INTERMEDIATES
C
C
      CALL FORMWL(ICORE(i0),MAXCOR,IUHF) 
      CALL WTWTW(ICORE(i0),MAXCOR,IUHF)

      TERM1 = .TRUE.
      TERM2 = .TRUE.
      TERM3 = .TRUE.
      TERM4 = .TRUE.
      TERM5 = .TRUE.
      TERM6 = .TRUE.

      DO_TAU = .TRUE.
      CALL FORMW4(ICORE(i0),MAXCOR,IUHF,TERM1,TERM2,TERM3,
     &            TERM4,TERM5,TERM6,DO_TAU,.TRUE.)

      TERM1=.TRUE.
      TERM2=.TRUE.
      TERM3=.TRUE.
      TERM4=.TRUE.
      TERM5=.TRUE.
      TERM6=.TRUE.
      CALL FORMW5(ICORE(i0),MAXCOR,IUHF,TERM1,TERM2,TERM3,
     &            TERM4,TERM5,TERM6,DO_TAU,.TRUE.)

c************************************************************
c from now on lists 33-39 (newly created taus) are on 233-239
c*************************************************************

      ICYCLE = 1
 100  CONTINUE

      RLECYC=RLECYC+1
C
C FILL LAMBDA(2) INCREMENTS WITH <IJ||AB> INTEGRALS
C
C AND ZERO THE LAMBDA(1) INCREMENTS (ONLY CCSD AND QCISD)
C
      CALL INITIN(ICORE(i0),MAXCOR,IUHF)
CSSS      CALL INITSN(ICORE(i0),MAXCOR,IUHF)
C
C GENERATE V AND G INTERMEDIATES.
C
      call check_leom(icore(i0),Maxcor,Iuhf)

      CALL LGENINT(ICORE(i0),MAXCOR,IUHF)

      call check_leom(icore(i0),Maxcor,Iuhf)
C
C COMPUTE F AND G INTERMEDIATE CONTRIBUTION TO LAMBDA(2) INCREMENT.
C
      CALL F1INL2(ICORE(i0),MAXCOR,IUHF)
      CALL F2INL2(ICORE(i0),MAXCOR,IUHF)

      Write(6,*) "The lambda residuals before GinL2"
      call check_leom(icore(i0),Maxcor,Iuhf)

      CALL PCCD_G1INL2U(ICORE(i0),MAXCOR,IUHF)
      CALL PCCD_G2INL2U(ICORE(i0),MAXCOR,IUHF)

      Write(6,*) "The lambda residuals after GinL2"
      call check_leom(icore(i0),Maxcor,Iuhf)
C
C DO W INTERMEDIATE CONTRIBUTION TO T2 EQUATION (SAME AS THIRD ORDER CODE).
C
      CALL PCCD_L2INL2U(ICORE(i0),MAXCOR,IUHF,0)
C
      Write(6,*) "The lambda residuals after l2lad"
      call check_leom(icore(i0),Maxcor,Iuhf)
      IF (pCCD)  THEN
C Ajith Perera, 10/2021.
         Write(6,"(20x,a)") " ----------Warning------------"
         Write(6,"(3a)")  "  A pCCD or pCCD like calculation is being",
     +                   " perforemed and the off-diagonal"
         write(6,"(a)")   " blocks of T2 is set to zero"
         Write(6,"(20x,a)") " -----------------------------"
         call pccd_reset_vcc(Icore(i0),Maxcor,Iuhf,63)
      ENDIF 
      Write(6,*) "The lambda residuals after l2inl2"
      call check_leom(icore(i0),Maxcor,Iuhf)
C
C DENOMINATOR WEIGHT T2 INCREMENTS TO FORM NEW T2.
C
      CALL NEWT2(ICORE(i0),MAXCOR,IUHF)

      Write(6,*) "The lambda residuals after newT2"
      call check_leom(icore(i0),Maxcor,Iuhf)
      CALL DRTSTS(ICORE(i0),MAXCOR,ICYCLE,IUHF,ICONVG,ICONTL,
     &              SING,100,'L')
      CALL CMPENG2(ICORE(i0),MAXCOR,60,2,ECORR,ENERGY(ICYCLE+1,1),
     &               ENERGY(ICYCLE+1,2),IUHF,IPRTENG)

      IF (IFLAGS(21) .EQ.1) THEN
        CALL PCCD_DODIIS0(ICORE(i0),MAXCOR/IINTFP,IUHF,1,ICYCLE,
     &                    ICONVG,ICONTL,.FALSE.,144,61,190,0,90,
     &                    2,70)
      ENDIF

      IF (ICONVG .NE. 0)  THEN
         CALL DRMOVE(ICORE(i0),MAXCOR,IUHF,100,SING)
         CALL RNABIJ(ICORE(i0),MAXCOR,IUHF,'L')
      ENDIF 

      IF(ICONVG.EQ.0) THEN

         CALL CMPENG2(ICORE(i0),MAXCOR,60,2,ECORR,ENERGY(ICYCLE+1,1),
     &                ENERGY(ICYCLE+1,2),IUHF,IPRTENG)
         CALL AMPSUM(ICORE(i0),MAXCOR,IUHF,0,.FALSE.,'L')

CSSS         CALL FINISH(ICYCLE+1,PCCD,.FALSE.,.FALSE.)
         CALL PCCD_UFINISH(ICYCLE+1,"L")
         CALL PCCD_DDMPTGSS(ICORE(i0), MAXCOR/IINTFP, IUHF, 0,
     &                      'TGUESS  ')
         WRITE(6,1020)
 1020    FORMAT(/,77('-'),/,32X,'Exiting pccd_ldriver',/,77('-'),/)
         RETURN
      ENDIF

      ICYCLE = ICYCLE +  1
      IF (ICYCLE .LE. NCYCLE) GOTO 100

      IF (ICONVG .NE. 0) THEN
         CALL AMPSUM(ICORE(i0),MAXCOR,IUHF,0,.FALSE.,'T')
         Write(6,"(2a)") "The pCCD equations did not converge",
     &                     " in alloted number of cyclces."
      ENDIF

      RETURN
      END
