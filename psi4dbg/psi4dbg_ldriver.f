










      SUBROUTINE PSI4DBG_LDRIVER(WORK,MAXCOR,IUHF)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION WORK(MAXCOR)
      DIMENSION ECORR(3)
      INTEGER RELCYC
     
      LOGICAL PCCD,CCD,LCCD
      LOGICAL NONHF
      LOGICAL CIS,EOM
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

      COMMON/METH/PCCD,CCD,LCCD
      COMMON /ENERGY/ ENERGY(500,2)
      COMMON/EXCITE/EOM,CIS
      COMMON/NHFREF/NONHF
      COMMON/EXTRAPO/RLE

C Initializations of various intermediate lists etc.
C
      Write(6,*)
      Write(6,"(2x,a)") "--------Entering the pccd_ldriver-------------"
      INEXT   = 1
      RLECYC  = 0
      IREF    = 1
      IBWPASS = 0
      IRREPX  = 1
      ICONTL  = IFLAGS(4)
      NCYCLE  = IFLAGS(7)
      NKEEP   = IFLAGS(12)
      IEVERY  = IFLAGS(13)

      CALL PCCD_SETLST(WORK(INEXT),MAXCOR,IUHF,"L")
      CALL PCCD_INIT(WORK(INEXT),MAXCOR,IUHF,"L")

      IF (LCCD) THEN
         Write(6,*)
         Write(6,"(a)") "    Lambda is constructed for LCCD."
         CALL MVSYM(WORK(INEXT),MAXCOR,IUHF,IRREPX,46,146)
         CALL PCCD_RNABIJ(WORK(INEXT),MAXCOR,IUHF,'L')
         CALL GETLST(WORK(INEXT),1,1,1,1,90)
         CALL PUTLST(WORK(INEXT),1,1,1,1,190)
         CALL PCCD_LGENINT(WORK(INEXT),MAXCOR,IUHF)
         WRITE(6,1020)
         RETURN
      ENDIF
     
      IF (IFLAGS(21) .EQ.1) THEN
         CALL PCCD_DIISLST(1,IUHF,.FALSE.,"L")
      END IF
      CALL PCCD_FORMWL(WORK(INEXT),MAXCOR,IUHF)
      CALL PCCD_WTWTW(WORK(INEXT),MAXCOR,IUHF)
      CALL PCCD_FORMW4(WORK(INEXT),MAXCOR,IUHF,.TRUE.,.TRUE.,
     &                 .TRUE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.,
     &                 .TRUE.)

      CALL PCCD_FORMW5(WORK(INEXT),MAXCOR,IUHF,.TRUE.,.TRUE.,
     &                 .TRUE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.,
     &                 .TRUE.)
      ICYCLE = 1
 100  CONTINUE

      RELCYC = RLECYC + 1
      CALL PCCD_INITIN(WORK(INEXT),MAXCOR,IUHF)
      CALL PCCD_INITSN(WORK(INEXT),MAXCOR,IUHF,ICYCLE,"L")

      CALL PCCD_CMPENG(WORK(INEXT),MAXCOR,146,ECORR,ENERGY(1,1))

      CALL PCCD_LGENINT(WORK(INEXT),MAXCOR,IUHF)
C
       call pccd_checkintms(work(inext),maxcor,iuhf,0)
       call pccd_check_t2(work(inext),Maxcor,Iuhf)
C
C Evaluate P(ab)Sum_e T(ij,ae){F(b,e)-1/2 sum_m T(m,b)F(m,e)}
C contribution to T2 and Sum_e T(i,e)F(a,e) contribution to T1.
C
      CALL PCCD_F1INL2(WORK(INEXT),MAXCOR,IUHF)
      CALL PCCD_F2INL2(WORK(INEXT),MAXCOR,IUHF)

       call pccd_check_t2(work(inext),Maxcor,Iuhf)

      CALL PCCD_G1INL2(WORK(INEXT),MAXCOR,IUHF)
      CALL PCCD_G2INL2(WORK(INEXT),MAXCOR,IUHF)

       call pccd_check_t2(work(inext),Maxcor,Iuhf)

C   +P(ij)P(ab) Sum_me T(i,e)T(m,a) W(mb,ej) (DRRNG)
C   +1/2 Sum_mn Tau(mn,ab)W(mn,ij) +
C    1/2 Sum_ef Tau(ij,ef)W(ab,ef)           (DRLAD)
C
      CALL PCCD_L2INL2(WORK(INEXT),MAXCOR,IUHF,0)

       call pccd_check_t2(work(inext),Maxcor,Iuhf)

C Denominator weigh the T2 increments. 
C
      IF (PCCD) CALL PCCD_RESET(WORK(INEXT),MAXCOR,IUHF,63)
      CALL PCCD_NEWT2(WORK(INEXT),MAXCOR,IUHF)

       call gwinl1(work(Inext),Maxcor,Iuhf)
       call pccd_l2inl1(work(inext),Maxcor,Iuhf)
       call pccd_check_t2(work(inext),Maxcor,Iuhf)
C
      CALL PCCD_DRTSTS(WORK(INEXT),MAXCOR,ICYCLE,IUHF,ICONVG,ICONTL,
     &            .FALSE.,0,'T')

      IF (IFLAGS(21) .EQ.1)  THEN
         CALL PCCD_DODIIS0(WORK(INEXT),MAXCOR,IUHF,1,ICYCLE,
     &                    ICONVG,ICONTL,.FALSE.,144,61,190,0,90,2,70,
     &                    '     ',DAMP_PARAMETER)
      END IF

      CALL PCCD_CMPENG(WORK(INEXT),MAXCOR,63,ECORR,ENERGY(ICYCLE+1,1))

      IF (ICONVG.NE.0) THEN

         CALL PCCD_DRMOVE(WORK(INEXT),MAXCOR,IUHF,100)
         CALL PCCD_RNABIJ(WORK(INEXT),MAXCOR,IUHF,'L')

      END IF
      IF (ICONVG.EQ.0) THEN

         CALL PCCD_CMPENG(WORK(INEXT),MAXCOR,146,ECORR,ENERGY(1,1))
         CALL PCCD_AMPSUM(WORK(INEXT),MAXCOR,IUHF,0,.FALSE.,'T')

         CALL PCCD_FINISH(ICYCLE+1)
         CALL PCCD_DDMPTGSS(WORK(INEXT), MAXCOR/IINTFP, IUHF, 0,
     &                      'TGUESS  ')
         WRITE(6,1020)
         RETURN
 1020 FORMAT(/,77('-'),/,32X,'Exiting pccd_ldriver',/,77('-'),/)

      ENDIF 
C
      ICYCLE = ICYCLE +  1
      IF (ICYCLE .LE. NCYCLE) GOTO 100
       
      IF (ICONVG .NE. 0) THEN
         CALL PCCD_AMPSUM(WORK(INEXT),MAXCOR,IUHF,0,.FALSE.,'T')
         Write(6,"(2a)") "The pCCD equations did not converge",
     &                     " in alloted number of cyclces."
      ENDIF 
C
      RETURN
      END
