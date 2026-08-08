










      SUBROUTINE PCCD_E4SENG(W,D,DEN,E,NTSIZ,POP,VRT,NOCC,NVRT,
     +                       SE4,ISPIN)
C
C     THIS ROUTINE CALCULATES THE NEW T1 AMPLITUDES AND
C     FOR FOURTH ORDER PERTURBATION THEORY THE SINGLE
C     CONTRIBUTION TO THE ENERGY
C
C     THIS ROUTINE HAS TO BE CALLED SEPARATELY FOR ALPHA AND BETA
C     CASE
cjp
cjp in bwcc, this routine is called twice. in the first pass it
cjp collects matrix elements of Heff, which correspond to 
cjp single excitations (uhf only), in the second pass, called after
cjp Heff diagonalalization is complete, the amplitudes are finally updated
cjp
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINCC,CICALC
      INTEGER DIRPRD,POP,VRT
cjp
      integer itrue,atrue,intptr
      DIMENSION W(NTSIZ),E(NOCC+NVRT),
     &          D(NTSIZ),DEN(NTSIZ),POP(8),VRT(8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /LINEAR/ LINCC,CICALC
      COMMON /CORENG/ ELAST
C
      EQUIVALENCE(IFLAGS(2),METHOD)
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


C
C     GET ORBITAL ENERGIES FROM DISK
C
      IF(ISPIN.EQ.1)THEN
       CALL GETREC(20,'JOBARC','SCFEVALA',IINTFP*(NOCC + NVRT),
     &             E)
      ELSE
       CALL GETREC(20,'JOBARC','SCFEVALB',IINTFP*(NOCC + NVRT),
     &             E)
      ENDIF
cjp fix diagonal fock elements
      if(isbwcc) then
        do i=1,nocc+nvrt
          e(i)=e(i)+fockcd(i,iref,ispin)
        enddo
      endif
C
C     LOOP OVER IRREPS 
C
      IND=0
      INDI=0
      INDA0=0
C
      DO 100 IRREP=1,NIRREP
C
      NOCCI=POP(IRREP)
      NVRTI=VRT(IRREP)
c
c
C
cjp
cjp we do this both independently on ibwpass, although it might be splitted 
cjp in the 2 passes
cjp identify also the internal amplitudes
cjp and check them if they come out zero;
cjp set the denominator to blowdenom
cjp the following is also completely wrong if symmetry is on
cjp NOTE for parallelization ... parts of these arrays are meaningless
      intptr=1
      DO  10 I=1,NOCCI
      INDI=INDI+1
      INDA=INDA0
      DO  10 IA=1,NVRTI
      INDA=INDA+1
      IND=IND+1
      if(isbwcc) then
           itrue=invhnum(indi,iref,ispin)
           atrue=invpnum(inda,iref,ispin)
           if(isactive(itrue,ispin).and.isactive(atrue,ispin))then
             den(ind)=denomblow
cjp
cjp collect matrix elements of effective hamiltonian
cjp
          if(ibwpass.eq.1.and.nref.gt.1) then
            if(intptr.le.internnum1(iref,ispin) .and.
     +      internindex1(intptr,iref,ispin).eq.ind) then
cjp             this is offdiag element, but maybe such one which does not belong to current iref
                if(iref.eq.internfrom1(intptr,iref,ispin)) then
cjp                     relevant one
                        heff(internto1(intptr,iref,ispin),iref)=w(ind)
                endif
            intptr=intptr+1
cjp         they are sorted acc. to increasing internindex, look simply at next
           else
c
cjp if such internal excitation is not included in the model space
cjp notice the user about it
             if(abs(w(ind)).gt.1d-8) then
        write(6,*)'@E4SENG-W: ref. ',iref,
     &  ' nonzero internal T1 amplitude:',
     &               itrue,atrue,w(ind),e(itrue)-e(atrue),w(ind)
     &    /(e(itrue)-e(atrue)),
     & '  ADDITIONAL REFERENCES MIGHT BE IMPORTANT'

             bwwarning(2)=.true.
             bwwarntext(2)='@E4SENG-W: Nomnzero internal T1 amplitudes
     & encountered, additional references might be important'
        if(useeq429) stop 'algorithm to treat  nonzero 
     & internal T1s not activated'
        endif
              endif
cjp             intptr.le.internnum1(iref,ispin) .and. ...
            endif
cjp             ibwpass.eq.1.and.nref.gt.1
cjp end collect matrix elements of effective hamiltonian
cjp
cjp    was not internal excitation
           else
             den(ind)=e(itrue)-e(atrue)
           endif
      else
      DEN(IND) = E(INDI) - E(NOCC + INDA)
      endif
   10 CONTINUE
      INDA0=INDA0+NVRTI
100   CONTINUE
C
C     CALCULATE NEW AMPLITUDES
C
cjp do not do it if ibwpass ==1
        if(ibwpass.eq.1) then
cjp store the intermediates as if they were final amplitudes
cjp they will be needed anyway in the general term of eq 4.28 in Hubac paper
        do i=1,ntsiz
        d(i)=w(I)
        enddo
        return
        endif
      FACT=0.0
      IF(CICALC)FACT=ELAST
cjp BW-CC shift of the denominator
      if(isbwcc) fact=ecorrbw
      DO 20 I=1,NTSIZ
C
      tmp=1.d0/(DEN(I)+FACT)
cjp intruder prevention
      if(isbwcc .and. abs(tmp).gt.intruder) tmp=0.d0
      D(I)=W(I)*tmp
20    CONTINUE

      SE4=SE4+SDOT(NTSIZ,D,1,W,1)

      RETURN
      END
