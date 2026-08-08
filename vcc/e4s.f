










      SUBROUTINE E4S(ICORE,MAXCOR,IUHF,SE4)
C
C  THIS ROUTINE IS THE DRIVER FOR THE CONTRIBUTION OF THE DOUBLES
C  TO THE SINGLES.
C  IT CALCULATES IN ADDITION THE FOURTH-ORDER SINGLE ENBERGY                
C  CONTRIBUTION AND CALCULATES THE NEW T1-AMPLITUDES
C
CEND
C
C CODED JUNE/90  JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,POP1,POP2,VRT1,VRT2
      LOGICAL LAMBDA,NONHF,AOBASIS
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,QCISD,UCC
      EQUIVALENCE(IFLAGS(2),METHOD)
      DIMENSION ICORE(MAXCOR)
      COMMON/SWITCH/MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &              QCISD,UCC
      COMMON/SYM/ POP1(8),POP2(8),VRT1(8),VRT2(8),NTAA,NTBB,NF1AA,
     &            NF1BB,NF2AA,NF2BB
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /NHFREF/ NONHF
C
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &DIRPRD(8,8)
      DATA ONE /1.0/
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
C
C  FLAG IF LAMBDA OR CC EQUATIONS ARE SOLVED
C
      WRITE(6,*) '@E4S-IUHF-DEBUG: IUHF=',IUHF
      LAMBDA=.FALSE.
      AOBASIS=IFLAGS(93).EQ.2
C
      NOCCA=NOCCO(1)
      NOCCB=NOCCO(2)
      NVRTA=NVRTO(1)
      NVRTB=NVRTO(2)
      NWAA=NTAA
      NWBB=NTBB
C
      ITAA=MAXCOR+1-NWAA*IINTFP
      MXCOR=MAXCOR-NWAA*IINTFP
      IF(IUHF.EQ.0) THEN
       ITBB=ITAA
      ELSE
       ITBB=ITAA-NWBB*IINTFP
       MXCOR=MXCOR-NWBB*IINTFP
      ENDIF
C
cjp    
      WRITE(6,*) '@E4S-IBWPASS-DEBUG: ibwpass=',ibwpass,' isbwcc=',
     &           isbwcc
      if(ibwpass.lt.2) then
      IF((METHOD.EQ.3.OR.METHOD.EQ.4).AND.(.NOT.AOBASIS)) THEN
       CALL ZERO(ICORE(ITAA),NWAA)
       IF(IUHF.EQ.1) CALL ZERO(ICORE(ITBB),NWBB)
      ELSE IF(METHOD.EQ.3.OR.METHOD.EQ.4.AND.AOBASIS) THEN
       CALL GETLST(ICORE(ITAA),1,1,1,1,90)
       IF(IUHF.EQ.1) CALL GETLST(ICORE(ITBB),1,1,1,2,90)
      ELSE
       CALL GETLST(ICORE(ITAA),1,1,1,3,90)
       IF(IUHF.EQ.1) CALL GETLST(ICORE(ITBB),1,1,1,4,90)
      ENDIF
C
C  WITHIN THE SPIN ADAPTED RHF CODE NO CALL TO T2t1AA2 IS NECESSARY
C
      WRITE(6,*) '@E4S-PRECHECK-DEBUG: IUHF=',IUHF,' METHOD=',METHOD
      IF(IUHF.EQ.1) THEN
       CALL T2T1AA2(ICORE(ITAA),ICORE,MXCOR,POP1,VRT1,1,LAMBDA,0)
C
       CALL T2T1AA2(ICORE(ITBB),ICORE,MXCOR,POP2,VRT2,2,LAMBDA,0)
      ENDIF
      call checksum_vcc_debug('E4S-AA2-A',ICORE(ITAA),NWAA)
      call checksum_vcc_debug('E4S-AA2-B',ICORE(ITBB),NWBB)
C
C
      CALL T2T1AB2(ICORE(ITAA),ICORE,MXCOR,POP1,POP2,VRT1,VRT2,1,
     &             IUHF,LAMBDA,0)
C
      IF(IUHF.EQ.1) THEN
      CALL T2T1AB2(ICORE(ITBB),ICORE,MXCOR,POP2,POP1,VRT2,VRT1,2,
     &             IUHF,LAMBDA,0)
      ENDIF
      call checksum_vcc_debug('E4S-AB2-A',ICORE(ITAA),NWAA)
      call checksum_vcc_debug('E4S-AB2-B',ICORE(ITBB),NWBB)
C
C  WITHIN THE SPIN ADAPTED RHF CODE NO CALL TO T2T1AA1 IS NECESSARY
C
      IF(.NOT.AOBASIS) THEN
       IF(IUHF.EQ.1) THEN
        CALL T2T1AA1(ICORE(ITAA),ICORE,MXCOR,POP1,VRT1,1,LAMBDA,0)
C
        CALL T2T1AA1(ICORE(ITBB),ICORE,MXCOR,POP2,VRT2,2,LAMBDA,0)
       ENDIF
      call checksum_vcc_debug('E4S-AA1-A',ICORE(ITAA),NWAA)
      call checksum_vcc_debug('E4S-AA1-B',ICORE(ITBB),NWBB)
C
C
       CALL T2T1AB1(ICORE(ITAA),ICORE,MXCOR,POP1,POP2,VRT1,VRT2,1,
     &              IUHF,LAMBDA,0)
C
       IF(IUHF.EQ.1) THEN
       CALL T2T1AB1(ICORE(ITBB),ICORE,MXCOR,POP2,POP1,VRT2,VRT1,2,
     &              IUHF,LAMBDA,0)
       ENDIF
      ENDIF
      call checksum_vcc_debug('E4S-AB1-A',ICORE(ITAA),NWAA)
      call checksum_vcc_debug('E4S-AB1-B',ICORE(ITBB),NWBB)
      IF(NONHF)THEN
       CALL GETLST(ICORE,1,1,1,3,93)
       CALL SAXPY(NWAA,ONE,ICORE,1,ICORE(ITAA),1)
       IF(IUHF.NE.0)THEN
        CALL GETLST(ICORE,1,1,1,4,93)
        CALL SAXPY(NWBB,ONE,ICORE,1,ICORE(ITBB),1)
       ENDIF
      ENDIF
C
C     CALCULATE NEW AMPLITUDES AND MBPT(4) SINGLE ENERGY
C
cjp ibwpass1
      endif
cjp
      SE4=0.0D0
      IDEN=1
      IE=IDEN+IINTFP*NWAA
      ID=IE+IINTFP*(NOCCA+NVRTA)
      IEND=ID+IINTFP*NWAA
      IF(IEND.LT.MXCOR) THEN
cjp prepare for pass 2
      if(ibwpass.eq.2) then
cjp read the intermediate result instead of the previous computation
cjp do it transparently for E4SENG 
        call getlst(icore(itaa),1,1,1,3,90)
        if(iuhf.eq.1) call getlst(icore(itbb),1,1,1,4,90)
      endif
C
       CALL E4SENG(ICORE(ITAA),ICORE(ID),ICORE(IDEN),ICORE(IE),
     &             NWAA,POP1,VRT1,NOCCA,NVRTA,SE4,1)
       IF((METHOD.EQ.3.OR.METHOD.EQ.4).AND.(.NOT.AOBASIS)) THEN
        CALL UPDMOI(1,NWAA,1,90,0,0)
        CALL PUTLST(ICORE(ID),1,1,1,1,90)
       ELSE IF(METHOD.EQ.3.OR.METHOD.EQ.4.AND.AOBASIS) THEN
        CALL PUTLST(ICORE(ID),1,1,1,1,90)
       ELSE
cjp here the new t1 amplitudes are written
cjp old t1's are not overwritten!
C
        CALL PUTLST(ICORE(ID),1,1,1,3,90)

       ENDIF
      ENDIF

      IF(IUHF.EQ.0) THEN
       SE4=2.0D0*SE4
      ELSE
       IDEN=1
       IE=IDEN+IINTFP*NWBB
       ID=IE+IINTFP*(NOCCB+NVRTB)
       IEND=ID+IINTFP*NWBB
       IF(IEND.LT.MXCOR) THEN
         call checksum("T1BB-increment", ICORE(ITBB), NWBB)
C
        CALL E4SENG(ICORE(ITBB),ICORE(ID),ICORE(IDEN),ICORE(IE),
     &              NWBB,POP2,VRT2,NOCCB,NVRTB,SE4,2)
        IF((METHOD.EQ.3.OR.METHOD.EQ.4).AND.(.NOT.AOBASIS)) THEN
         CALL UPDMOI(1,NWBB,2,90,0,0)
         CALL PUTLST(ICORE(ID),1,1,1,2,90)
        ELSE IF(METHOD.EQ.3.OR.METHOD.EQ.4.AND.AOBASIS) THEN
         CALL PUTLST(ICORE(ID),1,1,1,2,90)
        ELSE
         CALL PUTLST(ICORE(ID),1,1,1,4,90)
        ENDIF
       ENDIF
      ENDIF
      SING1=.TRUE.
         call Check_T2_VCC(Icore,Maxcor,Iuhf)
      RETURN
      END
