






































































































































































































      SUBROUTINE FORMT2(W,DENOM,IUHF)
C
C THIS ROUTINE FORMS THE FIRST ORDER T2 VECTOR IN SYMMETRY 
C  PACKED FORM AND WRITES IT OUT.  ALSO COMPUTE THE CONTRIBUTION
C  TO THE CORRELATION ENERGY.
C
C      ARRAYS : W - USED TO HOLD INTEGRAL ARRAYS.  MUST
C                   BE DIMENSIONED TO LARGEST IRREP OF THIS
C                   ARRAY.
C
C               DENOM - USED FOR DENOMINATOR ARRAYS.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION W(1),DENOM(1),EAA,EAB,ENERGY,SDOT,ETAA,ESCF
      DOUBLE PRECISION ESING
      LOGICAL NONHF
      LOGICAL POLSCF
      CHARACTER*2 SPCASE(3)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFEA(2),NFMI(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),D(18)
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /SHIFT/  ISHIFT,NDRGEO
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


      DATA SPCASE /'AA','BB','AB'/
      ETAA=0.0
cjp avoid printing of senseless numbers in the case of bwcc
      if(isbwcc) write(6,*)
     &  '@FORMT2-I BW-CC: pseudo MBPT(2) energies not printed'
      if(.not.isbwcc) WRITE(6,1010)
1010  FORMAT(/,T3,'@FORMT2-I, Second-order MBPT correlation energies:')
cjp
      NONHF=(IFLAGS(38)+IFLAGS(77).NE.0).or.isbwcc
      IF(NONHF)THEN
       if(.not.isbwcc) WRITE(6,1011)
1011   FORMAT(T3,'@FORMT2-I, Singles contribution will be calculated.')
       IF(IFLAGS(39).EQ.0)THEN
        if(.not.isbwcc) WRITE(6,1012)
1012    FORMAT(T3,'@FORMT2-W, MBPT(2) energies are correct only ',
     &           'for semicanonical orbitals.')
       ENDIF
      ENDIF 
      CALL GETREC(20,'JOBARC','SCFENEG ',IINTFP,ESCF)
      if(.not.isbwcc) WRITE(6,1005)
      if(.not.isbwcc) WRITE(6,1004)ESCF
cjp here initial t value is computed as integrals.1/denominators
cjp listd .. .denoms; listt ... amplitudes; listw ... integrals
cjp and they are immediatelly used to comp. initial correlation energy
      DO 5 ISPIN=1,1+IUHF
       LISTW=13+ISPIN + ISHIFT
       LISTD=47+ISPIN
       LISTT=43+ISPIN
       EAA=0.0 
       DO 10 IRREP=1,NIRREP
        NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
        DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))
        Print*, NUMDIS,DISSIZ
        CALL GETLST(W,1,NUMDIS,2,IRREP,LISTW)
        CALL GETLST(DENOM,1,NUMDIS,1,IRREP,LISTD)
        DO 20 IELEM=1,DISSIZ*NUMDIS
         DENOM(IELEM)=DENOM(IELEM)*W(IELEM)
20      CONTINUE
        CALL PUTLST(DENOM,1,NUMDIS,1,IRREP,LISTT)
        CALL CHECKSUM("@-INTPRC:T2", DENOM, NUMDIS*DISSIZ)
        EAA=EAA+SDOT(NUMDIS*DISSIZ,W,1,DENOM,1)
10     CONTINUE
       if(.not.isbwcc) WRITE(LUOUT,1000)SPCASE(ISPIN),EAA
       IF(IUHF.EQ.0)EAA=2.0*EAA
       ETAA=ETAA+EAA
5     CONTINUE
      LISTW=16 + ISHIFT
      LISTD=50
      LISTT=46
      EAB=0.0
      DO 110 IRREP=1,NIRREP
       NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
       DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))
cjp order of integrals and denominators same, since the integrals were sorted
cjp wrto the true fermi vacuum and the denoms were constructed so that
cjp the order corresponds to the same f vacuum
cjp check it anyway in the dijab and in the integral symmetry packing routine
       CALL GETLST(W,1,NUMDIS,2,IRREP,LISTW)
       CALL GETLST(DENOM,1,NUMDIS,1,IRREP,LISTD)
       DO 120 IELEM=1,DISSIZ*NUMDIS
        DENOM(IELEM)=DENOM(IELEM)*W(IELEM)
120    CONTINUE
       CALL PUTLST(DENOM,1,NUMDIS,1,IRREP,LISTT)
       write(6,"(a)") " The T2-ABAB"
CSSS       call output(Denom,1,DISSIZ,1,NUMDIS,DISSIZ,NUMDIS,1)
       CALL CHECKSUM("@-INTPRC:T2-ABAB", DENOM, NUMDIS*DISSIZ)
       EAB=EAB+SDOT(NUMDIS*DISSIZ,W,1,DENOM,1)
110    CONTINUE
       if(.not.isbwcc) WRITE(LUOUT,1000)SPCASE(3),EAB
       ENERGY=EAB+ETAA
C
C CALCULATE SINGLES CONTRIBUTION TO SECOND-ORDER ENERGY
C  FOR SEMICANONICAL ORBITALS
C
       ESING=0.0
       IF(NONHF)THEN
        DO 130 ISPIN=1,1+IUHF
         I000=1
         I010=I000+NT(ISPIN)
         I020=I010+NT(ISPIN)
         CALL GETLST(W(I000),1,1,1,ISPIN,90)
         CALL GETLST(W(I010),1,1,1,ISPIN+2,93)
         ESING=ESING+SDOT(NT(ISPIN),W(I000),1,W(I010),1)
         If ((Iflags(2) .EQ. 45) .OR. 
     &       (Iflags(2) .EQ. 46)) 
     &       ESING = 2.0D0*ESING 
   
130     CONTINUE
        if(.not.isbwcc) WRITE(LUOUT,1003)ESING
       ENDIF
C
       if(.not.isbwcc) WRITE(LUOUT,1002)ENERGY
       if(.not.isbwcc) WRITE(LUOUT,1001)ENERGY+ESCF+ESING
C
C WRITE TOTAL ENERGY TO JOBARC (THIS IS AT LEAST NEEDED FOR
C OPTIMIZATIONS USING NUMERICALLY EVALUATED GRADIENTS
C
CJDW   8/30/95
C
C      Do not write out the second-order energy in scf analytical hessian
C      polyrate runs.
C Not just for polyrate, but forall analytic SCF Hessian calcualtions
C 12/09, Ajith Perera
C
       POLSCF = .FALSE.
       POLSCF = IFLAGS(2).EQ.0 .AND. IFLAGS(3).EQ.2 
cjp
cjp do this only for reference no. 1; 
       IF(POLSCF.or.isbwcc.and.nref.gt.1)THEN
        if(.not.isbwcc) WRITE(6,1020)
       ELSE
        CALL PUTREC(20,'JOBARC','TOTENERG',IINTFP,ENERGY+ESCF+ESING)
       ENDIF
cjp       if(.not.isbwcc) print *,'in formt2 energies', ENERGY+ESCF+ESING,ENERGY,ESCF,ESING
C
       WRITE(6,1005)
1000   FORMAT(T15,'E2(',A2,')',T37,'=',F18.12,' a.u.')
1001   FORMAT(T15,'Total MBPT(2) energy',T37,'=',F18.12,' a.u.')
1002   FORMAT(T15,'E2(TOT)',T37,'=',F18.12,' a.u.')
1003   FORMAT(T15,'E2(SINGLE) ',T37,'=',F18.12,' a.u.')
1004   FORMAT(T15,'E(SCF)',T37,'=',F18.12,' a.u.')
1005   FORMAT(T15,46('-'))
cjp
1020   FORMAT(T15,' @FORMT2-I, MR-BWCC or SCF Hessian Polyrate Run.',
     &            ' MBPT(2) energy not written. ')
       RETURN
       END
