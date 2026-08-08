











      SUBROUTINE NEWT2(ICORE,MAXCOR,IUHF)
C
C THIS ROUTINE PICKS UP THE FINAL T2 OR T2 INCREMENTS, DENOMINATOR
C  WEIGHTS THEM AND THEN OVERWRITES THE T2 INCREMENT LIST WITH THE NEW VALUES.
C  GOES ONE IRREP AT A TIME OVER ALL SPIN CASES.  FOR RHF CASES,
C  IT ALSO CALLS A ROUTINE WHICH FORMS THE AA AMPLITUDES FROM THE 
C  AB VALUES.
cjp
c This routine is modified for mr-bw-CC in following way:
c
c 1. ecorrbw variable (contains epsilon0-Heff(mm)) is additional
c     parameter of subroutine
c 2. offdiagonal Heff elements are identified as the internal amplitudes
c 3. 1/D(ij,ab) is replaced by 1/(D(ij,ab)+ecorrbw): instead of vecprd
c     routine a new routine bwvecprd is called
c 4. terms from disconnected diagrams arred added to the RHS of bwcc T2 
c
c This task is split in two passes - the routine will be called twice
c In the first pass the offdiagonal Heff elements are gathered
c The second pass must be dine after Heff diagonalization is finished!
c since until then ecorrbw is not known.
c For standard cc, the routine performs everything in one pass
c
cjp

C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ELAST,FACT
      DIMENSION ICORE(MAXCOR)
      LOGICAL LINCC,CICALC
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/ POP1(8),POP2(8),VRT1(8),VRT2(8),NTAA,NTBB,NF1AA,
     &             NF1BB,NF2AA,NF2BB
      COMMON /LINEAR/ LINCC,CICALC
cjp needed common switch
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,QCISD
      LOGICAL NONHF,TRIPIT,TRIPNI,TRIPNI1,T3STOR,PRESNT,INIT,
     &        DORESET,UCC,RESTART,BRUECK
      COMMON /SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC
      COMMON /CORENG/ ELAST
cjp
      COMMON /INFO/ NOCCO(2),NVRTO(2)
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


        integer intptr
         real*8 getheff
         real*8 scalfactor
cjp help variables to pass uniquelly through all 3 spin cases
        integer ntxx,nocvrxx,ispinxx
        dimension ntxx(3),nocvrxx(3),ispinxx(3)
        
cjp
        nocca=nocco(1)
      noccb=nocco(2)
      nvrta=nvrto(1)
      nvrtb=nvrto(2)
      nocvrxx(1)=nocca+nvrta
      nocvrxx(2)=nocca+nvrta
      nocvrxx(3)=nocca+nvrta
      ntxx(1)=ntaa
      ntxx(2)=ntbb
      ntxx(3)=ntaa
      ispinxx(1)=1
      ispinxx(2)=2
      ispinxx(3)=1

      if(isbwcc.and.bwgossip) 
     +    write(6,*)'entering bwnewt2, ecorrbw= ',ecorrbw
      IBOT=1
      IF(IUHF.EQ.0)IBOT=3
      I000=1
      DO 5 ISPIN=IBOT,3
       LSTDEN=47+ISPIN
       LSTINC=60+ISPIN
       DO 10 IRREP=1,NIRREP
        NDSSYM=IRPDPD(IRREP,ISYTYP(2,LSTINC))
        DISSYM=IRPDPD(IRREP,ISYTYP(1,LSTINC))
        NSIZE=NDSSYM*DISSYM
        I010=I000+NSIZE*IINTFP
        I020=I010+NSIZE*IINTFP
        IF(I020.GT.MAXCOR)CALL INSMEM('NEWT2',I020,MAXCOR)
cjp read amplitude increments
        CALL GETLST(ICORE,1,NDSSYM,1,IRREP,LSTINC)
cjp
cjp identify and store Heff offdiagonal elements
cjp skip this if ibwpass ne 1
        if(isbwcc .and. nref.gt.1 .and. ibwpass.eq.1) then
cjp FOR LATER NOTE THAT WE ARE HERE INSIDE THE ISPIN AND IRREP LOOPS!!!
cjp
cjp@@@ NOTE: we should also check if some internal T2 amplitudes
cjp are nonzero, as it is done in e4seng for t1. It might be done
cjp also in bwvecprd - but would probably slow down the code
cjp
        intptr=1
        do i=1,nsize
        if(intptr.le.internnum(iref,ispin) .and. 
     +      internindex(intptr,iref,ispin).eq.i) then
cjp             this is offdiag element, but maybe such one which does not belong to current iref
                if(iref.eq.internfrom(intptr,iref,ispin)) then
cjp                     relevant one
                        heff(internto(intptr,iref,ispin),iref)=
     +                   getheff(icore(I000),i)
                endif
                intptr=intptr+1
cjp             they are sorted acc. to increasing internindex, look simply at next
        endif
        enddo
        endif
cjp end heff gathering
cjp we are done with first pass in bwcc - continue irrep and spin loops
        if(ibwpass.eq.1) goto 10
c
c
c
c
cjp
cjp add unlinked contributions according to Hubac eq. 4.28
cjp in the size-ext correction, these unlinked contributions will be omitted
cjp
        if(isbwcc.and.sing1) then
cjp we need core for temp storage of t1
        i030=i020+ntxx(ispin)*iintfp
        if(ispin.eq.3 .and. iuhf.ne.0) then
        i025=i030
        i030=i025+ntbb*iintfp
        else    
        i025=i020
        endif
        if(i030.gt.maxcor)call insmem('newt2',i030,maxcor)
cjp read old t1 amplitudes (aa case only for rhf)
        call getlst(icore(i020),1,1,1,ispinxx(ispin),90)
        if(ispin.eq.3 .and. iuhf.ne.0)
     +     call getlst(icore(i025),1,1,1,2,90)
c
c
c
cjp the use of useeq429=true is already obsolete
        if(useeq429.and.iuhf.eq.0) then
cjp simplified implementation for comparison with masik's' program
cjp add the t(i,a)t(j,b)-t(j,a)t(i,b) term scaled by ecorrbw
        call ftau(icore(I000),icore(i020),icore(i020),dissym,ndssym,
     +  pop1,pop2,vrt1,vrt2,irrep,ispin,ecorrbw)
        else
c
c
cjp add 2 terms of RHS of 4.28 with single excitation contributions
cjp first some allocation
        i040=i030+ntxx(ispin)*iintfp
        if(ispin.eq.3 .and. iuhf.ne.0) then
        i035=i040
        i040=i035+ntbb*iintfp
        else
        i035=i030
        endif
        i050=i040+ntxx(ispin)*iintfp
        if(ispin.eq.3 .and. iuhf.ne.0) then
        i045=i050
        i050=i045+ntbb*iintfp
        else
        i045=i040
        endif
        i060=i050+nocvrxx(ispin)*iintfp
        if(ispin.eq.3 .and. iuhf.ne.0) then
        i055=i060
        i060=i055+(noccb+nvrtb)*iintfp
        else
        i055=i050
        endif
        if(i060.gt.maxcor)call insmem('newt2',i060,maxcor)
c
cjp read t1 right hand sides and fix them
cjp for this is essential that second pass of e4s is called AFTER second pass of newt2
        call getlst(icore(i030),1,1,1,2+ispinxx(ispin),90)
        if(ispin.eq.3 .and. iuhf.ne.0) 
     +     call getlst(icore(i035),1,1,1,4,90)
c
cjp but these are not the <...> terms of 4.28 yet, denomin(not zeroed for t1)*t(i,a) must be added to them!
c
        if(ispin.eq.1 .or. ispin.eq.3) 
     +      call fixrhs(icore(i030),icore(i020),icore(i040),icore(i050),
     +         ntaa,pop1,vrt1,nocca,nvrta,1,iuhf)
        if(ispin.eq.2) 
     +      call fixrhs(icore(i030),icore(i020),icore(i040),icore(i050),
     +         ntbb,pop2,vrt2,noccb,nvrtb,2,iuhf)
        if(ispin.eq.3 .and. iuhf.ne.0) 
     +  call fixrhs(icore(i035),icore(i025),icore(i045),icore(i055),
     +         ntbb,pop2,vrt2,noccb,nvrtb,2,iuhf)
c
cjp NOTE: at ecorrbw=0, fixrhs result vanishes for all but internal amplitudes
c
c
c
cjp add rhs contributions to t2 equations, permute Pij Pab
cjp therefore use gtau procedure
        if(correctiontype.eq.0) then
          if(ihubaccorr.lt.2) then
             scalfactor=hfakt
          else
             scalfactor=0.0
          endif
        else
          scalfactor=1.0
        endif
       if(scalfactor.ne.0.0) call gtau(icore(I000),icore(i020),
     +    icore(i030),icore(i025),icore(i035),
     +    dissym,ndssym,pop1,pop2,vrt1,vrt2,irrep,ispin,scalfactor)
c
cjp these are monoexcitations on the LHS of 4.28
cjp this is OK, there is only one permutation and both quantities are t1
cjp so we can use simply ftau
c
cjp this is skipped in a posteriori correction and scaled in iterative one
cjp (ecorrbw is already scaled appropriatelly)
        if(ihubaccorr.lt.2) 
     +      call ftau(icore(I000),icore(i020),icore(i025),
     +         dissym,ndssym,pop1,pop2,vrt1,vrt2,irrep,ispin,-ecorrbw)
        endif
cjp
        endif
cjp
cjp denominator-weight the T2 amplitude increments
cjp
cjp read denominators
        CALL GETLST(ICORE(I010),1,NDSSYM,1,IRREP,LSTDEN)
        call checksum("@-NEWT2=D2 :",ICORE(I010),NDSSYM*DISSYM)
        call checksum("@-NEWT2-T2 :", ICORE(I000), NSIZE)   
        IF(CICALC)CALL CIDENOM(NSIZE,ELAST,ICORE(I010))
cjp
        if(isbwcc) then
cjp ICORE(I000) are t2 increments
        call bwvecprd(ecorrbw,ICORE(I000),ICORE(I010),ICORE(I000),NSIZE)
        else
        CALL VECPRD(ICORE(I000),ICORE(I010),ICORE(I000),NSIZE)
        endif

cjp store new amplitudes
        CALL PUTLST(ICORE(I000),1,NDSSYM,1,IRREP,LSTINC)
        call checksum("@-NEWT2-T2:", ICORE(I000), NSIZE)   
        write(*,*)
10     CONTINUE
c end of the ispin loop
5     CONTINUE

cjp done with pass one
       if(ibwpass.eq.1) return
C
C FORM T2(AA) AMPLITUDES AND DUMP THEM.
cjp in spin-restricted, the aa and bb amplitudes are redundant 
cjp and can be generated from the ab amplitudes
C
      IF(IUHF.EQ.0)THEN
       DO 200 IRREP=1,NIRREP
        NUMAB=IRPDPD(IRREP,ISYTYP(2,63))
        DSZAB=IRPDPD(IRREP,ISYTYP(1,63))
        NUMAA=IRPDPD(IRREP,ISYTYP(2,61))
        DSZAA=IRPDPD(IRREP,ISYTYP(1,61))
        NSIZ1=DSZAB*NUMAA
        NSIZ2=DSZAB*NUMAB
        I000=1
        I010=I000+IINTFP*NSIZ1
        I020=I010+IINTFP*NSIZ2
cjp new t2 amplitudes
        CALL GETLST(ICORE(I010),1,NUMAB,1,IRREP,63)
        CALL ASSYM(IRREP,POP1,DSZAB,DSZAB,ICORE(I000),ICORE(I010))
        CALL SQSYM(IRREP,VRT1,DSZAA,DSZAB,NUMAA,ICORE(I010),ICORE(I000))
        CALL PUTLST(ICORE(I010),1,NUMAA,1,IRREP,61)
200    CONTINUE 
      ENDIF
      RETURN
      END
