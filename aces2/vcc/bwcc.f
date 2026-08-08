











c MR-BW-CCSD extension added by Jiri Pittner (1998-99)

      subroutine bwcc(pCCD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISTSZ,RLECYC
      INTEGER POP(8,2),VRT(8,2)
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,QCISD,LINCC
      LOGICAL CICALC,NONHF,TRIPIT,TRIPNI,TRIPNI1,T3STOR,PRESNT,INIT,
     &        DORESET,UCC,BRUECK,BRKCNV,READT
      LOGICAL ROHF4,ITRFLG,INITL1
      LOGICAL CIS,EOM,DO_HBAR_4LCCSD
      LOGICAL pCCD
      DIMENSION ECORR(3)
      COMMON / / ICORE(1)
      COMMON /ISTART/ I0,ICRSIZ
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /FLAGS/  IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC
      COMMON /TRIPLES/TRIPNI,TRIPNI1,TRIPIT,T3STOR
cAP - 500 is the maximum number of CC iterations (not basis functions)
      COMMON /ENERGY/ ENERGY(500,2),IXTRLE(500)
      COMMON /LINEAR / LINCC,CICALC
      COMMON /NHFREF/ NONHF
      COMMON /ROHF/ ROHF4,ITRFLG
      COMMON /CORENG/ ELAST
      COMMON /BRUECKNER/ BRUECK
      COMMON /SYMLOC/ ISYMOFF(8,8,25)
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
C
      EQUIVALENCE (IFLAGS(2),METHOD)
C
c!!! not necessary to repeat, in fact multiple call of popcorn wastes memory
c MOINTS will be just opened in crapsi, need such one which will exist
      iref = 1
      CALL CRAPSI(ICORE,IUHF,0)

c do not allow crazy attempts
      if (method.ne.8.and.method.ne.10
c     &    .and.method.ne.22
     &    ) then
         print *,
     &      'MR-BW-CC available only within CCD and CCSD approximations'
         call aces_exit(1)
      end if

c we assume that record numbers are the same in all MO files
c this could be done parallelized, but is negligible
      if (nref.gt.1) then
         do iref = 1, nref
            call storemoio
         end do
      end if
      if (.not.masik) call bwprep(nocco,nvrto,iuhf)
      nonhf=.true.

      CALL ZERO(ENERGY,1000)
      CALL IZERO(IXTRLE,500)
      CALL SETMET_VCC(.FALSE.,.FALSE.,.FALSE.,.FALSE.)
      CALL SETLOG
      ICONVG=1
      DO_HBAR_4LCCSD = .FALSE.

      do i = 1, nref
         ibwconvg(i) = 1
      end do

      BRUECK = (IFLAGS(22).EQ.1)
      ITRFLG = .TRUE.
      ROHF4  = .FALSE.

c Nevin added to insure maxcor is aligned from the top as well as the bottom
c      MAXCOR=ICRSIZ
      MAXCOR=ICRSIZ-iand(icrsiz,1)

c ----------------------------------------------------------------------

      do 1234 iref = 1, nref
      if (nref.gt.1) call reopenmo

C INITIALIZE DOUBLES LISTS
      CALL INMBPT(ICORE(I0),MAXCOR,IUHF)

C INITIALIZE SINGLE LISTS
      CALL INSING(METHOD,IUHF)

C SET UP SOME INFORMATION BEFORE CALCULATION.
      IEVERY = IFLAGS(13)
      NCYCLE = IFLAGS(7)
      NKEEP  = IFLAGS(12)
      ICONTL = IFLAGS(4)

      WRITE(6,801)
  801 FORMAT(T3,' Initial T amplitudes: ')
      CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,SING1,'T')
      WRITE(6,800)
  800 FORMAT(T3,' Correlation energies computed from initial ',
     &       'T amplitudes: ')

      call cmpeng(icore(i0),maxcor,43,0,ecorr,edummy1,edummy2,iuhf,1)
      energy(1,1) = 0.d0
      energy(1,2) = 0.d0
      ELAST = ENERGY(1,2)

C BLOCK OF CODE FOR FINITE ORDER CALCULATIONS.
      if (nref.gt.1) then
         call aces_cache_flush
         call aces_cache_reset
         call storemoio
      end if

1234  continue

c ----------------------------------------------------------------------

      do 2341 iref=1,nref
      if (nref.gt.1) call reopenmo

      CALL INMBPT(ICORE(I0),MAXCOR,IUHF)

C SG 7/24/98
C Try to read in a T guess.  If we can, set SING1 to be true.
      CALL DRDTGSS(ICORE(I0),MAXCOR/IINTFP,IUHF,'TGUESS  ',0,READT)

C Note that we need to test to see whether we are doing a CCD
C calculation before we turn the SING1 on. This was not tested
C originaly and CCD optimizations were failing. 01/2001, Ajith Perera
      IF ((IFLAGS(2).NE.5.AND.IFLAGS(2).NE.8).AND.READT) SING1 = .TRUE.

      CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')

C GENERATE W LISTS FOR INTERMEDIATES IF THIS IS NOT MBPT(3) OR LCCD.
      CALL SETLST(ICORE(I0),MAXCOR,IUHF)

      RLECYC = 0
      IF (MOD(IFLAGS(21),2).EQ.0) THEN
         CALL DRRLE(ICORE(I0),MAXCOR,IUHF,RLECYC,.FALSE.)
      END IF
      IF (IFLAGS(21).EQ.1) THEN
c diis probably should be iref specific
         CALL DIISLST(1,IUHF,METHOD.GE.6.AND.METHOD.NE.8)
      END IF
c@@@ here we have a problem: if we allow some stuff to be in core,
c       if the flushcache will really synchronize it into  file!
      CALL INCOR(I0,ICRSIZ,IUHF)
      MAXCOR = ICRSIZ
      if (nref.gt.1) then
         call aces_cache_flush
         call aces_cache_reset
         call storemoio
      end if

2341  continue

c ----------------------------------------------------------------------

      idiisstop = 0
      icycle = 1
555   continue
      RLECYC = RLECYC + 1

      do 3412 iref = 1, nref
      if (nref.gt.1) call reopenmo

C FILL T2 INCREMENTS WITH <IJ||AB> INTEGRALS
C AND ZERO THE T1 INCREMENTS (ONLY CCSD AND QCISD)
      CALL INITIN(ICORE(I0),MAXCOR,IUHF)
      IF (METHOD.GE.6.AND.METHOD.NE.8) THEN
         CALL INITSN(ICORE(I0),MAXCOR,IUHF)
      END IF

C GENERATE W AND F INTERMEDIATES.
      IF (LINCC) THEN
         INTTYP = 3
      ELSE
         INTTYP = 2
      END IF
      CALL GENINT(ICORE(I0),MAXCOR,IUHF,INTTYP,DO_HBAR_4LCCSD)

C COMPUTE F INTERMEDIATE CONTRIBUTION TO T2 INCREMENT.
      CALL FEACONT(ICORE(I0),MAXCOR,IUHF)
      CALL FMICONT(ICORE(I0),MAXCOR,IUHF)
      IF ((METHOD.GT.9.AND.SING1).OR.(METHOD.EQ.6.AND.SING1)) THEN
         CALL FMECONT(ICORE(I0),MAXCOR,IUHF,1)
         IF (IUHF.NE.0) CALL FMECONT(ICORE(I0),MAXCOR,IUHF,2)
         CALL T1INT2A(ICORE(I0),MAXCOR,IUHF)
         CALL T1INT2B(ICORE(I0),MAXCOR,IUHF)
         CALL T1INT1(ICORE(I0),MAXCOR,IUHF,1)
         IF (IUHF.NE.0) CALL T1INT1(ICORE(I0),MAXCOR,IUHF,2)
      END IF

C DO W INTERMEDIATE CONTRIBUTION TO T2 EQUATION
      CALL DRE3EN(ICORE(I0),MAXCOR,IUHF,0)

C COMPUTE THE CONTRIBUTION OF DOUBLES TO T1 AND
C DENOMINATOR WEIGHT THE NEW INCREMENTS
C (ONLY CCSD AND QCISD METHODS)
      IF (METHOD.GE.9.OR.METHOD.EQ.6.OR.METHOD.EQ.7) THEN
c set flag to perform only first part of processing in e4s, e4seng
         ibwpass=1
         CALL E4S(ICORE(I0),MAXCOR,IUHF,EDUMMY)
      END IF

C DENOMINATOR WEIGHT T2 INCREMENTS TO FORM NEW T2.
C
c NEWT2 uses common bwcc and conditionals
c set flag to perform only first part of processing
      ibwpass=1
      CALL NEWT2(ICORE(I0),MAXCOR,IUHF)

      if (nref.gt.1) then
         call aces_cache_flush
         call aces_cache_reset
         call storemoio
      end if

3412  continue

c communicate columns of heff between processors
c in parallel calculation
c diagonalize heff
      call heffdiag(iuhf,icycle)
      write(6,*)'@VCC: MR-BW-CC CORRELATION ENERGY = ',epsilon0
      write(6,*)

c store correctly in energy array, ecorr variable
      ixcycle = -1
      do 4123 iref = 1, nref
      if (nref.gt.1) call reopenmo

c setup BW shift of denominator
      if (ihubaccorr.lt.2) then
c ordinary BWCC value
         ecorrbw=epsilon0-heff(iref,iref)
      else
c FOR a posteriori size extenzivity correction make it like in RSPT-CC
         ecorrbw=0.d0
      end if
      ecorrbw0=epsilon0-heff(iref,iref)
c homotopical transition between BW and RS
      hfakt=1.0
      if (ihomotop.gt.0) then
         if (icycle.gt.ihomotop) then
c first such iteration ... switch off diis
            if (icycle.eq.ihomotop+1) idiisstop = 1
            hfakt=lambdahomotop**(icycle-ihomotop)
c restart diis eventually
            if (hfakt.lt.diishonset.and.
     &          idiisstop.eq.1.and.
     &          icycle.gt.ihomotop+1) then
c postpone this in order to call diislst for other references, too
c but set a flag that in this cycle, dodiis0 should already be called
               ixcycle=icycle
               if (iref+1.gt.nref) idiisstop=0
c remember, we are still inside iref loop
               IF (IFLAGS(21).EQ.1) THEN
                  write(6,*)'DIIS RESTART',icycle,hfakt,diishonset
                  CALL DIISLST(1,IUHF,METHOD.GE.6.AND.METHOD.NE.8)
               END IF
            end if
         end if
         ecorrbw0=epsilon0-heff(iref,iref)
         ecorrbw=hfakt*ecorrbw0
         if (icycle.gt.ihomotop) then
            write(6,*) 'homotopic solution: ',
     &             'cycle,hfakt,scaled ecorrbw=',icycle,hfakt,ecorrbw
         end if
      end if
c here we must finish the work of e4seng, newt2
      ibwpass=2
      call newt2(icore(i0),maxcor,iuhf)
c NOTE: PUTLST(ICORE(ID),1,1,1,3,90) contains elements needed for  eq 4.28
c which will be destroyed and replaced by final t1 amplitudes in e4s call
c thus newt2 is here called first, in contrast to original code
c@@@!!!note they correspond to NEW t1 amplitudes!
      call e4s(icore(i0),maxcor,iuhf,edummy)
C
C WRITE OUT RMS AND MAX DIFFERENCES FOR NEW T2 AND PUT NEW ERROR VECTOR
      call drtsts(icore(i0),maxcor,icycle,iuhf,
     &            ibwconvg(iref),icontl,sing1,0,'t')
c this computes already new diagonal heff elements to be used in next iter
c NOTE this overwrites heff
c in energy array store the true mr-cc energy!
      call cmpeng(icore(i0),maxcor,60,2,ecorr,
     &            edummy1,edummy2,iuhf,1)
      energy(icycle+1,1)=epsilon0
      energy(icycle+1,2)=epsilon0
      ELAST = ENERGY(ICYCLE+1,2)
C
C IF ICONVG IS EQUAL ZERO, CONVERGENCE HAS BEEN ACHIEVED, EXIT
C VIA FINISH
C
c switch off diis for a posteriori correction
      IF (IFLAGS(21).EQ.1.and.ihubaccorr.lt.2) THEN
c bwcc case
c NOTE: strictly taken, we do illegal thing by calling diis INDEPENDENTLY
c for each iref,
c we should better build a common vector for all irefs and do ONE diis
c however, in the cases test up to now it did not matter too much, perhaps
c in last iterations with strict conv. criterion is wastes 1-2 iterations
         if (idiisstop.eq.0 .or. icycle.eq.ixcycle) then
            call dodiis0(icore(i0),maxcor/iintfp,iuhf,1,icycle,
     &                   ibwconvg(iref),icontl,sing1,44,61,90,0,90,2,
     &                   70,  '     ',1.0D0)
         end if
      END IF

c NOTE this iref loop has been split in order to run the amplitude
c update consistently for all references
      if (nref.gt.1) then
         call aces_cache_flush
         call aces_cache_reset
         call storemoio
      end if

4123  continue
c
c TEST OF CONVERGENCE SHOULD BE OUTSIDE THE IREF LOOP
      if (ihomotop.gt.0.and.hfakt.gt.hfaktmax) then
c force further iterations if homotopy is not completed yet
         iconvg=1
      else
c test now on all processors
         iconvg=0
      end if
c iconvg zero, if all of them zero for all configs
      do iref = 1, nref
         if (ibwconvg(iref).ne.0) iconvg=1
      end do

c here is a place to exit from mr-bwcc loop
      if (iconvg.eq.0) then
         write(6,*)'@VCC: the MR-BW-CC calculation has converged'
      end if

c implementation of the Hubac size extenzivity correction for BWCC
      if (ihubaccorr.eq.1.and.iconvg.eq.0.or.
     &    ihubaccorr.gt.1) then
c pass one ... switch off diis for the future, compute another iteration in
c which RSPT denominators will be used ... pass 2
c in pass 3, the new amplitudes obtained using RSPT denominators will be used
c to obtain new Heff and the final energy by its diagonalization
         if (ihubaccorr.eq.3) then
            iconvg=0
         end if
         if (ihubaccorr.eq.2) then
            iconvg=1
            ihubaccorr=3
         end if
         if (ihubaccorr.eq.1) then
            write(6,*) 'Computing the size-extenzivity correction'
            iconvg=1
            ihubaccorr=2
         end if
      end if

c start homotopy if convergence has been reached
      if (ihomotop.gt.0.and.iconvg.eq.0.and.
     &    icycle.lt.ihomotop) then
         ihomotop=icycle
         iconvg=1
      end if

      do 4321 iref=1,nref
      if (nref.gt.1) call reopenmo
C
CJDW 10/24/95.
C    Do not update if we have converged. This allows us to maintain a
C    correspondence between the amplitudes and existing Hbar elements.
C    This is important for MN's finite-order EOM stuff.
C
c for mrbwcc, we also do not overwrite the amplitudes, from which the current
c energy was computed
c use the all-iref iconvg here!
      IF (ICONVG.NE.0) THEN
         CALL DRMOVE(ICORE(I0),MAXCOR,IUHF,0,SING1)
         CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')
      END IF
CJDW END
c when convergence was achieved
      IF (ICONVG.EQ.0) THEN
c the cmpeng call is probably obsolete here
c  call cmpeng(icore(i0),maxcor,43,0,ecorr,edummy1,edummy2,iuhf,1)
         energy(icycle+1,1)=epsilon0
         energy(icycle+1,2)=epsilon0
      END IF
C
      IF (ICONVG.EQ.0) THEN

CJDW 3/20/96. Print final amplitudes.
         CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,SING1,'T')

C  PUT CONVERGED AMPLITUDES AS INITIAL GUESS FOR SOLVING THE
C  LAMBDA EQUATIONS ON THE GAMLAM FILE
         WRITE(6,4000)
 4000    FORMAT(T3,' The CC iterations have converged.',/)
         IF (GRAD.OR.
     &       (IFLAGS(87).EQ.3 .AND. (TRIPIT.OR.TRIPNI1) ) ) THEN
            CALL INIT2(IUHF,PCCD)
            CALL DRMOVE(ICORE(I0),MAXCOR,IUHF,100,SING1)
         END IF
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
c this should be called for all iref
c leave it here in the iref loop 5
            IF(IUHF.EQ.0)CALL RESET(ICORE(I0),MAXCOR,IUHF)
         END IF
C
         IF (IFLAGS(87).EQ.3.AND.(TRIPIT.OR.TRIPNI)) THEN
            CALL RESRNG(ICORE(I0),MAXCOR,IUHF)
         END IF
C-----------------------------------------------------------------------
c this has sense only to be called once!
         if (iref+1.gt.nref) then
            CALL FINISH(ICYCLE+1)
         end if
C SG 7/21/98
         CALL DDMPTGSS(ICORE(I0), MAXCOR/IINTFP, IUHF, 0, 'TGUESS  ')

         IF (TRIPNI.AND.BRKCNV) THEN
            CALL TRPS(ICORE(I0),MAXCOR,IUHF,EDUMMY)
         END IF
         IF (IFLAGS(35).NE.0) CALL ACES_AUXCACHE_FLUSH

         if (nref.gt.1) then
            call aces_cache_flush
            call aces_cache_reset
            call storemoio
         end if

c go through whole iref loop
         if (iref+1.gt.nref) then
c ... for eneronly geometry optimizations - either in each program
c     gfname must know that MOINTS are MOINTS.n
c     or simpler and at the moment sufficient, MOINTS will be deleted
            istat = ishell('rm -f MOINTS.* MOABCD.*')
            write(6,*) '@VCC-I: files MOINTS.* MOABCD.* removed'
            WRITE(6,1020)
            call aces_fin
            stop
         end if

c     END IF (ICONVG.EQ.0)
      END IF

C DO THE RLE EXTRAPOLATION AND GET THE UPDATED RLE ENERGY.
C
c no EXTRAPOLATION for the a posteriori correction
      IF (MOD(IFLAGS(21),2).EQ.0 .and. ihubaccorr.lt.2) THEN
         CALL DRRLE(ICORE(I0),MAXCOR,IUHF,RLECYC,.FALSE.)
         CALL DRMOVE(ICORE(I0),MAXCOR,IUHF,0,SING1)
         CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')
      END IF
C
C PRINT AMPLITUDES EVERY IEVERY CYCLES.
C
c more default print for debug
c print the true converged amplitudes if correction is to be calculated
c since they will be overwritten by amplitudes obtained using the
c modified denominator
      IF (iconvg.ne.0 .and. ihubaccorr.ge.2 .or.
     &    (ihomotop.gt.0).and.iconvg.ne.0.or.
     &    IFLAGS(1).GE.10 .and. iconvg.ne.0) THEN
         IF (ihubaccorr.ge.2 .or. MOD(ICYCLE,IEVERY).EQ.0) THEN
            if (ihubaccorr.eq.2) then
               write(6,*) ' Ordinary BWCC converged amplitudes'
            end if
            if (ihubaccorr.eq.3) write(6,*) ' Corrected BWCC amplitudes'
            CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,SING1,'T')
         END IF
      END IF
c
      if (nref.gt.1) then
         call aces_cache_flush
         call aces_cache_reset
         call storemoio
      end if

4321  continue

      icycle=icycle+1
      if (icycle.le.ncycle) goto 555

      WRITE(*,*)
     &         '@VCC: The Coupled-Cluster equations did not converge!!!'
      do iref = 1, nref
         if (nref.gt.1) call reopenmo
         CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,SING1,'T')
         if (nref.gt.1) then
            call aces_cache_flush
            call aces_cache_reset
            call storemoio
         end if
      end do
      call aces_fin
      call aces_exit(1)

 1020 FORMAT(/,77('-'),/,32X,'Exiting xvcc',/,77('-'),/)
      WRITE(6,1020)

      END

