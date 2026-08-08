










      SUBROUTINE RCC_FINISH(ICYCLE)
C
C RUN-DOWN ROUTINE FOR CC CALCULATIONS.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*6 STAR
      CHARACTER*13 ITCALC(50)
      LOGICAL LINCC,CICALC

      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
cjp
      logical o
      COMMON /LINEAR/ LINCC,CICALC
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /ENERGY/ ECORR(500,2),IXTRLE(500)
      COMMON /FLAGS/ IFLAGS(100)
      DATA ITCALC /
     & '             ','             ','             ','             ',
     & '  LCCD       ','  LCCSD      ',' UCCSD(4)    ','   CCD       ',
     & '  UCC(4)     ','  CCSD       ','  CCSD       ','  CCSD       ',
     & '  CCSDT-1a   ','  CCSDT-1b   ','  CCSDT-2    ','  CCSDT-3    ',
     & '  CCSDT-4    ','  CCSDT      ','  LCCSDT     ','   CCD       ',
     & '  QCISD      ','  CCSD       ','  QCISD      ','   CID       ',
     & '   CISD      ','  QCISD      ','  CCSD       ','  CCSD       ',
     & '  CCSDT      ','  CCSDT      ','  CCSD       ','  CCSD       ',
     & '     CC3     ',' CCSDT-T1T2  ','  CCSDTQ-1   ','  CCSDTQF-1  ',
     & '  CCSDTQ-2   ',' CCSDTQ-3    ','  CCSDTQ     ','  ACCSD      ',
     & '             ',' ACCSD(T)    ','  CCSD(TQf)  ','  CCSDT(Qf)  ',
     & ' OO-MP(2)    ',' OO-MBPT(2)  ','  CC2        ','  rCCD       ',
     & ' drCCD       ','             '/
C
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
cjp in parallel run, assure that this printout will appear 
cjp at the end of output
c
cjp output only from one processor
      o=.true.

      CALL GETREC(20,'JOBARC','SCFENEG ',IINTFP,ESCF)
      if(o)WRITE(6,99)
99    FORMAT(T3,'     Summary of iterative solution of CC equations ')
      if(o)WRITE(6,100)
100   FORMAT(6X,59('-'))
      if(o)WRITE(6,101)
101   FORMAT(T24,'Correlation',T46,'Total',/,T8,'Iteration',T26,
     &       'Energy',T46,'Energy')
      if(o)WRITE(6,100)
      DO 20 I=1,ICYCLE
       STAR='JACOBI'
       IF(IXTRLE(I).EQ.1)STAR=' RLE  '
       IF(IFLAGS(21).EQ.1)STAR=' DIIS '
cjp write a note about correction of BWCC
       if(isbwcc.and.I.EQ.ICYCLE.and.ihubaccorr.eq.3) STAR=' CORR '
cjp write a note about homotopic correction
       if(isbwcc.and.i-1.gt.ihomotop.and.ihomotop.gt.0) then
        hfakt=lambdahomotop**(i-1-ihomotop)
        if(hfakt .lt. diishonset .and. IFLAGS(21).EQ.1 ) then
           STAR=' HDIIS'
        else
           STAR='HJACOB'
        endif
       endif
       IF(I.EQ.ICYCLE.AND.UCC)THEN
        if(o)WRITE(6,1001)
1001    FORMAT(T15,' Adding <T2^+T2^+WT2> energy contribution.')
       ENDIF
cjp info about homotopicsolution
       if(isbwcc .and. i-1.gt.ihomotop.and.ihomotop.gt.0) then
          if(o)WRITE(6,1002)I-1,ECORR(I,1),ECORR(I,1)+ESCF,STAR,hfakt
       else
          if(o)WRITE(6,1000)I-1,ECORR(I,1),ECORR(I,1)+ESCF,STAR
       endif
1000   FORMAT(T10,I4,T19,F18.12,T39,F19.12,T59,A6)
1002   format(T10,I4,T19,F18.12,T39,F19.12,T59,A6,1x,e9.3)
20    CONTINUE
      if(o)WRITE(6,100)
C
CJDW 10/1/96. Write out energy for iterative part of method.
C
cjp
       if(isbwcc ) then
         if(ihubaccorr.eq.3) then
        if(o)WRITE(6,1010) 'cMR-BW-'//ITCALC(IFLAGS(2)),
     &     ECORR(ICYCLE,1)+ESCF
         else
          if(ihomotop.gt.0) then
        if(o)WRITE(6,1010) 'hMR-BW-'//ITCALC(IFLAGS(2)),
     &    ECORR(ICYCLE,1)+ESCF
          else
         if(o) WRITE(6,1010) 'MR-BW-'//ITCALC(IFLAGS(2)),
     &ECORR(ICYCLE,1)+ESCF
          endif
         endif
        else
       if(o)WRITE(6,1010) ITCALC(IFLAGS(2)),ECORR(ICYCLE,1)+ESCF
        endif
C
      IF(.NOT.CICALC)THEN
cjp
       if(isbwcc) then
        if(o)write(6,210)
        else
       if(o)WRITE(6,202)
        endif
      ELSE
       if(o)WRITE(6,203)
       if(o)WRITE(6,205)
       if(o)WRITE(6,204)
      ENDIF
cjp summarize serious warnings
      do i=1,maxbwwarnings
      if(bwwarning(i).and.o)write(6,*) bwwarntext(i)
      enddo
c write total energy to jobarc
      CALL PUTREC(1,'JOBARC','TOTENERG',IINTFP,ECORR(ICYCLE,1)+ESCF)
cmn again for vibronic purposes:
      CALL PUTREC(1,'JOBARC','TOTENER2',IINTFP,ECORR(ICYCLE,1)+ESCF)
      CALL PUTREC(1,'JOBARC','PARENERG',IINTFP,ECORR(ICYCLE,1)+ESCF)
c
210   FORMAT(T7,'A multi-miracle has come to pass. ',
     &          'The MR-BW-CC iterations have converged.')
202   FORMAT(T7,'A miracle has come to pass. ',
     &          'The CC iterations have converged.')
203   FORMAT(T7,'A miracle has come to pass. ',
     &          'The CI iterations have converged.')
205   FORMAT(/,T12, ' ****************WARNING!!!!************* ')
204   FORMAT(T3,'For the most part, the authors of this program ',
     &         'feel that CI is ',/
     &         ,T3,'a bad method.  Please run CC calculations in the',
     &         ' future!',/)
 1010 FORMAT(/,T14,A13,' energy is ',F20.12,' a.u. ',/)
      RETURN
      END
