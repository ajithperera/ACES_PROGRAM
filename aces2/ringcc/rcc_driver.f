










      SUBROUTINE RCC_DRIVER(WORK,MAXCOR,IUHF)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION WORK(MAXCOR)
      DIMENSION ECORR(3)
      INTEGER RELCYC
     
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,QCISD,
     +        UCC,CC2,RCCD,DRCCD,RPA_T2VECS
C
C A ring CCD (rCCD) and direct ring CCD (drCCD) as described in JCP 139, 
C 104113 (2013) is implemented. 
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

      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
      COMMON /ENERGY/ ENERGY(500,2),IXTRLE(500)
      COMMON/T2_SOURCE/RPA_T2VECS

C Note at the first-iteration, T1 vector is zero and T2 is MBPT(2).

      CALL RCC_SETMET 
C
C Initializations of various intermediate lists etc.
C
      INEXT   = 1
      RLECYC  = 0
      IREF    = 1
      IBWPASS = 0
      ICONTL  = IFLAGS(4)
      NCYCLE  = IFLAGS(7)
      NKEEP   = IFLAGS(12)
      IEVERY  = IFLAGS(13)

C For UHF we have to work with Coulomb only integrals and symmetric 
C amplitudes. The entire ACES II structure built on antisymetric 
C integals and amplitudes break down. So, AAAA and BBBB lists 
C and routines that supports operation on them have been readjusted 
C to work with symmetric integrals and amplitudes. This of course
C matters only for UHF. RHF calcs (as long as they can be solely 
C determined by ABAB combination) can proceed normally. 

      CALL RCC_SETLST(WORK(INEXT),MAXCOR,IUHF)

C Check whether r/dr-CCD vectors are on lists 44,45 and 46 from
C the RPA code.

      IF (RPA_T2VECS) THEN
         CALL RCC_RNABIJ(WORK(INEXT),MAXCOR,IUHF,'T')
         CALL RCC_GENINT(WORK(INEXT),MAXCOR,IUHF)
         RETURN
      ENDIF 

      IF (IUHF .NE. 0) THEN 
         CALL RCC_DRIVE_D2T2(WORK(INEXT),MAXCOR,IUHF)
         CALL RCC_RNABIJ(WORK(INEXT),MAXCOR,IUHF,'T')
      ELSE
         IF (RCCD) THEN
            NSIZE = IDSYMSZ(1,ISYTYP(1,46),ISYTYP(2,46))
            CALL GETALL(WORK(INEXT),NSIZE,1,46)
            CALL PUTALL(WORK(INEXT),NSIZE,1,44)
            CALL RCC_ABIJ_4RHF(WORK(INEXT),MAXCOR,IUHF)
         ENDIF 
         CALL RCC_RNABIJ(WORK(INEXT),MAXCOR,IUHF,'T')
      ENDIF 

      CALL RCC_AMPSUM(WORK(1),MAXCOR,IUHF,0,.FALSE.,'T')
CSSS      call check_t2(work(1),maxcor,iuhf)
CSSS      CALL RCC_SETLST(WORK(INEXT),MAXCOR,IUHF)

      CALL RCC_CMPENG(WORK(INEXT),MAXCOR,43,0,ECORR,ENERGY(1,1),
     &                ENERGY(1,2),IUHF,1)
      ELAST = ENERGY(1,2)

      IF (MOD(IFLAGS(21),2).EQ.0) THEN
         CALL DRRLE(WORK(INEXT),MAXCOR,IUHF,RLECYC,.FALSE.)
      END IF
      IF (IFLAGS(21) .EQ.1) THEN
         CALL DIISLST(1,IUHF,.FALSE.)
      END IF

      ICYCLE = 1
      DAMP_PARAMETER = IFLAGS2(170)*0.010D0

      Write(*,*)
      write(*,"(3x,a,F6.4)") "The RCC damp parameter = ", DAMP_PARAMETER
      Write(*,*)

 100  CONTINUE
     
      RELCYC = RLECYC + 1

      IF (DRCCD) THEN
         CALL DRCC_INITIN(WORK(INEXT),MAXCOR,IUHF)
      ELSE IF (RCCD) THEN
         CALL RCC_INITIN(WORK(INEXT),MAXCOR,IUHF)
      ENDIF 

C Construct W(mb,ej)= <mb||ej> -1/2T2(jn,fb)<mn||ef> intermediate

      CALL RCC_GENINT(WORK(INEXT),MAXCOR,IUHF)
C
C Construct +P(ij)P(ab) Sum_me T2(im,ae) W(mb,ej) (DRRNG)
C
      CALL RCC_DORNG(WORK(INEXT),MAXCOR,IUHF)
C
C Denominator weigh the T2 increments. 
C
      CALL RCC_NEWT2(WORK(INEXT),MAXCOR,IUHF)

      IF (DAMP_PARAMETER .NE. 1.0D0) THEN
         CALL DAMP_CC_RESIDUAL(WORK(INEXT),MAXCOR,IUHF,.FALSE.,
     &                         44,61,90,90,DAMP_PARAMETER)
      ENDIF
C
      CALL DRTSTS(WORK(INEXT),MAXCOR,ICYCLE,IUHF,ICONVG,ICONTL,
     &            .FALSE.,0,'T')
      CALL RCC_CMPENG(WORK(INEXT),MAXCOR,60,2,ECORR,
     &                ENERGY(ICYCLE+1,1),ENERGY(ICYCLE+1,2),IUHF,1)

      ELAST = ENERGY(ICYCLE+1,2)
C
      IF (IFLAGS(21) .EQ.1)  THEN
         CALL RCC_DODIIS0(WORK(INEXT),MAXCOR,IUHF,1,ICYCLE,
     &                    ICONVG,ICONTL,.FALSE.,44,61,90,0,90,2,70,
     &                    '     ',DAMP_PARAMETER)
      END IF

      IF (ICONVG.NE.0) THEN

         CALL DRMOVE(WORK(INEXT),MAXCOR,IUHF,0,.FALSE.)
         CALL RCC_RNABIJ(WORK(INEXT),MAXCOR,IUHF,'T')

      END IF

      IF (ICONVG.EQ.0) THEN

         IF (DRCCD) CALL RCC_ABTOAA(WORK(INEXT),MAXCOR,IUHF,61)
         CALL RCC_CMPENG(WORK(INEXT),MAXCOR,43,0,ECORR,
     &                   ENERGY(ICYCLE+1,1),ENERGY(ICYCLE+1,2),
     &                   IUHF,1)
     
         CALL RCC_AMPSUM(WORK(INEXT),MAXCOR,IUHF,0,.FALSE.,'T')
         IF (IUHF.NE.0) CALL S2PROJ(WORK(INEXT),MAXCOR,IUHF,.FALSE.)
C
         CALL RCC_FINISH(ICYCLE+1)
         CALL DDMPTGSS(WORK(INEXT), MAXCOR/IINTFP, IUHF, 0, 'TGUESS  ')
         CALL ACES_FIN 
         WRITE(6,1020)
 1020 FORMAT(/,77('-'),/,32X,'Exiting xvcc',/,77('-'),/)

         CALL ACES_EXIT(0)

      ENDIF 

      IF (MOD(IFLAGS(21),2).EQ.0) THEN
         CALL DRRLE(WORK(INEXT),MAXCOR,IUHF,RLECYC,.FALSE.)
         CALL DRMOVE(WORK(INEXT),MAXCOR,IUHF,0,.FALSE.)
         CALL RCC_RNABIJ(WORK(INEXT),MAXCOR,IUHF,'T')
      END IF
C
      ICYCLE = ICYCLE +  1
      IF (ICYCLE .LE. NCYCLE) GOTO 100
       
      IF (ICONVG .NE. 0) THEN
         CALL RCC_AMPSUM(WORK(INEXT),MAXCOR,IUHF,0,.FALSE.,'T')
         Write(6,"(2a)") "The RCCD or DRCCD equations did not converge",
     &                     " in alloted number of cyclces."
         call aces_exit(1)
      ENDIF 
C
      RETURN
      END
