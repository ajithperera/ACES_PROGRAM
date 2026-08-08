




























































































































































































      SUBROUTINE PCCD_URDRIVER(ICORE,MAXCOR,IUHF)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ICORE(MAXCOR)
      DIMENSION ECORR(3)
      INTEGER RELCYC

      LOGICAL NONHF
      LOGICAL PCCD,CCD,LCCD
      LOGICAL CIS,EOM
      LOGICAL BRUECK,BRKCNV
      LOGICAL DORESET,OO_CC
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,QCISD
      LOGICAL UCC,CC2
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
c files.com : begin
      integer        luout, moints
      common /files/ luout, moints
c files.com : end

      COMMON/CALC/PCCD,CCD,LCCD
      COMMON /ENERGY/ ENERGY(500,2)
      COMMON/EXCITE/EOM,CIS
      COMMON/NHFREF/NONHF
      COMMON/EXTRAPO/RLE
CSSS      COMMON /SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
CSSS     &                QCISD,UCC,CC2

      CALL ZERO(ENERGY,1000)
      CALL PCCD_SETMET()
      I0      = 1
      ICONTL  = IFLAGS(4)
      NCYCLE  = IFLAGS(7)
      NKEEP   = IFLAGS(12)
      IEVERY  = IFLAGS(13)
      BRUECK  = (IFLAGS(22) .EQ. 1)
      METHOD  =  IFLAGS(2)
      DAMP_PARAMETER = IFLAGS2(170)*0.010D0
      SING1   = .FALSE.
      OO_CC   = (IFLAGS(2) .EQ. 52 .OR.
     &           IFLAGS(2) .EQ. 53 .OR.
     &           IFLAGS(2) .EQ. 54)
      LUOUT   = 6
      MOINTS  = 50
C
      CALL SETLST(ICORE(I0),MAXCOR,IUHF)
      CALL INSING(10,IUHF)
      ICYCLE = 1
      RLECYC = 0
      CALL INMBPT(ICORE(I0),MAXCOR,IUHF)
      IF (IFLAGS(21) .EQ.1) THEN
        CALL PCCD_DIISLST(1,IUHF,.FALSE.,"R")
      END IF
      CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')

      ICONVG = 1

      WRITE(6,801)
  801 FORMAT(T3,' Initial T amplitudes: ')
      CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,SING1,'T')
      WRITE(6,800)
  800 FORMAT(T3,' Correlation energies computed from initial ',
     &          'T amplitudes: ')

      CALL CMPENG(ICORE(I0),MAXCOR,43,0,ECORR,ENERGY(1,1),ENERGY(1,2),
     &            IUHF,1)
      ELAST = ENERGY(1,2)

 100  CONTINUE
      RLECYC = RLECYC + 1
C
C Start the T2 lists by adding <ab||ij>, This is actually a
C contribution to T2.
C
      CALL INITIN(ICORE(I0),MAXCOR,IUHF)

       call check_t2(icore(i0),Maxcor,Iuhf) 
      IF (METHOD.GE.6.AND.METHOD.NE.8) THEN
         CALL INITSN(ICORE(I0),MAXCOR,IUHF)
      END IF
c
C GENERATE W AND F INTERMEDIATES.
C
      IF (LINCC) THEN
         INTTYP = 3
      ELSE
         INTTYP = 2
      END IF
C 
C Generate W(mbej), W(mnij) W(mbej, F(ea), F(ij) and F(me)
C intermediates.  
C
      Write(6,"(a)") " The checksums of all the intms before genint"
      Write(6,"(a,2l)") " Flags singles and non-hf:",SING1,NONHF
      Call checkintms(icore(i0),maxcor,Iuhf,1)
      write(6,"(a)") " The T1/T2 residuals" 
      call check_t2(icore(i0),Maxcor,Iuhf) 
C
      CALL GENINT(ICORE(I0),MAXCOR,IUHF,INTTYP,DO_HBAR_4LCCSD)
C
C COMPUTE F INTERMEDIATE CONTRIBUTION TO T2 INCREMENT.
C
C Evaluate P(ab)Sum_e T(ij,ae){F(b,e)-1/2 sum_m T(m,b)F(m,e)}
C contribution to T2 and Sum_e T(i,e)F(a,e) contribution to T1.
C
      CALL FEACONT(ICORE(I0),MAXCOR,IUHF)
       call check_t2(icore(i0),Maxcor,Iuhf)
C
C Evaluate -P(ij)Sum_e T(im,ab){F(mj)-1/2 sum_m T(j,e)F(m,e)}
C contribution to T2 and Sum_e T(m,a)F(m,i) contribution to T1.
C
      CALL FMICONT(ICORE(I0),MAXCOR,IUHF)
       call check_t2(icore(i0),Maxcor,Iuhf)
C
      IF ((METHOD.GT.9.AND.SING1).OR.(METHOD.EQ.6.AND.SING1)) THEN

C Evaluate Sum_e T(im,ae)F(m,e) contribution to T1. 
C
       CALL FMECONT(ICORE(I0),MAXCOR,IUHF,1)
       IF (IUHF.NE.0) CALL FMECONT(ICORE(I0),MAXCOR,
     &                                  IUHF,2)
C
       call check_t2(icore(i0),Maxcor,Iuhf)
C 
C Evaluate -P(ab)Sum_m T(m,a)<mb||ij> contribution to T2.
C 
         CALL T1INT2A(ICORE(I0),MAXCOR,IUHF)
C
C Evaluate +P(ij)Sum_m T(i,e)<ab||ej> contribution to T2.
C
         CALL T1INT2B(ICORE(I0),MAXCOR,IUHF)
C
C Evaluate -Sum_nf T(n,f)<na||if> contribution to T1.
C
         CALL T1INT1(ICORE(I0),MAXCOR,IUHF,1)
         IF (IUHF.NE.0) CALL T1INT1(ICORE(I0),MAXCOR,
     &                             IUHF,2)
C
       call check_t2(icore(i0),Maxcor,Iuhf)
      END IF
C 
C DO W INTERMEDIATE CONTRIBUTION TO T2 EQUATION
C 
C Evaluate the following three contributions to T2:
C
C   -P(ij)P(ab) Sum_me T(i,e)T(m,a) <mb||ej> (T12INT2)
C   +P(ij)P(ab) Sum_me T(i,e)T(m,a) W(mb,ej) (DRRNG)
C   +1/2 Sum_mn Tau(mn,ab)W(mn,ij) + 
C    1/2 Sum_ef Tau(ij,ef)W(ab,ef)           (DRLAD)
C
       call check_t2(icore(i0),Maxcor,Iuhf) 
      CALL DRE3EN(ICORE(I0),MAXCOR,IUHF,0)
C 
C
C COMPUTE THE CONTRIBUTION OF DOUBLES TO T1 AND
C DENOMINATOR WEIGHT THE NEW INCREMENTS
C (ONLY CCSD AND QCISD METHODS)
C 
      IF ((METHOD.GE.9.OR.METHOD.EQ.6.OR.METHOD.EQ.7)
     &    .AND. .NOT. OO_CC) THEN
C
C Evaluate the following three contributions to T1:
C
C   -1/2 Sum_mef T(im,ef)<ma||ef>  (T2T1AA1, T2T1AB1)
C   -1/2 Sum_men T(nm,ei)<nm||ei>  (T2T1AA2, T2T1AB2)
C
       call check_t2(icore(i0),Maxcor,Iuhf) 
C Also, Do the T1 = T1/{f(i,i) - f(a,a)} to get a new T1
C
         CALL E4S(ICORE(I0),MAXCOR,IUHF,EDUMMY)
      END IF
      CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,SING1,'T')

C Also the call E4S for pCCD (OO-CCD) and Brueckner. This will do T2 into T1
C terms. Ajith Perera,10/2021. 

       IF (OO_CC .AND. BRUECK) THEN
          CALL E4S(ICORE(I0),MAXCOR,IUHF,EDUMMY)
       ENDIF 
C
C DENOMINATOR WEIGHT T2 INCREMENTS TO FORM NEW T2.
C 
C Do the T2 = T2/{f(i,i) + f(j,j) - f(a,a) - f(b,b)}
C to get a new T2.
C
      IF (pCCD) THEN
C Ajith Perera, 10/2021.
         Write(6,"(20x,a)") " ----------Warning------------"
         Write(6,"(2a)")  " A pCCD or pCCD like calculation is being",
     +                    " performed and the off-diagonal"
         Write(6,"(a)")   " blocks of T2 is set to zero"
         Write(6,"(20x,a)") " -----------------------------"
         call pccd_reset_vcc(icore(i0),Maxcor,Iuhf,63)
      ENDIF 
       call check_t2(icore(i0),Maxcor,Iuhf)
      CALL NEWT2(ICORE(I0),MAXCOR,IUHF)

      CALL CMPENG(ICORE(I0),MAXCOR,60,2,ECORR,ENERGY(ICYCLE+1,1),
     &            ENERGY(ICYCLE+1,2),IUHF,1)
      ELAST = ENERGY(ICYCLE+1,2)
C
C IF ICONVG IS EQUAL ZERO, CONVERGENCE HAS BEEN ACHIEVED, EXIT
C VIA FINISH
C
      IF (IFLAGS(21) .EQ.1) THEN
         CALL PCCD_DODIIS0(ICORE(I0),MAXCOR/IINTFP,IUHF,1,ICYCLE,
     &                ICONVG,ICONTL,.FALSE.,44,61,90,0,90,2,70)
      END IF

      IF (ICONVG.NE.0) THEN
         CALL DRMOVE(ICORE(I0),MAXCOR,IUHF,0,SING1)
         CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')
      END IF

      IF (ICONVG.EQ.0) THEN
         CALL CMPENG(ICORE(I0),MAXCOR,43,0,ECORR,ENERGY(ICYCLE+1,1),
     &               ENERGY(ICYCLE+1,2),IUHF,1)
         CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,SING1,'T')
C
CSSS         CALL FINISH(ICYCLE+1,PCCD,.FALSE.,.FALSE.)
         CALL PCCD_UFINISH(ICYCLE+1,"R")
         CALL INIT2(IUHF,PCCD.OR.CCD.OR.LCCD)
         CALL DRMOVE(ICORE(I0),MAXCOR,IUHF,100,SING1)
         CALL DDMPTGSS(ICORE(I0), MAXCOR/IINTFP, IUHF, 0, 'TGUESS  ')
C
         IF (DORESET) THEN
            IF(IUHF.EQ.0)CALL RESET(ICORE(I0),MAXCOR,IUHF)
         END IF 
        WRITE(6,1020)
 1020 FORMAT(/,77('-'),/,32X,'Exiting pccd_urdriver',/,77('-'),/)

         RETURN
      ENDIF 
C
      ICYCLE = ICYCLE +  1
      IF (ICYCLE .LE. NCYCLE) GOTO 100

      IF (ICONVG .NE. 0) THEN
         CALL PCCD_AMPSUM(ICORE(I0),MAXCOR,IUHF,0,.FALSE.,'T')
         IF (pCCD) THEN
         Write(6,"(2a)") "The pCCD equations did not converge",
     &                     " in alloted number of cyclces."
         ELSE IF(CCD) THEN
         Write(6,"(2a)") "The CCD equations did not converge",
     &                     " in alloted number of cyclces."
         ELSE IF(LCCD) THEN
         Write(6,"(2a)") "The LCCD equations did not converge",
     &                     " in alloted number of cyclces."
         ENDIF
      ENDIF



      RETURN
      END

