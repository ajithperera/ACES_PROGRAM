










      SUBROUTINE RCC_DODIIS0(SCR,MAXCOR,IUHF,IRREPX,ICYCLE,
     &                       ICONVG,ICONTL,SING,LISTT20,LISTT21,
     &                       LISTT10,IOFF0,LISTT11,IOFF1,LSTERR,
     &                       TYPE,DAMP_PARAMETER)

C
C DRIVER FOR DIIS CONVERGENCE ACCELERATION FOR DERIVATIVE T AMPLITUDE
C EQUATIONS
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL SING
C
      DIMENSION SCR(MAXCOR)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/FLAGS/IFLAGS(100)
cYAU - This is the same in DIISLST and should probably be an include file.
cjp
cjp separate data for individual references
cjp probably would be better to treat all amplitudes togerher
cjp as one large vector, but it seems to work like this nicely


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


      integer ld_diis
      parameter (ld_diis=25)
      double precision r(ld_diis,ld_diis,maxref)
      common/diisdat1/r,idimdiis(maxref)
      save  /diisdat1/
cjp
C
      DATA ONE,ONEM/1.D0,-1.D0/
C
      INDXF(I,J,N)=I+(J-1)*N
C
cjp   initialize
      lenvec=idimdiis(iref)
      ndimr=min(icycle,iflags(12))
      idimr=1+mod(icycle-1,iflags(12))
cjp
      TOL=10.D0**(-ICONTL)

C
C DO BOOKKEEPING AND COMPUTE RESIDUAL VECTOR
C
      MAXDM1=ld_diis+2
      I000=1
      I010=I000+MAX(LENVEC,MAXDM1*MAXDM1)
      I020=I010+MAX(LENVEC,MAXDM1*MAXDM1)
      IEND=I020+MAXDM1
      IF(IEND.GE.MAXCOR) CALL INSMEM('DODIIS0',IEND,MAXCOR)
C
C Symmetrize the T2 amplitudes for RHF cases.
C
      IF (IUHF.EQ.0) THEN
         CALL SYMETRIZE_T2(SCR, MAXCOR, 1, LISTT20+2)
         CALL SYMETRIZE_T2(SCR, MAXCOR, 1, LISTT21+2)
      END IF
C
      IOFFT=I000
      IOFFHT=I010
      IF(SING)THEN
       DO 10 ISPIN=1,1+IUHF
        CALL GETLST(SCR(IOFFT) ,1,1,1,IOFF0+ISPIN,LISTT10)
        CALL GETLST(SCR(IOFFHT),1,1,1,IOFF1+ISPIN,LISTT11)
        IOFFT=IOFFT+IRPDPD(IRREPX,8+ISPIN)
        IOFFHT=IOFFHT+IRPDPD(IRREPX,8+ISPIN)
10     CONTINUE
      ENDIF
       DO 11 ISPIN=3,1,IUHF-2
       LISTT=LISTT20-1+ISPIN
       NSIZE=IDSYMSZ(IRREPX,ISYTYP(1,LISTT),ISYTYP(2,LISTT))
       CALL GETALL(SCR(IOFFT) ,NSIZE,IRREPX,LISTT20-1+ISPIN)
       CALL GETALL(SCR(IOFFHT),NSIZE,IRREPX,LISTT21-1+ISPIN)

       IOFFT=IOFFT+NSIZE
       IOFFHT=IOFFHT+NSIZE
11    CONTINUE
C
C SAVE OLD T AS NEW APPROXIMATE SOLUTION VECTOR ON 1,LSTERR
C
      CALL PUTLST(SCR(I000),IDIMR,1,1,1,LSTERR)
      CALL SAXPY (LENVEC,ONEM,SCR(I000),1,SCR(I010),1)

      IF (DAMP_PARAMETER .NE. 1.0D0) THEN

         CALL DSCAL(LENVEC, 1.0D0/DAMP_PARAMETER, SCR(I010), 1)
         RESMITER=SCR(I010-1+ISAMAX(LENVEC,SCR(I010),1))
         CALL DSCAL(LENVEC, DAMP_PARAMETER, SCR(I010), 1)
         
      ElSE

         RESMITER=SCR(I010-1+ISAMAX(LENVEC,SCR(I010),1))

      ENDIF 
C
C CALCULATE ERROR MATRIX UP TO THIS POINT
C
      DO I=1,IDIMR-1
         CALL GETLST(SCR(I000),I,1,1,2,LSTERR)
         X=SDOT(LENVEC,SCR(I000),1,SCR(I010),1)
         R(I,IDIMR,iref)=X
         R(IDIMR,I,iref)=X
      END DO
C
C SAVE T(NEW)-T(OLD) AS CORRESPONDING ERROR VECTOR ON 2,LSTERR
C
      CALL PUTLST(SCR(I010),IDIMR,1,1,2,LSTERR)
      X=SDOT(LENVEC,SCR(I010),1,SCR(I010),1)
      
      R(IDIMR,IDIMR,iref)=X
C
C CALCULATE ERROR MATRIX UP TO THIS POINT (cont.)
C
      DO I=IDIMR+1,NDIMR
         CALL GETLST(SCR(I000),I,1,1,2,LSTERR)
         X=SDOT(LENVEC,SCR(I000),1,SCR(I010),1)
         R(I,IDIMR,iref)=X
         R(IDIMR,I,iref)=X
      END DO
C
C CALCULAtE DIIS EXPANSION COEFFICIENTS
C
cjp
      I=ld_diis
      CALL DODIIS(R(1,1,iref),SCR(I000),SCR(I020),I,NDIMR)
C
      IF(isbwcc.or.IFLAGS(1).GE.10)THEN
       WRITE(6,1000)
       CALL PRVECR(SCR(I020),NDIMR)
      ENDIF
C
C EVALUATE DIIS RESIDUAL
C
      fact = scr(i020-1+idimr)
      do i = 0, lenvec-1
         scr(i000+i) = fact*scr(i010+i)
      end do
cYAU - fix this to load as many vectors as possible and use dgemm
      IOFF=I020-1+NDIMR
      DO I=NDIMR,IDIMR+1,-1
         CALL GETLST(SCR(I010),I,1,1,2,LSTERR)
         CALL SAXPY (LENVEC,SCR(IOFF),SCR(I010),1,SCR(I000),1)
         IOFF=IOFF-1
      END DO
      IOFF=IOFF-1
      DO I=IDIMR-1,1,-1
         CALL GETLST(SCR(I010),I,1,1,2,LSTERR)
         CALL SAXPY (LENVEC,SCR(IOFF),SCR(I010),1,SCR(I000),1)
         IOFF=IOFF-1
      END DO
      RESMDIIS=SCR(ISAMAX(LENVEC,SCR(I000),1))
C
C WRITE OUT RESIDUALS AND CHECK CONVERGENCE
C
      WRITE(6,1001)ICYCLE
      WRITE(6,1002)RESMITER
      WRITE(6,1003)RESMDIIS
      IF(ABS(RESMITER).LT.TOL)THEN
       ICONVG=0
      WRITE(6,1004)ICYCLE
C
C This is commented to avoid problems where DIIS converge
C faster than the amplitudes. Be/STO-3G, Ajith /0308/2000 
C 
C      ELSEIF(ABS(RESMDIIS).LT.TOL)THEN
C      ICONVG=0
C       WRITE(6,1004)ICYCLE
      ELSE
       ICONVG=1
      ENDIF
C SG 7/23/98
      IF (IFLAGS(22) .EQ. 1) THEN
        CALL PUTREC(20, 'JOBARC', 'T2CNVCRT', IINTFP,
     &     MIN(ABS(RESMITER), ABS(RESMDIIS)))
      ENDIF
C
C EXTRAPOLATE, SKIP ONLY IN THE FIRST ITERATION
C
      IF(NDIMR.GE.2) THEN
COLD
C
C PUT DIIS RESIDUAL ON DISK - THIS CORRESPONDS TO DIIS EXTRAPOLANT
C
c       CALL PUTLST(SCR(I000),1,1,1,2,LSTERR)
COLD
C
C USE (DIIS EXTRAPOLANT + JACOBI UPDATE) AS NEXT VECTOR
C
C GENERATE UPDATED GUESS, THE JACOBI UPDATE IS ALREADY ON SCR(I000)
C ADD THE SOLUTION VECTORS WITH APPROPRIATE COEFFICIENTS
C
c       CALL ZERO(SCR(I000),LENVEC)
       IOFF=I020
       DO 501 I=1,NDIMR
        FACT=SCR(IOFF)
        CALL GETLST(SCR(I010),I,1,1,1,LSTERR)
        CALL SAXPY (LENVEC,FACT,SCR(I010),1,SCR(I000),1)
        IOFF=IOFF+1
501    CONTINUE
COLD
c       CALL PUTLST(SCR(I000),1,1,1,1,LSTERR)
C
C GENERATE UPDATED JACOBI GUESS
C
c       CALL GETLST(SCR(I010),1,1,1,2,LSTERR)
c       CALL SAXPY (LENVEC,ONE,SCR(I010),1,SCR(I000),1)
COLD
C
C WRITE IT OVER CURRENT INCREMENT LISTS
C
       IOFFT=I000
       IF(SING)THEN
        DO 110 ISPIN=1,1+IUHF
         CALL PUTLST(SCR(IOFFT),1,1,1,IOFF1+ISPIN,LISTT11)
         IOFFT=IOFFT+IRPDPD(IRREPX,8+ISPIN)
110     CONTINUE
       ENDIF
       DO 111 ISPIN=3,3-2*IUHF,-1
        LISTT=LISTT20-1+ISPIN
        NSIZE=IDSYMSZ(IRREPX,ISYTYP(1,LISTT),ISYTYP(2,LISTT))
        CALL PUTALL(SCR(IOFFT),NSIZE,IRREPX,LISTT21-1+ISPIN)
        IOFFT=IOFFT+NSIZE
111    CONTINUE
C
C Create the AA list from AB list for RHF calcualtions. For DRCCD
C Calculations, we do not need to do this for every iteration. After
C conergence, we can do this once to get the correct energies. 
C
CSSS       IF (IUHF.EQ.0) CALL RCC_ABTOAA(SCR, MAXCOR*IINTFP, 
CSSS     &                                IUHF, LISTT21)
C
      ENDIF
C
      RETURN
1000  FORMAT(T3,'Current DIIS expansion coefficients : ')
1001  FORMAT(T3,'Convergence information after ',I5,' iterations: ')
1002  FORMAT(T3,'Largest element of residual vector : ',E15.8,'.')
1003  FORMAT(T3,'Largest element of DIIS residual   : ',E15.8,'.')
1004  FORMAT(T3,'Amplitude equations converged in ',I5, 
     &          ' iterations.')
      END
