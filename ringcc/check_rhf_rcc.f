










      Subroutine check_rhf_rcc(Work,Maxcor,Iuhf)

      Implicit Double precision (A-H,O-Z)
      Dimension Work(Maxcor)
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







c This common block contains the IFLAGS and IFLAGS2 arrays for JODA ROUTINES
c ONLY! The reason is that it contains both arrays back-to-back. If the
c preprocessor define MONSTER_FLAGS is set, then the arrays are compressed
c into one large (currently) 600 element long array; otherwise, they are
c split into IFLAGS(100) and IFLAGS2(500).

c iflags(100)  ASVs reserved for Stanton, Gauss, and Co.
c              (Our code is already irrevocably split, why bother anymore?)
c iflags2(500) ASVs for everyone else

      integer        iflags(100), iflags2(500)
      common /flags/ iflags,      iflags2
      save   /flags/






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

C Form spin-adapted integral list and write rhem to lists 117 and 118.

      Imode  = 0
      Irrepx = 1
      CALL INIPCK(1,19,21,197,IMODE,0,1)
      CALL INIPCK(1,19,21,198,IMODE,0,1)
C
      List_aa  = 61
      List_ab  = 63 
      Nsize_aa = Idsymsz(Irrepx,Isytyp(1,LIst_aa),Isytyp(2,List_aa))
      Nsize_ab = Idsymsz(Irrepx,Isytyp(1,LIst_ab),Isytyp(2,List_ab))
      If (Nsize_aa .NE. Nsize_ab) Then
         Write(6,"(a)") "Different AAAA and BBBB blocks"
         call Errex
      Endif 

      I000 = 1
      I010 = I000 + Nsize_aa
      I020 = I010 + Nsize_ab
      I030 = I020 + Nsize_ab

      Call Getall(Work(I000),Nsize_AA,Irrepx,List_AA)
      Call Dcopy(Nsize_aa,Work(I000),1,Work(I010),1)
      Call Getall(Work(I020),Nsize_AA,Irrepx,List_AB)

c Form T2(AAAA) + T2(ABAB)
     
       Call Daxpy(Nsize_aa,1.0D0,Work(I020),1,Work(I000),1)

      call checksum("2J-K:",Work(I000),NSIZE_ab)
c Form T2(AAAA) - T2(ABAB) 

       Call Daxpy(Nsize_aa,-1.0D0,Work(I020),1,Work(I010),1)
CSSS       Call Dcopy(Nsize_aa,Work(I020),1,Work(I010),1)
      call checksum("-K  :",Work(I010),NSIZE_ab)
       Nsize_ajbi = Idsymsz(Irrepx,Isytyp(1,21),Isytyp(2,21))
       Nsize_abij = Idsymsz(Irrepx,Isytyp(1,16),Isytyp(2,16))

      I030  = I020 + NSIZE_AJBI 
      I040  = I030 + NSIZE_ABIJ 
      I050  = I040 + NSIZE_ABIJ
      IEND  = I050

      If ((Nsize_aa .NE. Nsize_ab)     .or. 
     +    (NSIZE_AJBI .ne. NSIZE_ABIJ) .or.
     +    (Nsize_aa .NE. NSIZE_AJBI))   Then
         Write(6,"(a)") "Different AAAA and BBBB blocks"
         call Errex
      Endif 
      If (Iend .GE. MAXcor) Call Insmem("check_rhf_rcc",Iend,
     +                                   Maxcor)

C Read <Aj|bi> integrals

      CALL GETALL(Work(I020),NSIZE_AJBI,IRREPX,21)

C <Aj|bi> -> <Ab|ji>

      CALL SSTGEN(Work(I020),Work(I030),NSIZE_ABIJ,VRT(1,1),
     +            POP(1,2),VRT(1,2),POP(1,1),Work(I040),IRREPX,
     +           "1324")

C These are exchange integrals of <Ab|Ij> (RHF only)

C Form the triplet part of the spin-adapted <Ab|ji> by -<Ab|ji>

      CALL DSCAL(NSIZE_ABIJ,-1.0D0,Work(I030),1)
      CALL PUTALL(Work(I030),NSIZE_ABIJ,IRREPX,198)

      call checksum("List-198:",Work(I030),NSIZE_ABIJ)

c Form the singlet part of the spin-adapted <Ab|ij> by 2<Ab|Ij>-<Ab|jI>

      CALL GETALL(Work(I020),NSIZE_ABIJ,IRREPX,16)
      CALL DSCAL(NSIZE_ABIJ,2.0D0,Work(I020),1)
      CALL DAXPY(NSIZE_ABIJ,1.0D0,Work(I030),1,Work(I020),1)
      CALL PUTALL(Work(I020),NSIZE_ABIJ,IRREPX,197)
      call checksum("List-197:",Work(I020),NSIZE_ABIJ)

c Compute the triplet and singlet spin-adapted energies.

      Call Getall(Work(I020),Nsize_abij,Irrepx,197)
      E_s = Ddot(Nsize_abij,Work(I020),1,Work(I000),1)

      Call Getall(Work(I020),Nsize_abij,Irrepx,198)
      E_t = Ddot(Nsize_abij,Work(I020),1,Work(I010),1)
     
      Write(6,"(a,F15.9)") "Singlet-spin adapted energy: ", E_s
      Write(6,"(a,F15.9)") "Triplet-spin adapted energy: ", E_t

      E_total = 0.75D0*E_t + 0.25D0*E_s
      
      Write(6,"(a,F15.9)") "Total energy: ", E_total

      Return
      End
