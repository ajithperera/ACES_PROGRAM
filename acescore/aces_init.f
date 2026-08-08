











c This routine initializes the ACES environment for use by any program.
c It or aces_init_rte should be called before any other ACES routine.

c INPUT
c logical BALLOCMEM : do (not) allocate memory according to IFLAGS(h_IFLAGS_mem)
c                     NOTE: If memory is allocated, then aces_init will
c                           initialize the I/O subsystem and the cache
c                           subsystem.

c INPUT/OUTPUT
c int ICORE(*) : a reference integer (address) used to anchor a malloc'ed heap,
c                the value of which is set to the byte address of the heap
c                (not valid on 32-bit builds with 64-bit pointers)

c OUTPUT
c int ICORENDX : the integer index that addresses the malloc'ed heap
c int ICOREDIM : the number of addressable integers at ICORE(ICORENDX)
c int IUHF : the spin state of the reference wavefunction
c            = 0; spin-  restricted (RHF)
c            = 1; spin-unrestricted (UHF)








c#define _DEBUG_ACES_INIT

      subroutine aces_init(iCore,iCoreNdx,iCoreDim,iUHF,bAllocMem)
      implicit none

c ARGUMENTS
      integer iCore(*), iCoreNdx, iCoreDim, iUHF
      logical bAllocMem

c EXTERNAL FUNCTIONS
      INTEGER c_adr
      external c_adr

c PARAMETERS
      integer iMemMin, iMemInc
      parameter (iMemMin=1*1024*1024,iMemInc=1*1024*1024)

c INTERNAL VARIABLES
      double precision dTmp
      integer i0, iMem, iTmp
c i0 must survive across reentrant calls: on a reentrant call (xFlag.ne.0)
c the malloc block below that sets it is skipped, but iCoreNdx=i0 further
c down still runs unconditionally -- without SAVE, i0 is undefined on
c that path and corrupts the caller's ICORE indexing (istart.com's I0).
      save i0

c TIMING VARIABLES
      integer stime_sec, stime_usec
      integer utime_sec, utime_usec
      integer rtime_sec, rtime_usec

c COMMON BLOCKS
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end


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



c aces_time.com : begin

c These six times hold the timing data from aces_init (*_in) to
c aces_fin (*_out). utime and stime count the total number of user
c and system seconds since the start of the process. rtime counts
c the total number of real seconds since 1 January 1970. ame_timed
c is a logical flag set in aces_init telling aces_fin to print out
c a timing summary.

      external aces_bd_aces_time

      double precision   ame_utime_in,  ame_stime_in,  ame_rtime_in,
     &                   ame_utime_out, ame_stime_out, ame_rtime_out
      logical            ame_timed
      common /aces_time/ ame_utime_in,  ame_stime_in,  ame_rtime_in,
     &                   ame_utime_out, ame_stime_out, ame_rtime_out,
     &                   ame_timed
      save   /aces_time/

c aces_time.com : end
c jobarc.com : begin

c This data tracks the contents of the JOBARC file. 'physical' records refer
c to direct I/O while 'logical' records refer to the ACES archive elements.





      external aces_bd_jobarc

c marker(i) : the name of logical record i
c rloc(i)   : the integer index in JOBARC that starts logical record i
c rsize(i)  : the integer-length of logical record i
c nrecs  : the number of physical records in the JOBARC file
c irecwd : the integer-length of a physical record
c irecln : the    recl-length of a physical record

      character*8     marker(1000)
      integer         rloc  (1000),
     &                rsize (1000),
     &                nrecs, irecwd, irecln
      common /jobarc/ marker,
     &                rloc,
     &                rsize,
     &                nrecs, irecwd, irecln
      save   /jobarc/

c bJAUp  : a flag for bombing in get/putrec if aces_ja_init has not been called
c bJAMod : a flag for updating JAINDX in aces_ja_fin

      logical           bJAUp, bJAMod
      common /ja_flags/ bJAUp, bJAMod
      save   /ja_flags/

c jobarc.com : end
c lists.com : begin

c These common blocks contain global information about the arrays in storage.
c Elements prepended with "bw" are for storing the file metadata while working
c on multiple references.






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



      external aces_bd_lists

c moio  (iGrp,iFam) : the physical record that contains the first element
c                     of the array (iGrp,iFam)
c moiowd(iGrp,iFam) : the integer-word index of the first element
c moiods(iGrp,iFam) : the number of columns in the array
c moiosz(iGrp,iFam) : the number of rows    in the array
c moiofl(iGrp,iFam) : the external file unit that contains the array

      integer          moio  (10,500),
     &                 moiowd(10,500),
     &                 moiosz(10,500),
     &                 moiods(10,500),
     &                 moiofl(10,500),
     &               bwmoio  (10,500,maxref),
     &               bwmoiowd(10,500,maxref),
     &               bwmoiosz(10,500,maxref),
     &               bwmoiods(10,500,maxref),
     &               bwmoiofl(10,500,maxref)
      common /lists/   moio,   moiowd,   moiosz,   moiods,   moiofl,
     &               bwmoio, bwmoiowd, bwmoiosz, bwmoiods, bwmoiofl
      save   /lists/

c moiomxsz(iGrp,iFam) : the original length of a one-dimensional array
c                       (This is shameful. Arrays should not be re-dimensioned
c                        at will during a job.)

      integer               moiomxsz(10,500),
     &                    bwmoiomxsz(10,500,maxref)
      common /lists_mxsz/   moiomxsz,
     &                    bwmoiomxsz
      save   /lists_mxsz/

c pRec(i)    : the index of the physical record in file i containing free space
c              (i is the internal unit number of the storage file.)
c iIntOff(i) : the integer offset from the beginning of the physical record
c              needed to address the free space

      integer            pRec   (5),
     &                   iIntOff(5),
     &                 bwpRec   (5,maxref),
     &                 bwiIntOff(5,maxref)
      common /io_ptrs/   pRec,   iIntOff,
     &                 bwpRec, bwiIntOff
      save   /io_ptrs/

c bIOUp  : a flag for bombing in get/putlst if aces_io_init has not been called
c bIOMod : a flag for updating the records in aces_io_fin

      logical           bIOUp, bIOMod
      common /io_flags/ bIOUp, bIOMod
      save   /io_flags/

c lists.com : end
c BWCC needs aces_init to be reentrant, so an initialized common value must
c exist. The blockdata routine is appended to the end of this file.
      external aces_bd_aces_reflag
      integer              xFlag
      common /aces_reflag/ xFlag
      save   /aces_reflag/

c ----------------------------------------------------------------------


c   o initialize timing information
      call c_rutimes(utime_sec,utime_usec,stime_sec,stime_usec)
      call c_gtod   (rtime_sec,rtime_usec)
      ame_stime_in = (1.d-6 * stime_usec) + stime_sec
      ame_utime_in = (1.d-6 * utime_usec) + utime_sec
      ame_rtime_in = (1.d-6 * rtime_usec) + rtime_sec
      ame_timed = .true.


c   o setup stdout buffering to see output immediately
      call bufferlines

c   o initialize the runtime environment (held in common blocks)
      call aces_init_rte

c   o gather parallel statistics (needed by gfname in ja_init)
      call aces_com_parallel_aces

c   o initialize the job archive subsystem (getrec/putrec)
c     In pyaces's single-process model, aces_init can be called again by
c     legacy routines (e.g. CRAPSI) after the job archive is already up
c     from an earlier phase in the SAME process (separate processes each
c     called this exactly once in the traditional architecture). Cycle it
c     closed-then-reopened instead of hard-erroring; aces_ja_fin flushes
c     the in-memory record index to disk first, so nothing already
c     written is lost, and the following aces_ja_init re-reads it fresh.
      if (bJAUp) call aces_ja_fin
      call aces_ja_init
c???      call getrec(1,'JOBARC','IENDSTAT',1,iStat)

c   o load the ACES State Variables
      call getrec(1,'JOBARC','IFLAGS', 100,iflags)
      call getrec(1,'JOBARC','IFLAGS2',500,iflags2)

c   o allocate memory and initialize the I/O subsystem
      if (bAllocMem) then

c      o load the requested core size
         iMem = iflags(36)
         WRITE(6,*) '@ACES_INIT-DEBUG: iflags(36)=', iflags(36),
     &              ' iMem=', iMem

c      o allocate core memory
         if (xFlag.eq.0) then
            iCore(1) = 0
            do while ((iCore(1).eq.0).and.(iMem.gt.iMemMin))
               call aces_malloc(iMem,iCore,i0)
               if (iCore(1).eq.0) iMem = iMem - iMemInc
            end do
            if (iMem.lt.iflags(36)) then
               print *, '@ACES_INIT: MEMORY WARNING!'
               print *, '            requested ',iflags(36),' integers'
               print *, '            allocated ',iMem,' integers'
            end if
            if (iCore(1).ne.0) then
               WRITE(6,*) '@ACES_INIT-DEBUG: base addr=',
     &                    c_adr(iCore(i0)), ' end addr=',
     &                    c_adr(iCore(i0+iMem))-1, ' i0=', i0,
     &                    ' iMem=', iMem
               xFlag = 1
            else
c               print *, '@ACES_INIT: Request for ',iMem,
c     &                  ' integers of memory failed.'
               print *, '@ACES_INIT: unable to allocate at least ',
     &                  iMemMin,' integers of memory'
               call aces_exit(1)
            end if
         end if

c      o initialize the I/O subsystem ('T' creates and initializes the cache)
c        Same reentrancy issue as the job archive above: cycle it
c        closed-then-reopened if a legacy routine calls aces_init again
c        later in the same process; aces_io_fin flushes pending list
c        metadata to JOBARC first, so nothing already written is lost.
         if (bIOUp) call aces_io_fin
         call aces_io_init(iCore,i0,iMem,.true.)

c      o transfer the iCore statistics
         iCoreNdx = i0
         iCoreDim = iMem

c     else if (.not.bAllocMem)
      else

         iCore(1) = 0
         iCoreNdx = 1
         iCoreDim = 1

c     end if (bAllocMem)
      end if

c   o assign iUHF
      if (iflags(11).eq.0) then
         iUHF = 0
      else
         iUHF = 1
      end if

c   o initialize the chemical system (held in common blocks)
      call aces_init_chemsys


      return
c     end subroutine aces_init
      end

c ----------------------------------------------------------------------
      blockdata aces_bd_aces_reflag
      integer              xFlag
      common /aces_reflag/ xFlag
      save   /aces_reflag/
      data                 xFlag /0/
      end
c ----------------------------------------------------------------------

