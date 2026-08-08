

















      Subroutine Getmem_a2(Icore,iCoreNdx,iCoreDim,iMem,bAllocMem)

      implicit none

c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end


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



      integer              xFlag
      common /aces_reflag/ xFlag
      save   /aces_reflag/
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
 
c   o allocate memory and initialize the I/O subsystem
      if (bAllocMem) then

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
               xFlag = 1
            else
c               print *, '@ACES_INIT: Request for ',iMem,
c     &                  ' integers of memory failed.'
               print *, '@ACES_INIT: unable to allocate at least ',
     &                  iMemMin,' integers of memory'
               call aces_exit(1)
            end if
         end if

c      o initialize the I/O subsystem ('T' creates and initializes the
cache)
CSSS         call aces_io_init(iCore,i0,iMem,.true.)

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

      Return
      End
