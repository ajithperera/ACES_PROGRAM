











c This routine retrieves NCOLS columns starting at column IFIRST
c of a stored double precision 2-dimensional array.

c OUTPUT
c double IDEST : the destination array
c                NOTE: Internally, this is treated as an integer array
c                      since a counter runs across both IDEST and ICORE.

c INPUT
c int IFIRST : the index of the first column to retrieve
c int NCOLS  : the number of columns to retrieve
c int XCACHE : (OBSOLETE) the cache buffer for the I/O operation
c              = 1; T amplitudes
c              = 2; integrals
c int ILEFT  : the left  (row) address index of the array [MOIO(ILEFT,IRIGHT)]
c int IRIGHT : the right (col) address index of the array [MOIO(ILEFT,IRIGHT)]

      subroutine A2_getlst(iDest,iFirst,nCols,xCache,iLeft,iRight)
      implicit none

c ARGUMENTS
      integer iDest(*), iFirst, nCols, xCache, iLeft, iRight

c EXTERNAL FUNCTIONS
      double precision dnrm2, dTmp
      integer idamax

c INTERNAL VARIABLES
      integer iFileNum, iStat
      integer nRows
      integer iRec, iRecNdx, iTmp
      integer nLeft, nGet, iOff
      integer iPos

c COMMON BLOCKS
c icore.com : begin

c icore(1) is an anchor in memory that allows subroutines to address memory
c allocated with malloc. This system will fail if the main memory is segmented
c or parallel processes are not careful in how they allocate memory.

      integer icore(1)
      common / / icore

c icore.com : end
c lists.com : begin
c These common blocks contain global information about the arrays in storage.
c Elements prepended with "bw" are for storing the file metadata while working
c on multiple references.

c moio  (iGrp,iFam) : the physical record that contains the first element
c                     of the array (iGrp,iFam)
c moiowd(iGrp,iFam) : the integer-word index of the first element
c moiods(iGrp,iFam) : the number of columns in the array
c moiosz(iGrp,iFam) : the number of rows    in the array
c moiofl(iGrp,iFam) : the external file unit that contains the array

      integer*8 MAX_IO_GRPS
      integer*8 MAX_IO_FAMS
      integer*8 FAMS_PER_FILE
      integer*8 FIRST_IO_LUN
      integer*8 MAX_IO_LUNS
      integer*8 maxref 

      parameter(MAX_IO_GRPS=10)
      parameter(MAX_IO_FAMS=500)
      parameter(FAMS_PER_FILE=100)
      parameter(FIRST_IO_LUN=50)
      parameter(MAX_IO_LUNS=5)
      parameter(MAXREF=1)

      integer*8        moio  (MAX_IO_GRPS,MAX_IO_FAMS),
     &                 moiowd(MAX_IO_GRPS,MAX_IO_FAMS),
     &                 moiosz(MAX_IO_GRPS,MAX_IO_FAMS),
     &                 moiods(MAX_IO_GRPS,MAX_IO_FAMS),
     &                 moiofl(MAX_IO_GRPS,MAX_IO_FAMS),
     &               bwmoio  (MAX_IO_GRPS,MAX_IO_FAMS,maxref),
     &               bwmoiowd(MAX_IO_GRPS,MAX_IO_FAMS,maxref),
     &               bwmoiosz(MAX_IO_GRPS,MAX_IO_FAMS,maxref),
     &               bwmoiods(MAX_IO_GRPS,MAX_IO_FAMS,maxref),
     &               bwmoiofl(MAX_IO_GRPS,MAX_IO_FAMS,maxref)
      common /lists/   moio,   moiowd,   moiosz,   moiods,   moiofl,
     &               bwmoio, bwmoiowd, bwmoiosz, bwmoiods, bwmoiofl
      save   /lists/

c moiomxsz(iGrp,iFam) : the original length of a one-dimensional array
c                       (This is shameful. Arrays should not be re-dimensioned
c                        at will during a job.)

      integer*8           moiomxsz(MAX_IO_GRPS,MAX_IO_FAMS),
     &                    bwmoiomxsz(MAX_IO_GRPS,MAX_IO_FAMS,maxref)
      common /lists_mxsz/   moiomxsz,
     &                    bwmoiomxsz
      save   /lists_mxsz/

c pRec(i)    : the index of the physical record in file i containing free space
c              (i is the internal unit number of the storage file.)
c iIntOff(i) : the integer offset from the beginning of the physical record
c              needed to address the free space

      integer*8          pRec   (MAX_IO_LUNS),
     &                   iIntOff(MAX_IO_LUNS),
     &                 bwpRec   (MAX_IO_LUNS,maxref),
     &                 bwiIntOff(MAX_IO_LUNS,maxref)
      common /io_ptrs/   pRec,   iIntOff,
     &                 bwpRec, bwiIntOff
      save   /io_ptrs/

c bIOUp  : a flag for bombing in get/putlst if aces_io_init has not been called
c bIOMod : a flag for updating the records in aces_io_fin

      logical           bIOUp, bIOMod
      common /io_flags/ bIOUp, bIOMod
      save   /io_flags/

c lists.com : end
c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer*8       iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end
c filspc.com : begin

c This common block contains the dimensions of the physical records used by the
c MOIO storage files.

c iprcln : the    byte-length of a physical record
c iprcwd : the integer-length of a physical record
c iprcrl : the    recl-length of a physical record

      integer*8       iprcln, iprcwd, iprcrl
      common /filspc/ iprcln, iprcwd, iprcrl
      save   /filspc/

c filspc.com : end
c auxcache.com : begin

c The auxiliary cache is a programmer-controlled list cache. The dimensions are
c the same as those of the MOIO arrays. The programmer may load lists into icore
c memory and set quikget(?,?) to their icore addresses. When getlst and putlst
c operate on ANY list, quikget is checked to see if the list lives in icore.
c If so, the operation is performed on the in-core data instead of hitting the
c storage file(s).

c WARNING
c    There is no automatic updating of storage files with data in the auxiliary
c cache. If the memory-resident data is altered and must be stored to disk, then
c the quikget value must be destroyed (zeroed) and the actual location of the
c data must then be passed to putlst. The aces_auxcache_flush routine
c systematically stores the quikget values, calls putlst, and restores the
c quikget values.

      integer*8 quikget(MAX_IO_GRPS,MAX_IO_FAMS)
      common /auxcache/ quikget
      save   /auxcache/

c auxcache.com : end

c ----------------------------------------------------------------------

      iTmp = 0
c   o assert I/O subsystem is up
      if (.not.bIOUp) then
         print *, '@GETLST: Assertion failed.'
         print *, '   bIOUp = ',bIOUp
         iTmp = 1
      end if
c   o assert iFirst > 0 and nCols >= 0
      if ((iFirst.lt.1).or.(nCols.lt.0)) then
         print *, '@GETLST: Assertion failed.'
         print *, '   iFirst = ',iFirst
         print *, '   nCols  = ',nCols
         iTmp = 1
      end if
c   o assert iLeft and iRight are properly bound
      if ((iLeft .lt.1).or.(MAX_IO_GRPS.lt.iLeft ).or.
     &    (iRight.lt.1).or.(MAX_IO_FAMS.lt.iRight)    ) then
         print *, '@GETLST: Assertion failed.'
         print *, '   iRight = ',iRight
         print *, '   iLeft  = ',iLeft
         iTmp = 1
      end if
c   o assert the list was touched
      if (moio(iLeft,iRight).lt.1) then
         print *, '@GETLST: Assertion failed.'
         print *, '   List (',iLeft,',',iRight,') does not exist.'
         iTmp = 1
      end if
c   o temporary trap to gauge list reads
c
c      if (iRight .EQ. 17 .or. iRight .eq. 18 .or.
c     &    iRight .EQ. 19 .or. iRight .eq. 20 .or.
c     &    iRight .EQ. 21 .or. iRight .eq. 22) then
c           Print *, 'Reading Redundent lists'
c           iTmp = 1
c      endif 
      if (iTmp.ne.0) call aces_exit(iTmp)

c ----------------------------------------------------------------------

c   o look up the length of one column
      nRows = moiosz(iLeft,iRight)

      if ((nRows.lt.1).or.(nCols.lt.1)) return

c   o make sure the column range requested is properly bound
      iTmp = moiods(iLeft,iRight)
      if ((iFirst.lt.1).or.(iTmp.lt.iFirst).or.
     &    (iTmp.lt.(iFirst-1+nCols))           ) then
         print *, '@GETLST: Error reading list (',iLeft,',',iRight,')'
         print *, '         cols available = ',iTmp
         print *, '         cols requested = ',nCols
         print *, '         range start    = ',iFirst
         call aces_exit(1)
      end if

c   o do a fast in-core fetch from the auxiliary cache
      if (quikget(iLeft,iRight).ne.0) then
         iPos = quikget(iLeft,iRight) + nRows*(iFirst-1)*iintfp
         call xcopy(nRows*nCols,icore(iPos),1,iDest,1)
         return
      end if

c   o determine which file this array lives in
      iFileNum = moiofl(iLeft,iRight)


c   o find the first record and integer index that point to the first element
      iRecNdx = moiowd(iLeft,iRight) + nRows*(iFirst-1)*iintfp
      iTmp    = (iRecNdx-1)/iprcwd
      iRec    = moio(iLeft,iRight) + iTmp
      iRecNdx = iRecNdx            - iTmp*iprcwd

c   o do the first partial record and then loop over whole records
c     which contain the remaining data
      nLeft = nRows*nCols*iintfp
      nGet  = min(nLeft,iprcwd+1-iRecNdx)
      call getlst_io(iDest,iFileNum,iRec,iRecNdx,nGet)
      iOff  = 1     + nGet
      nLeft = nLeft - nGet
      do while (nLeft.ne.0)
         nGet = min(nLeft,iprcwd)
         iRec = iRec + 1
         call getlst_io(iDest(iOff),iFileNum,iRec,1,nGet)
         iOff  = iOff  + nGet
         nLeft = nLeft - nGet
      end do



      return

c   o I/O error
 666  print *, '@GETLST: I/O error'
      print *, '         list = [',iLeft,',',iRight,']'
      print '(/)'
      call aces_io_error('GETLST',iFileNum,iStat)

c     end subroutine getlst
      end

