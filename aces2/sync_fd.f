











c This routine gathers all the data needed for parallel finite differences.
c The goto routine for what should be done is symcor/upd_fd.F





































































































































































































      subroutine sync_fd
      implicit none

!
!
!  (C) 1993 by Argonne National Laboratory and Mississipi State University.
!      All rights reserved.  See COPYRIGHT in top-level directory.
!
!
! user include file for MPI programs, with no dependencies
!
! It really is not possible to make a perfect include file that can
! be used by both F77 and F90 compilers, but this is close.  We have removed
! continuation lines (allows free form input in F90); systems whose
! Fortran compilers support ! instead of just C or * for comments can
! globally replace a C in the first column with !; the resulting file
! should work for both Fortran 77 and Fortran 90.
!
! If your Fortran compiler supports ! for comments, you can run this
! through sed with
!     sed -e 's/^C/\!/g'
!
! We have also removed the use of contractions (involving the single quote)
! character because some users use .F instead of .f files (to invoke the
! cpp preprocessor) and further, their preprocessor is determined to find
! matching single quote pairs (and probably double quotes; given the
! different rules in C and Fortran, this sounds like a disaster).  Rather than
! take the position that the poor users should get a better system, we
! have removed the text that caused problems.  Of course, the users SHOULD
! get a better system...
!
! return codes
      INTEGER MPI_SUCCESS,MPI_ERR_BUFFER,MPI_ERR_COUNT,MPI_ERR_TYPE
      INTEGER MPI_ERR_TAG,MPI_ERR_COMM,MPI_ERR_RANK,MPI_ERR_ROOT
      INTEGER MPI_ERR_GROUP
      INTEGER MPI_ERR_OP,MPI_ERR_TOPOLOGY,MPI_ERR_DIMS,MPI_ERR_ARG
      INTEGER MPI_ERR_UNKNOWN,MPI_ERR_TRUNCATE,MPI_ERR_OTHER
      INTEGER MPI_ERR_INTERN,MPI_ERR_IN_STATUS,MPI_ERR_PENDING
      INTEGER MPI_ERR_REQUEST, MPI_ERR_LASTCODE
      PARAMETER (MPI_SUCCESS=0,MPI_ERR_BUFFER=1,MPI_ERR_COUNT=2)
      PARAMETER (MPI_ERR_TYPE=3,MPI_ERR_TAG=4,MPI_ERR_COMM=5)
      PARAMETER (MPI_ERR_RANK=6,MPI_ERR_ROOT=7,MPI_ERR_GROUP=8)
      PARAMETER (MPI_ERR_OP=9,MPI_ERR_TOPOLOGY=10,MPI_ERR_DIMS=11)
      PARAMETER (MPI_ERR_ARG=12,MPI_ERR_UNKNOWN=13)
      PARAMETER (MPI_ERR_TRUNCATE=14,MPI_ERR_OTHER=15)
      PARAMETER (MPI_ERR_INTERN=16,MPI_ERR_IN_STATUS=17)
      PARAMETER (MPI_ERR_PENDING=18,MPI_ERR_REQUEST=19)
      PARAMETER (MPI_ERR_LASTCODE=1073741823)
!
      INTEGER MPI_UNDEFINED
      parameter (MPI_UNDEFINED = (-32766))
!
      INTEGER MPI_GRAPH, MPI_CART
      PARAMETER (MPI_GRAPH = 1, MPI_CART = 2)
      INTEGER  MPI_PROC_NULL
      PARAMETER ( MPI_PROC_NULL = (-1) )
!
      INTEGER MPI_BSEND_OVERHEAD
      PARAMETER ( MPI_BSEND_OVERHEAD = 512 )

      INTEGER MPI_SOURCE, MPI_TAG, MPI_ERROR
      PARAMETER(MPI_SOURCE=2, MPI_TAG=3, MPI_ERROR=4)
      INTEGER MPI_STATUS_SIZE
      PARAMETER (MPI_STATUS_SIZE=4)
      INTEGER MPI_MAX_PROCESSOR_NAME, MPI_MAX_ERROR_STRING
      PARAMETER (MPI_MAX_PROCESSOR_NAME=256)
      PARAMETER (MPI_MAX_ERROR_STRING=512)
      INTEGER MPI_MAX_NAME_STRING
      PARAMETER (MPI_MAX_NAME_STRING=63)
!
      INTEGER MPI_COMM_NULL
      PARAMETER (MPI_COMM_NULL=0)
!
      INTEGER MPI_DATATYPE_NULL
      PARAMETER (MPI_DATATYPE_NULL = 0)

      INTEGER MPI_ERRHANDLER_NULL
      PARAMETER (MPI_ERRHANDLER_NULL = 0)

      INTEGER MPI_GROUP_NULL
      PARAMETER (MPI_GROUP_NULL = 0)

      INTEGER MPI_KEYVAL_INVALID
      PARAMETER (MPI_KEYVAL_INVALID = 0)

      INTEGER MPI_REQUEST_NULL
      PARAMETER (MPI_REQUEST_NULL = 0)
!
      INTEGER MPI_IDENT, MPI_CONGRUENT, MPI_SIMILAR, MPI_UNEQUAL
      PARAMETER (MPI_IDENT=0, MPI_CONGRUENT=1, MPI_SIMILAR=2)
      PARAMETER (MPI_UNEQUAL=3)
!
!     MPI_BOTTOM needs to be a known address; here we put it at the
!     beginning of the common block.  The point-to-point and collective
!     routines know about MPI_BOTTOM, but MPI_TYPE_STRUCT as yet does not.
!
!     MPI_STATUS_IGNORE and MPI_STATUSES_IGNORE are similar objects
!     Until the underlying MPI library implements the C version of these
!     (a null pointer), these are declared as arrays of MPI_STATUS_SIZE
!
!     The types MPI_INTEGER1,2,4 and MPI_REAL4,8 are OPTIONAL.
!     Their values are zero if they are not available.  Note that
!     using these reduces the portability of code (though may enhance
!     portability between Crays and other systems)
!
      INTEGER MPI_TAG_UB, MPI_HOST, MPI_IO
      INTEGER MPI_BOTTOM
      INTEGER MPI_STATUS_IGNORE(MPI_STATUS_SIZE)
      INTEGER MPI_STATUSES_IGNORE(MPI_STATUS_SIZE)
      INTEGER MPI_INTEGER, MPI_REAL, MPI_DOUBLE_PRECISION
      INTEGER MPI_COMPLEX, MPI_DOUBLE_COMPLEX,MPI_LOGICAL
      INTEGER MPI_CHARACTER, MPI_BYTE, MPI_2INTEGER, MPI_2REAL
      INTEGER MPI_2DOUBLE_PRECISION, MPI_2COMPLEX, MPI_2DOUBLE_COMPLEX
      INTEGER MPI_UB, MPI_LB
      INTEGER MPI_PACKED, MPI_WTIME_IS_GLOBAL
      INTEGER MPI_COMM_WORLD, MPI_COMM_SELF, MPI_GROUP_EMPTY
      INTEGER MPI_SUM, MPI_MAX, MPI_MIN, MPI_PROD, MPI_LAND, MPI_BAND
      INTEGER MPI_LOR, MPI_BOR, MPI_LXOR, MPI_BXOR, MPI_MINLOC
      INTEGER MPI_MAXLOC
      INTEGER MPI_OP_NULL
      INTEGER MPI_ERRORS_ARE_FATAL, MPI_ERRORS_RETURN
!
      PARAMETER (MPI_ERRORS_ARE_FATAL=119)
      PARAMETER (MPI_ERRORS_RETURN=120)
!
      PARAMETER (MPI_COMPLEX=23,MPI_DOUBLE_COMPLEX=24,MPI_LOGICAL=25)
      PARAMETER (MPI_REAL=26,MPI_DOUBLE_PRECISION=27,MPI_INTEGER=28)
      PARAMETER (MPI_2INTEGER=29,MPI_2COMPLEX=30,MPI_2DOUBLE_COMPLEX=31)
      PARAMETER (MPI_2REAL=32,MPI_2DOUBLE_PRECISION=33,MPI_CHARACTER=1)
      PARAMETER (MPI_BYTE=3,MPI_UB=16,MPI_LB=15,MPI_PACKED=14)

      INTEGER MPI_ORDER_C, MPI_ORDER_FORTRAN
      PARAMETER (MPI_ORDER_C=56, MPI_ORDER_FORTRAN=57)
      INTEGER MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_CYCLIC
      INTEGER MPI_DISTRIBUTE_NONE, MPI_DISTRIBUTE_DFLT_DARG
      PARAMETER (MPI_DISTRIBUTE_BLOCK=121, MPI_DISTRIBUTE_CYCLIC=122)
      PARAMETER (MPI_DISTRIBUTE_NONE=123)
      PARAMETER (MPI_DISTRIBUTE_DFLT_DARG=-49767)
      INTEGER MPI_MAX_INFO_KEY, MPI_MAX_INFO_VAL
      PARAMETER (MPI_MAX_INFO_KEY=255, MPI_MAX_INFO_VAL=1024)
      INTEGER MPI_INFO_NULL
      PARAMETER (MPI_INFO_NULL=0)

!
! Optional Fortran Types.  Configure attempts to determine these.
!
      INTEGER MPI_INTEGER1, MPI_INTEGER2, MPI_INTEGER4, MPI_INTEGER8
      INTEGER MPI_INTEGER16
      INTEGER MPI_REAL4, MPI_REAL8, MPI_REAL16
      INTEGER MPI_COMPLEX8, MPI_COMPLEX16, MPI_COMPLEX32
      PARAMETER (MPI_INTEGER1=0,MPI_INTEGER2=4)
      PARAMETER (MPI_INTEGER4=6)
      PARAMETER (MPI_INTEGER8=13)
      PARAMETER (MPI_INTEGER16=0)
      PARAMETER (MPI_REAL4=10)
      PARAMETER (MPI_REAL8=11)
      PARAMETER (MPI_REAL16=12)
      PARAMETER (MPI_COMPLEX8=23)
      PARAMETER (MPI_COMPLEX16=24)
      PARAMETER (MPI_COMPLEX32=0)

      COMMON /MPIPRIV/ MPI_BOTTOM,MPI_STATUS_IGNORE,MPI_STATUSES_IGNORE
!
!     Without this save, some Fortran implementations may make the common
!     dynamic!
!
!     For a Fortran90 module, we might replace /MPIPRIV/ with a simple
!     SAVE MPI_BOTTOM
!
      SAVE /MPIPRIV/

      PARAMETER (MPI_MAX=100,MPI_MIN=101,MPI_SUM=102,MPI_PROD=103)
      PARAMETER (MPI_LAND=104,MPI_BAND=105,MPI_LOR=106,MPI_BOR=107)
      PARAMETER (MPI_LXOR=108,MPI_BXOR=109,MPI_MINLOC=110)
      PARAMETER (MPI_MAXLOC=111, MPI_OP_NULL=0)
!
      PARAMETER (MPI_GROUP_EMPTY=90,MPI_COMM_WORLD=91,MPI_COMM_SELF=92)
      PARAMETER (MPI_TAG_UB=80,MPI_HOST=82,MPI_IO=84)
      PARAMETER (MPI_WTIME_IS_GLOBAL=86)
!
      INTEGER MPI_ANY_SOURCE
      PARAMETER (MPI_ANY_SOURCE = (-2))
      INTEGER MPI_ANY_TAG
      PARAMETER (MPI_ANY_TAG = (-1))
!
      INTEGER MPI_VERSION, MPI_SUBVERSION
      PARAMETER (MPI_VERSION    = 1, MPI_SUBVERSION = 2)
!
!     There are additional MPI-2 constants
      INTEGER MPI_ADDRESS_KIND, MPI_OFFSET_KIND
      PARAMETER (MPI_ADDRESS_KIND=4)
      PARAMETER (MPI_OFFSET_KIND=8)
!
!     All other MPI routines are subroutines
!     This may cause some Fortran compilers to complain about defined and
!     not used.  Such compilers should be improved.
!
!     Some Fortran compilers will not link programs that contain
!     external statements to routines that are not provided, even if
!     the routine is never called.  Remove PMPI_WTIME and PMPI_WTICK
!     if you have trouble with them.
!
      DOUBLE PRECISION MPI_WTIME, MPI_WTICK,PMPI_WTIME,PMPI_WTICK
      EXTERNAL MPI_WTIME, MPI_WTICK,PMPI_WTIME,PMPI_WTICK
!
!     The attribute copy/delete subroutines are symbols that can be passed
!     to MPI routines
!
!      EXTERNAL MPI_NULL_COPY_FN, MPI_NULL_DELETE_FN, MPI_DUP_FN

c INTERNAL PARAMETERS
      integer iMemInc
      parameter (iMemInc=1048576)

c INTERNAL VARIABLES
      integer iCore(2), i0, iCrSiz, iMemMin
      double precision dDipXYZ(3), dPolXYZ(3,3)
      integer nDbls, nAtom, nSize, nPoint, iOff, iLast, iTmp
      integer lSend, lRecv, lFree, lGeom, lGrdXYZ, lGrdInt, lIMAP
      character*4 szType
      integer iCount, iRoot, iError
      logical bGeomOpt, bAnGrad

c COMMON BLOCKS


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



c parallel_aces.com : begin

c This common block contains the MPI statistics for each MPI process. The values
c are initialized in the acescore library.

      external aces_bd_parallel_aces




      integer                nprocs, irank, icpuname

      character*(256) szcpuname

      common /parallel_aces/ nprocs, irank, icpuname,
     &                       szcpuname
      save   /parallel_aces/

c parallel_aces.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end

c ----------------------------------------------------------------------

c   o start the JOBARC subsystem
      call aces_ja_init

c   o load system size and number of displacements
      call getrec(1,'JOBARC','NUMPOINT',1,nPoint)
      call getrec(1,'JOBARC','NATOMS',1,nAtom)
      nSize = 3*nAtom

c   o set action flags
      bGeomOpt = iflags2(5).ne.0
      bAnGrad  = iflags2(138).eq.0

c   o estimate minimum memory requirements
      if (bAnGrad) then
c      - getting the derivatives
         iTmp = iintfp*nSize*3
         iTmp = iTmp + nAtom + iAnd(nAtom,1)
         iTmp = iTmp + iintfp*nSize*(1+nSize)
         iMemMin = iTmp
c      - updating the derivatives
         iTmp = iintfp*(nSize+(nSize+12)*nPoint)
         iMemMin = max(iMemMin,iTmp)
      else
         iMemMin = iintfp*nPoint*2
      end if

c   o allocate memory
      iCrSiz = iflags(36)
      iCore(1) = 0
      iCore(2) = 0
      do while ((iCore(1).eq.0.and.iCore(2).eq.0).and.
     &          (iCrSiz.gt.iMemMin))
         call aces_malloc(iCrSiz,iCore,i0)
         if (iCore(1).eq.0.and.iCore(2).eq.0) iCrSiz = iCrSiz - iMemInc
      end do
      if (iCore(1).eq.0.and.iCore(2).eq.0) then
         print *, '@SYNC_FD: unable to allocate at least ',
     &            iMemMin,' integers of memory'
         call aces_exit(1)
      end if
c      print *, '@SYNC_FD: Allocated ',iCrSiz,' integers @ iCore(',i0,')'

c   o initialize the iCore free space pointer
      lFree = i0

c   o find last displacement which was calculated
      call getrec(1,'JOBARC','FDCALCTP',nPoint,iCore(lFree))
      iOff = nPoint-1
      do while ((iCore(lFree+iOff).ge.0).and.(iOff.ge.0))
         iOff = iOff - 1
      end do
      iLast = 1+iOff
      if (iLast.eq.0.and.(bGeomOpt.or.irank.ne.0)) then
         print *, '@SYNC_FD: Assertion failed.'
         print *, '          There was no previous calculation.'
         call aces_exit(1)
      end if

c CONCEPT
c    There are 4 records that contain information needed to form numerical
c 1st or 2nd derivatives: ENGPOINT, GRDPOINT, DIPPOINT, and POLPOINT.
c    Numerical gradients require ENGPOINT. Technically, fully numerical
c Hessians also require the REFENERG record, but this point is always done by
c root process, which forms the Hessian, so we don't have to consider it here.
c Analytical gradients, which imply numerical Hessians, require GRD-, DIP-,
c and POL-POINT records.
c    sync_fd will marshall all the records into a continuous array and reduce
c it onto the master process.

      if (bAnGrad) then

c      o determine the scope of symmetry
         call getrec(-1,'JOBARC','DANGERUS',1,iTmp)
         if (iflags(79).eq.0.and.iTmp.eq.0) then
            szType='FULL'
         else
            szType='COMP'
         end if

c      o assign indices
         lGrdInt = lFree
         lFree = lFree + iintfp*nSize
         lGrdXYZ = lFree
         lFree = lFree + iintfp*nSize
         lGeom = lFree
         lFree = lFree + iintfp*nSize
         lIMAP = lFree
         lFree = lFree + nAtom + iand(nAtom,1)
         nDbls = (iCrSiz+i0-lFree)/iintfp

c      o load the original coordinates for the last displacement
         call getrec(1,'JOBARC','FDCOORDS',iintfp*nSize*iLast,
     &               iCore(lGrdInt))
         call c_memmove(iCore(lGeom),
     &                  iCore(lGrdInt+iintfp*nSize*iOff),
     &                  ifltln*nSize)

c      o load and rotate the derivatives
         call getgrd(nAtom,szType,iCore(lGeom),
     &               iCore(lGrdXYZ),iCore(lGrdInt),
     &               dDipXYZ,dPolXYZ,iCore(lIMAP),
     &               iCore(lFree),nDbls,.false.)

c      o reclaim scratch after GRDINT
         lSend = lGrdXYZ
         lRecv = lSend
         iCount = 0

c      o energy gradient
         nDbls = nSize*nPoint
         call getrec(1,'JOBARC','GRDPOINT',iintfp*nDbls,iCore(lRecv))
         call dcopy(nSize,iCore(lGrdInt),1,
     &                    iCore(lRecv+iintfp*nSize*iOff),1)
         lRecv = lRecv + iintfp*nDbls
         iCount = iCount + nDbls

c      o dipole
         nDbls = 3*nPoint
         call getrec(1,'JOBARC','DIPPOINT',iintfp*nDbls,iCore(lRecv))
         call dcopy(3,dDipXYZ,1,iCore(lRecv+iintfp*3*iOff),1)
         lRecv = lRecv + iintfp*nDbls
         iCount = iCount + nDbls

c      o polarizability
         nDbls = 9*nPoint
         call getrec(1,'JOBARC','POLPOINT',iintfp*nDbls,iCore(lRecv))
         call dcopy(9,dPolXYZ,1,iCore(lRecv+iintfp*9*iOff),1)
         lRecv = lRecv + iintfp*nDbls
         iCount = iCount + nDbls

         iRoot = 0
         call MPW_Reduce(iCore(lSend),iCore(lRecv),
     &                   iCount,MPI_DOUBLE_PRECISION,MPI_SUM,
     &                   iRoot,MPI_COMM_WORLD,iError)
         if (iError.ne.0) then
            print *, '@SYNC_FD: MPI reduction failed.'
            print *, '          return code = ',iError
            call aces_exit(1)
         end if

         if (irank.eq.0) then
            call putrec(1,'JOBARC',
     &                  'GRDPOINT',iintfp*nSize*nPoint,iCore(lRecv))
            lRecv = lRecv + iintfp*nSize*nPoint
            call putrec(1,'JOBARC',
     &                  'DIPPOINT',iintfp*3*nPoint,iCore(lRecv))
            lRecv = lRecv + iintfp*3*nPoint
            call putrec(1,'JOBARC',
     &                  'POLPOINT',iintfp*9*nPoint,iCore(lRecv))
         end if

c     else if (.not.bAnGrad) then
      else

         call getrec(1,'JOBARC','ENGPOINT',iintfp*nPoint,iCore(lFree))
         if (iLast.ne.0) then
            if (iflags(87).eq.0) then
               call getrec(1,'JOBARC','TOTENERG',iintfp,
     &                     iCore(lFree+iintfp*iOff))
            else
               call getrec(1,'JOBARC','TOTENER2',iintfp,
     &                     iCore(lFree+iintfp*iOff))
            end if
         end if
         iCount = nPoint
         iRoot  = 0
         if (bGeomOpt) then
            call MPW_Allreduce(iCore(lFree),iCore(lFree+iintfp*nPoint),
     &                         iCount,MPI_DOUBLE_PRECISION,MPI_SUM,
     &                         MPI_COMM_WORLD,iError)
         else
            call MPW_Reduce(iCore(lFree),iCore(lFree+iintfp*nPoint),
     &                      iCount,MPI_DOUBLE_PRECISION,MPI_SUM,
     &                      iRoot,MPI_COMM_WORLD,iError)
         end if
         if (iError.ne.0) then
            print *, '@SYNC_FD: MPI reduction failed.'
            print *, '          return code = ',iError
            call aces_exit(1)
         end if
         if (bGeomOpt.or.irank.eq.0) then
            call putrec(1,'JOBARC',
     &                  'ENGPOINT',iintfp*nPoint,
     &                  iCore(lFree+iintfp*nPoint))
         end if

c     end if (bAnGrad)
      end if

c   o free memory
      call c_free(iCore)

c   o finalize the JOBARC subsystem and free memory
      call aces_ja_fin

      return
c     end subroutine sync_fd
      end

