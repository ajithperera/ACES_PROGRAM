











c This routine drives the AO-based particle-particle ladder contraction(s)
c using a single-pass algorithm.




c#define _DEBUG_AOLAD2

      subroutine aolad2(iCore,iCoreDim,iUHF,bTau,irp_x,
     &                  list_ao,list_ao_inc)
      implicit none

      integer iCore(*), iCoreDim, iUHF, irp_x, list_ao, list_ao_inc
      logical bTau

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

      integer iUnit, iLnBuf
      integer nSizeT, iSizT2Mix
      integer iSpin, irrep, nMO, nAO, i
      integer irp_ij, max_ij
      integer irp_ab, max_ab
      integer iOff, ndx, iTmp, iTmp1, iTmp2, iTmp3, iTmp4
      integer iRoots(8), iBchOff(8), iBchLen(8)
      integer displs(256), recvcounts(256), iErr

      integer pFree,
     &        pdBf, piBf,
     &        pIK0, pIL0, pJK0, pJL0, pKI0, pKJ0,
     &        pVal, pSym, pTyp, pWhr, pAOS, pMap,
     &        p_T2, p_dZ
      character*80 szFileName
      integer       iFileName
      logical bInCore

      INTEGER        IAOPOP(8),IOFFAO(8),IOFFV(8,2),IOFFO(8,2),
     &               IRPDPDAO(8),IRPDPDAOS(8),ISTART(8,8),ISTARTMO(8,3)
      COMMON /AOSYM/ IAOPOP,   IOFFAO,   IOFFV,     IOFFO,
     &               IRPDPDAO,   IRPDPDAOS,   ISTART,     ISTARTMO


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
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end
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

c ----------------------------------------------------------------------

c   o assign a unit number to the AO integral file(s)
      iUnit = 10

c   o retrieve the numbers of MOs and AOs
      nMO = nocco(1) + nvrto(1)
      call getrec(1,'JOBARC','NBASTOT',1,nAO)

c   o fill in the AOSYM common block
      call getaoinf(iUHF,irp_x)

c   o set the buffer length
      iLnBuf = 600

c   o calculate the total number of back-transformed T2 amplitudes
      if (iUHF.eq.0) then
         nSizeT = 0
         do irp_ij = 1, nirrep
            irp_ab = dirprd(irp_ij,irp_x)
            nSizeT =   nSizeT
     &               + IRPDPDAO(irp_ab)*irpdpd(irp_ij,isytyp(2,46))
         end do
      else
         nSizeT = 0
         do irp_ij = 1, nirrep
            irp_ab = dirprd(irp_ij,irp_x)
            iTmp =   irpdpd(irp_ij,isytyp(2,44))
     &             + irpdpd(irp_ij,isytyp(2,45))
     &             + irpdpd(irp_ij,isytyp(2,46))
            nSizeT = nSizeT + IRPDPDAO(irp_ab)*iTmp
         end do
      end if

c ----------------------------------------------------------------------

c   o create the memory map
      pFree = 1

c   o allocate the integral and index buffers
      pdBf  = pFree
      piBf  = pdBf + max(NT(1),NT(2),iLnBuf,nAO*nMO)*iintfp
      pFree = piBf + max(NT(1),NT(2),iLnBuf,nAO*nMO*iintfp)

c   o allocate temp arrays for AO pair-indices
      pIK0  = pFree
      pIL0  = pIK0 + iLnBuf*8
      pJK0  = pIL0 + iLnBuf*8
      pJL0  = pJK0 + iLnBuf*8
      pKI0  = pJL0 + iLnBuf*8
      pKJ0  = pKI0 + iLnBuf*8
      pFree = pKJ0 + iLnBuf*8

c   o ???
      pVal  = pFree
      pSym  = pVal + iLnBuf*8*iintfp
      pTyp  = pSym + iLnBuf*8
      pWhr  = pTyp + iLnBuf*8
      pAOS  = pWhr + iLnBuf*8
c      pMap  = pAOS + max(iLnBuf*8,nAO)
      pMap  = pAOS + nAO
      pFree = pMap + nAO*nAO

c   o shift pFree to a double boundary (assuming iCore(1) is aligned)
      pFree = pFree + 1 - iand(pFree,1)

c   o restrict the implementation to incore (holding all T2 and Z arrays)
      iTmp = (iCoreDim+1-pFree)/iintfp
      if (iTmp.lt.nSizeT*2) then
         print *, '@AOLAD2: insufficient memory'
         print *, '         need ',nSizeT*2,' dbls, have ',iTmp
         print *, '         retry using AO_LADDERS=MULTIPASS'
         call aces_exit(1)
      end if
      bInCore = .true.
      p_T2  = pFree
      p_dZ  = p_T2 + nSizeT*iintfp
      pFree = p_dZ + nSizeT*iintfp

c ----------------------------------------------------------------------

c   o initialize iAOSym
      ndx = pAOS
      do irrep = 1, nirrep
         do i = 1, IAOPOP(irrep)
            iCore(ndx) = irrep
            ndx = ndx + 1
         end do
      end do

c   o initialize iMap
      call aosymvec(iCore(pMap),nAO)

c   o load the T2(XX,I<J) amplitudes and transpose them to T2(I<J,XX)
      iOff = 0
      ndx  = p_T2
      do irp_ij = 1, nirrep
         irp_ab = dirprd(irp_ij,irp_x)
         max_ij = irpdpd(irp_ij,isytyp(2,46))
         max_ab = IRPDPDAO(irp_ab)
         call getlst(iCore(p_dZ),1,max_ij,1,irp_ij,list_ao+3)
         call transp(iCore(p_dZ),iCore(ndx),max_ij,max_ab)
         ndx = ndx + iintfp*max_ab*max_ij
         iRoots(irp_ab) = iOff
         iOff = iOff + mod(max_ab,nprocs)
         iOff = mod(iOff,nprocs)
         call paces_batch_stat(irank,nprocs,iRoots(irp_ab),max_ab,
     &                         iBchOff(irp_ab),iBchLen(irp_ab))
      end do
      if (iUHF.ne.0) then
         do iSpin = 2, 1, -1
         do irp_ij = 1, nirrep
            irp_ab = dirprd(irp_ij,irp_x)
            max_ij = irpdpd(irp_ij,isytyp(2,43+iSpin))
            max_ab = IRPDPDAO(irp_ab)
            call getlst(iCore(p_dZ),1,max_ij,1,irp_ij,list_ao+iSpin)
            call transp(iCore(p_dZ),iCore(ndx),max_ij,max_ab)
            ndx = ndx + iintfp*max_ab*max_ij
         end do
         end do
      end if

c   o initialize the T2 increments
      call zero(iCore(p_dZ),nSizeT)

c   o open the integral files and process each

      call gfname('IIII',szFileName,iFileName)
      open(unit=iUnit,file=szFileName(1:iFileName),
     &     form='UNFORMATTED',status='OLD')
      call rdaoijkl2(iCore(p_T2),iCore(p_dZ),iCore(pdBf),iCore(piBf),
     &               iCore(pIK0),iCore(pIL0),iCore(pJK0),iCore(pJL0),
     &               iCore(pKI0),iCore(pKJ0),iCore(pVal),iCore(pSym),
     &               iCore(pTyp),iCore(pWhr),iCore(pAOS),iCore(pMap),
     &               iLnBuf,iUnit,iUHF,nAO,iSizT2Mix,irp_x,
     &               iBchOff,iBchLen)
      close(unit=iUnit,status='KEEP')

      if (nirrep.gt.1) then

      call gfname('IJIJ',szFileName,iFileName)
      open(unit=iUnit,file=szFileName(1:iFileName),
     &     form='UNFORMATTED',status='OLD')
      call rdaoijkl2(iCore(p_T2),iCore(p_dZ),iCore(pdBf),iCore(piBf),
     &               iCore(pIK0),iCore(pIL0),iCore(pJK0),iCore(pJL0),
     &               iCore(pKI0),iCore(pKJ0),iCore(pVal),iCore(pSym),
     &               iCore(pTyp),iCore(pWhr),iCore(pAOS),iCore(pMap),
     &               iLnBuf,iUnit,iUHF,nAO,iSizT2Mix,irp_x,
     &               iBchOff,iBchLen)
      close(unit=iUnit,status='KEEP')

      call gfname('IIJJ',szFileName,iFileName)
      open(unit=iUnit,file=szFileName(1:iFileName),
     &     form='UNFORMATTED',status='OLD')
      call rdaoijkl2(iCore(p_T2),iCore(p_dZ),iCore(pdBf),iCore(piBf),
     &               iCore(pIK0),iCore(pIL0),iCore(pJK0),iCore(pJL0),
     &               iCore(pKI0),iCore(pKJ0),iCore(pVal),iCore(pSym),
     &               iCore(pTyp),iCore(pWhr),iCore(pAOS),iCore(pMap),
     &               iLnBuf,iUnit,iUHF,nAO,iSizT2Mix,irp_x,
     &               iBchOff,iBchLen)
      close(unit=iUnit,status='KEEP')

      if (nirrep.gt.2) then

      call gfname('IJKL',szFileName,iFileName)
      open(unit=iUnit,file=szFileName(1:iFileName),
     &     form='UNFORMATTED',status='OLD')
      call rdaoijkl2(iCore(p_T2),iCore(p_dZ),iCore(pdBf),iCore(piBf),
     &               iCore(pIK0),iCore(pIL0),iCore(pJK0),iCore(pJL0),
     &               iCore(pKI0),iCore(pKJ0),iCore(pVal),iCore(pSym),
     &               iCore(pTyp),iCore(pWhr),iCore(pAOS),iCore(pMap),
     &               iLnBuf,iUnit,iUHF,nAO,iSizT2Mix,irp_x,
     &               iBchOff,iBchLen)
      close(unit=iUnit,status='KEEP')

c     end if (nirrep.gt.2)
      end if

c     end if (nirrep.gt.1)
      end if

      if (bInCore) then
         if (iUHF.eq.0) then
            ndx = p_dZ
            do irp_ij = 1, nirrep
               irp_ab = dirprd(irp_ij,irp_x)
               max_ij = irpdpd(irp_ij,isytyp(2,46))
               max_ab = IRPDPDAO(irp_ab)
               iTmp   = iintfp * max(max_ab,max_ij)
               iTmp1  = pdBf
               iTmp2  = iTmp1 + iTmp
               iTmp3  = iTmp2 + iTmp
               iTmp4  = iTmp3 + iTmp
               if (iTmp4.gt.p_T2) then
                  print *, '@AOLAD2: insufficient memory'
                  call aces_exit(1)
               end if
               call symrhf3(irp_ab,irp_ij,pop(1,1),IAOPOP,max_ij,
     &                      iCore(ndx),
     &                      iCore(iTmp1),iCore(iTmp2),iCore(iTmp3))
               call transp(iCore(ndx),iCore(p_T2),max_ab,max_ij)
               call putlst(iCore(p_T2),1,max_ij,1,irp_ij,list_ao_inc+3)
               ndx = ndx + iintfp*max_ab*max_ij
            end do
         else
            ndx = p_dZ
            do iSpin = 3, 1, -1
            do irp_ij = 1, nirrep
               irp_ab = dirprd(irp_ij,irp_x)
               max_ij = irpdpd(irp_ij,isytyp(2,43+iSpin))
               max_ab = IRPDPDAO(irp_ab)
               call transp(iCore(ndx),iCore(p_T2),max_ab,max_ij)
               call putlst(iCore(p_T2),1,max_ij,1,irp_ij,
     &                     list_ao_inc+iSpin)
               ndx = ndx + iintfp*max_ab*max_ij
            end do
            end do
         end if
c     end if (bInCore)
      end if

      return
c     end subroutine aolad2
      end

