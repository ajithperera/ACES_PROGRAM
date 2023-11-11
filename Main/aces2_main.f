














































































































































































































      Subroutine aces2_main()
      implicit none

      integer iuhf 
c COMMON BLOCKS
c icore.com : begin

c icore(1) is an anchor in memory that allows subroutines to address memory
c allocated with malloc. This system will fail if the main memory is segmented
c or parallel processes are not careful in how they allocate memory.

      integer icore(1)
      common / / icore

c icore.com : end
c istart.com : begin
      integer*8       i0, icrsiz
      common /istart/ i0, icrsiz
      save   /istart/
c istart.com : end
c flags.com : begin
      integer*8      iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer*8       iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
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

c ----------------------------------------------------------------------
         Call Joda_()
         Call Aces_ja_fin
         Call aces_init(icore,i0,icrsiz, iuhf, .true.)
CSSS         Call Test_py(icore(i0),Icrsiz,Scr,Ndim,Iuhf)
c ----------------------------------------------------------------------
      Return
      End 

