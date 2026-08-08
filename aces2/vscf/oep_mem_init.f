










      subroutine oep_mem_init(maxmem, irppmem, imemget)
c 
      implicit double precision (a-h, o-z)
c


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



c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
      maxmem  = iflags (36) - irppmem - iflags(44)*iflags(37)
      maxmem  = maxmem      + mod (maxmem, 2)
      imemget = maxmem      + irppmem + iflags(44)*iflags(37)
      imemget = imemget     + mod (imemget, 2)
     
      Return
      End

