











c This routine destroys column I of an array A(N,M) by shifting back the
c column range I+1 through M.

      subroutine reduce(a,n,i,m)
      implicit none
      integer n, i, m
      double precision a(n,m)


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



      if ((n.lt.1).or.(m.lt.i).or.(i.lt.1)) then
         print *, '@REDUCE: Assertion failed.'
         print *, '   rm col index = ',i
         print *, '   rows = ',n
         print *, '   cols = ',m
         call aces_exit(1)
      end if
      if (i.eq.m) return
      call c_memmove(a(1,i),a(1,i+1),n*(m-i)*ifltln)
      return
      end

