










      subroutine upcase(str)







c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler























cYAU - ACES3 stuff . . . we hope - #include <aces.par>




      character*(*) str
      integer l,i
c      callstack_curr='UPCASE'
      l=len(str)
      do 10 i=1,l
        if (str(i:i).ge.'a' .and. str(i:i).le.'z')
     &      str(i:i)=char(ichar(str(i:i))+ichar('A')-ichar('a'))
   10 continue
      return
      end

c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End:
