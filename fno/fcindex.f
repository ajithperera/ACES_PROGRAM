










      function fcindex(str,char)







c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler























cYAU - ACES3 stuff . . . we hope - #include <aces.par>




      integer fcindex
      character *(*) str
      character*1 char
      integer i,strlen,slen
c      callstack_curr='FCINDEX'
      i=1
      slen=strlen(str)
 10   if (i.gt.slen) then
        fcindex=0
        return
      end if
      if (str(i:i).eq.char) then
        fcindex=i
        return
      end if
      i=i+1
      goto 10
      end

c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End:
