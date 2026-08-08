










      character *(*) function readline(line,unit,eof)

c This reads the next non-comment non-blank line.  If none is found,
c eof is set to true.  Otherwise, the line is returned.  A comment is
c defined as a line starting with a "#".  linenumber is the line number
c in the file.  It should be initalized to zero before this is called
c for the first time.








c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler























cYAU - ACES3 stuff . . . we hope - #include <aces.par>






c This contains the global string for identifying the current subroutine
c or function (provided the programmer set it).  cf. tools/callstack.F
c BE GOOD AND RESET CURR ON EXIT!

      character*64                callstack_curr,callstack_prev
      common /callstack_curr_com/ callstack_curr,callstack_prev
      save   /callstack_curr_com/



      integer line,unit
      logical eof

      integer strlen

      callstack_curr='READLINE'

      eof=.false.
   10 read(unit,'(a)',end=900) readline
      line=line+1
      if (strlen(readline).eq.0 .or.  readline(1:1).eq.'#') goto 10
      return

  900 eof=.true.
      return
      end

c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End:
