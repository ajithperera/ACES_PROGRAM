










      function tfkef(roa,rob,grdaa,grdab,grdbb)
      
c Determine the value of the Thomas Fermi kinetic energy functional.
      
      implicit none







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



      double precision
     &    roa,rob,grdaa,grdab,grdbb

      double precision
     &    tfkef,ro,tfcnst,pi

      callstack_curr='TFKEF'
      pi = acos(-1.d0)
      ro=roa+rob
      tfcnst=(3.d0/10.d+00)*(3.d0*pi*pi)**(2.d0/3.d0)
      tfkef=tfcnst*(ro**(5.d0/3.d0))
     
      return
      end
