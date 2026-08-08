










      subroutine stretch(density,nso,nao)

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





      integer
     &   nso, nao, i, offset
      double precision
     &   density(nao*nao*2)

      offset=nao**2-nso**2
      do i=(nao**2+nso**2),(nao**2+1),(-1)
         density(i)=density(i-offset)
      end do

      return
      end
