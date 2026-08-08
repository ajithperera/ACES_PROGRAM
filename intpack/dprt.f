











      subroutine dprt(nrow,ncol,mat,title)







c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler























cYAU - ACES3 stuff . . . we hope - #include <aces.par>




      integer nrow,ncol
      double precision mat(nrow,ncol)
      character *(*) title
      integer i,j
      write (6,*) '*** ',title
      do i=1,nrow
        write(6,100) (mat(i,j),j=1,ncol)
      end do
      write(6,*)
  100 format (15(f9.6,1x))
      return
      end
