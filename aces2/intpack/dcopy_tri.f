










      subroutine dcopy_tri(lenght,old,new)







c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler























cYAU - ACES3 stuff . . . we hope - #include <aces.par>




       
      integer lenght

      double precision old (lenght,lenght,lenght),
     &       new (lenght,lenght,lenght)

      integer i,j,k
      
      do i=1,lenght
         do j=1,lenght
            do k=1,lenght
               new(i,j,k)=old(i,j,k)
            end do
         end do
      end do
      return
      end
