














      subroutine MPW_Finalize(ierror)
      implicit none

      integer ierror
      integer w_ierror

      call MPI_Finalize(w_ierror)
      ierror = w_ierror

      return
      end

