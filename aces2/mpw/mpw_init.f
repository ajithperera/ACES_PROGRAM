














      subroutine MPW_Init(ierror)
      implicit none

      integer ierror
      integer w_ierror

      call MPI_Init(w_ierror)
      ierror = w_ierror

      return
      end

