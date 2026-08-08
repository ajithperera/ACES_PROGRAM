














      subroutine MPW_Initialized(flag,ierror)
      implicit none

      logical flag
      integer ierror
      logical w_flag
      integer w_ierror

      call MPI_Initialized(w_flag,w_ierror)
      flag = w_flag
      ierror = w_ierror

      return
      end

