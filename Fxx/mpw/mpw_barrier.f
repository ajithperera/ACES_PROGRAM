














      subroutine MPW_Barrier(comm,ierror)
      implicit none

      integer comm, ierror
      integer w_comm, w_ierror

      w_comm = comm
      call MPI_Barrier(w_comm,w_ierror)
      ierror = w_ierror

      return
      end

