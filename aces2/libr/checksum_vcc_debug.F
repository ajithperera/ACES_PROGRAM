      subroutine checksum_vcc_debug(sTag,a,n)
      implicit none

      character*(*) sTag
      integer n
      double precision a(n)

      double precision sm, nm
      integer i
      sm=0.d0
      nm=0.d0
      do i=1,n
         sm=sm+a(i)
         nm=nm+a(i)*a(i)
      end do
      write(*,'(2a,1x,2e30.16)') '@CHECKSUM: ',sTag,sm,sqrt(nm)

      return
      end
