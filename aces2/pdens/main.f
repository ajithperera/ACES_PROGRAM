










      program pdens_main
      implicit none



































































































































































































c


c icore.com : begin

c icore(1) is an anchor in memory that allows subroutines to address memory
c allocated with malloc. This system will fail if the main memory is segmented
c or parallel processes are not careful in how they allocate memory.

      integer icore(1)
      common / / icore

c icore.com : end





c istart.com : begin
      integer         i0, icrsiz
      common /istart/ i0, icrsiz
      save   /istart/
c istart.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
      integer iuhf
c ----------------------------------------------------------------------
      call aces_init_rte
      call aces_com_parallel_aces
      call aces_ja_init
      call getrec(1,'JOBARC','IFLAGS', 100,iflags)
      call getrec(1,'JOBARC','IFLAGS2',500,iflags2)
c ----------------------------------------------------------------------

      icrsiz = iflags(36)
      icore(1) = 0
      do while ((icore(1).eq.0).and.(icrsiz.gt.1000000))
         call aces_malloc(icrsiz,icore,i0)
         if (icore(1).eq.0) icrsiz = icrsiz - 1000000
      end do
      if (icore(1).eq.0) then
         print *, '@MAIN: unable to allocate at least ',
     &            1000000,' integers of memory'
         call aces_exit(1)
      end if

c ----------------------------------------------------------------------

      call aces_io_init(icore,i0,icrsiz,.true.)

      if (iflags(11) .eq. 0) then
         iuhf = 0
      else
        iuhf = 1
      endif
c
         print*, "Entering PERT_DENS_MAIN" 
         print*, "The iuhf:", iuhf
      Write(6,*)
      Write(6,"(T2,a)") "Perturbed densities are computed"

         call aces_init_chemsys
         call pert_dens_main(icore(i0), icrsiz, iuhf)
c ----------------------------------------------------------------------

 9997 continue
      call aces_io_fin
 9998 continue
c      call c_free(icore)
 9999 continue
      call aces_ja_fin
C
      Write(6,*)
      Write(6,"(T2,a)") "Perturbed densities were successfully computed"
      call c_exit(0)
      end

