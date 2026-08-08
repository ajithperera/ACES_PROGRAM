











      program main
      implicit none



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
c method.com : begin
      integer iuhf
      logical scf, nonhf
      common /method/ iuhf, scf, nonhf
c method.com : end

      call crapsi(icore,iuhf,0)
      call nmr(icore(i0),icrsiz,iuhf)
      call aces_fin

      stop
      end

