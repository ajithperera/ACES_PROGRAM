










      program main   
      implicit none
      Integer iuhf, i0, ICoreDim



c icore.com : begin

c icore(1) is an anchor in memory that allows subroutines to address memory
c allocated with malloc. This system will fail if the main memory is segmented
c or parallel processes are not careful in how they allocate memory.

      integer icore(1)
      common / / icore

c icore.com : end






      call aces_init(iCore,i0,iCoreDim,iUHF,.true .)
      call lcc(iCore(i0),iCoreDim,iUHF)
      call aces_fin

      end

