










        Program main

      implicit none
      integer iuhf,ihfdftgrad
     


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
         call aces_init(icore,i0,icrsiz,iuhf,.true.)
         
         call getrec(20,'JOBARC','HFDFTGRA',1,ihfdftgrad)
           if (ihfdftgrad.eq.1) then
             call vhfksdint(icore(i0),icrsiz,iuhf)
             call vksdint
           else
             call vksdint
           end if
        call aces_fin
        end             
