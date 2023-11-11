










      Subroutine Runints(Iuhf)

      Implicit None
      Integer Iuhf

c icore.com : begin

c icore(1) is an anchor in memory that allows subroutines to address memory
c allocated with malloc. This system will fail if the main memory is segmented
c or parallel processes are not careful in how they allocate memory.

      integer icore(1)
      common / / icore

c icore.com : end
c istart.com : begin
      integer*8       i0, icrsiz
      common /istart/ i0, icrsiz
      save   /istart/
c istart.com : end
    
      Call Ints_(Icore(i0),Icrsiz,Iuhf)
      Call V2ja(icore(i0),Icrsiz,Iuhf)

      Return
      End

