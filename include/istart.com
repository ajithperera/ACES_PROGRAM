c istart.com : begin
      integer*8       i0, icrsiz
      common /istart/ i0, icrsiz
      save   /istart/
c istart.com : end
