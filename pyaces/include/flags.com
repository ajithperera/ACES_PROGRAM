c flags.com : begin
      integer*8      iflags(100)
      common /flags/ iflags
c flags.com : end
