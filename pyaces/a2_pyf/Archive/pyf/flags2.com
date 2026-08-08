c flags2.com : begin
      integer*8       iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
