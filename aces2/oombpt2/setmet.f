











      subroutine setmet(singles,include_fai,listv,cycle)



































































































































































































      implicit none
C Common blocks
      integer iflags(100)
      common /flags/iflags
      integer iflags2(500)
      common /flags2/iflags2
C Output variables
      integer listv,cycle
      logical singles, include_fai

      listv = 93
      if (iflags(2) .eq. 45) then

C        singles = .true.
C        include_fai = .false.
C
        singles = .false.
        include_fai = .true.

      else if (iflags(2) .eq. 46) then

        singles = .true.
        include_fai = .false.
C
C        singles = .false.
C        include_fai = .true.
      else
        write(6,*) 'Only OO-MP2 and OO-MBPT(2) should call xoombpt2'
        call errex
      endif
C Is this the first cycle?
      call getrec(-1,'JOBARC','OOMBPTIT',1,cycle)
      if (cycle .eq. 0) then
        if (iflags2(166) .gt. 0) then
          call putrec(20,'JOBARC','OOMBPTCY',1,
     &                                      iflags2(166))
        else
          call putrec(20,'JOBARC','OOMBPTCY',1,
     &                                      iflags(16))
        endif
      endif

      return
      end
