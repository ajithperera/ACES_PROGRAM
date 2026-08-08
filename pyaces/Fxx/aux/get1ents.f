










C*************************************************************
      subroutine Get1EInts(Label,Oneh,ldim)
C     
C     Get one electron integrals
C         
C*************************************************************
      implicit double precision (a-h,o-z)
c
      logical FileExist
      character*8 Name 
c
      Dimension buf(600),ibuf(600)
      Dimension Oneh(ldim)
c
      ilnbuf=600
      FileExist=.false.
      inquire(file='IIII',exist=FileExist)
      if (FileExist) then
         open(unit=10,file='IIII',form='UNFORMATTED',
     &        access='SEQUENTIAL')
         rewind 10
         
         call locate(10,Label)
         call zero(oneh,ldim)
         nut = 1
         do while (nut.gt.0)
            read(10) buf, ibuf, nut
            do int = 1, nut
               oneh(ibuf(int)) = buf(int)
            end do
         end do
      else
         write(*,*) 'Vmol IIII integral file does not exist'
      end if

      return
      end
