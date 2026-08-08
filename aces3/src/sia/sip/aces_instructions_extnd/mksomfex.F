c--------------------------------------------------------------------|
c   This is a superinstruction for SIAL code                         |
c   It construct exchange-like matrix elements for the SOMF          |
c                     Fock operator                                  |
c   Added by Denis Bokhan, Moscow Lomonosov State University,        |
c                  MOSCOW, RUSSIAN FEDERATION                        |
c--------------------------------------------------------------------|

      subroutine mksomfex(x,nindex, type, bval,eval,
     & bdim, edim, x2, nindex2, type2, bval2, eval2, bdim2, edim2)

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'spin_orb.h'

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)
      double precision x2(*)
      integer nindex2, type2(*), bval2(*), eval2(*)
      integer bdim2(*), edim2(*)

c        write(6,*) bval(1),eval(1),bval(2),eval(2),'|||||||'
c        write(6,*) bval(3),eval(3),bval(4),eval(4),'|||||||'
c        write(6,*) bval2(1),eval2(1),bval2(2),eval2(2),'++'


c-----------------------------------------------------------
      call do_mksomfex(x,bval(1),eval(1),type(1),
     &                   bval(2),eval(2),type(2),
     &                   bval(3),eval(3),type(3),
     &                   bval(4),eval(4),type(4), 
     &                x2,bval2(1),eval2(1),type2(1),
     &                   bval2(2),eval2(2),type2(2) )

c-----------------------------------------------------------
      return
      end

c-------------------------------------------------------------
c-------------------------------------------------------------  

      subroutine do_mksomfex(x1,a1,a2,type1,b1,b2,type2,
     & c1,c2,type6,d1,d2,type7,x2,i1,i2,type3,j1,j2,type4)
      implicit none
      include 'parallel_info.h'
      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
      include 'spin_orb.h'

      integer a1,a2,b1,b2,c1,c2,d1,d2,i1,i2,j1,j2,k1,k2
      integer i,j,k,l,a,b,l1,l2, msize,KK,ICK
      integer type1,type2,type3,type4,type5,type6,type7,type8
      double precision x1(a1:a2,b1:b2,c1:c2,d1:d2),x2(i1:i2,j1:j2)

       do k=c1,c2
       do l=d1,d2
       x2(k,l) = 0.d0
       enddo
       enddo

c--------------Do summation over occupied index--------------------
      do i=a1,a2
      do j=b1,b2
       do k=c1,c2
       do l=d1,d2
        IF(i.eq.j) then
c        write(6,*) x1(i,j,k,l), '--->',i,j,k,l 
        x2(k,l) = x2(k,l) + x1(i,j,k,l)
        endif
       enddo
       enddo
      enddo  
      enddo
c------------------------------------------------------------------
       do k=c1,c2
       do l=d1,d2
c       write(6,*) x2(k,l), '--->', k,l 
       enddo
       enddo


      return
      end
