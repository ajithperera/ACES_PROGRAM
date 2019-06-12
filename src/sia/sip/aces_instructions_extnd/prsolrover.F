c--------------------------------------------------------------------|
c   This is a superinstruction for SIAL code                         |
c   It prints out overlap between left and right eigenvectors        |  
c   Added by Denis Bokhan, Moscow Lomonosov State University,        |
c                  MOSCOW, RUSSIAN FEDERATION                        |
c--------------------------------------------------------------------|


       subroutine prsolrover(x,nindex, type, bval,eval,
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
      integer bdim2(*), edim2(*),matr_size

      matr_size = 1 + n_si + n_tr*3
c----------------------------------------------------------------
      call do_prsolrover(x,bval(1),eval(1),type(1),
     &                     bval(2),eval(2),type(2) )
c----------------------------------------------------------------
      return
      end
C=====================================================================
C=====================================================================
      subroutine do_prsolrover(x1,a1,a2,type1,b1,b2,type2)
      implicit none
      include 'parallel_info.h'
      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
      include 'spin_orb.h'

      integer a1,a2,b1,b2,c1,c2,d1,d2,i1,i2,j1,j2,k1,k2
      integer i,j,k,l,a,b,l1,l2, msize,KK,ICK
      integer type1,type2,type3,type4,type5,type6,type7,type8
      double precision x1(a1:a2,b1:b2)


       if(me.eq.0) then
       write(6,*) 
     & '| --- Overlap between left and right eigenvectors ---|'  


           do i=a1,a2
           do j=b1,b2
            write(6,*) x1(i,j),'      ', i, j
           enddo
           enddo

       write(6,*) 
     & '| -------------------------------------------------- |'

       endif

       return
       end
