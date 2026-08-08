c--------------------------------------------------------------------|
c   This is a superinstruction for SIAL code                         |
c   It prints out EOM-CC roots to the summary.out file               |
c   Added by Denis Bokhan, Moscow Lomonosov State University,        |
c                  MOSCOW, RUSSIAN FEDERATION                        |
c--------------------------------------------------------------------|


      subroutine print_eomso_roots(x,nindex, type, bval,eval, 
     & bdim, edim, x2, nindex2, type2, bval2, eval2, bdim2, edim2)

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)
      double precision x2(*)
      integer nindex2, type2(*), bval2(*), eval2(*)
      integer bdim2(*), edim2(*)
c----------------------------------------------------------------
      call do_print_eomso_roots(x, bval(1),eval(1),type(1),
     &                             bval(2),eval(2),type(2),
     &                         x2, bval2(1),eval2(1),type2(1),
     &                             bval2(2),eval2(2),type2(2) )
c----------------------------------------------------------------
      return
      end 
C=====================================================================
C=====================================================================
      subroutine do_print_eomso_roots(x1,a1,a2,type1,b1,b2,type2,
     &                                x2,i1,i2,type3,j1,j2,type4)
      implicit none
      include 'parallel_info.h'
      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'

      integer a1,a2,b1,b2,c1,c2,d1,d2,i1,i2,j1,j2,k1,k2
      integer i,j,k,l,a,b,l1,l2
      integer type1,type2,type3,type4,type5,type6,type7,type8
      double precision x1(a1:a2,b1:b2),
     & x2(i1:i2,j1:j2)


        if(me.eq.0) then
        open(66,file='summary.out')
        write(66,*) '--------------------------------------------------'

        do i = a1,a2
        do j= b1, b2

        if(x2(i,j).eq.1.d0) then 
        write(66,*) i, x1(i,j)*27.211385,'------------','(T)'
        endif 

        if(x2(i,j).eq.2.d0) then
        write(66,*) i, x1(i,j)*27.211385,'------------','(S)'
        endif

        enddo 
        enddo
        write(66,*) '--------------------------------------------------'
        close(66)
        endif
      return
      end  
