c---------------------------------------------------------------|
c  This subroutine is a superinstruction for the SIAL code      |
c  It performs calculation of a sum of type:                    |
c                    D = sum_{i} <i|d|i>                        |
c  Added by Denis Bokhan, Moscow Lomonosov State University,    |
c              Moscow, RUSSIAN FEDERATION                       |
c---------------------------------------------------------------|
      subroutine dipole_scalar(x,nindex, type, bval,eval, bdim, edim,
     *           x2, nindex2, type2, bval2, eval2, bdim2, edim2)


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

      call do_dipole_scalar(x, bval(1),eval(1),type(1),
     *                         bval(2),eval(2),type(2),
     *                      x2 )
c----------------------------------------------------------------
      return
      end

C----------------------------------------------------------------
C----------------------------------------------------------------

      subroutine do_dipole_scalar(x1,a1,a2,type1,b1,b2,type2,
     * x2)
      implicit none
      include 'parallel_info.h'
      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'

      integer a1,a2,b1,b2,c1,c2,d1,d2,i1,i2,j1,j2,k1,k2
      integer i,j,k,l,a,b,l1,l2
      integer type1,type2,type3
      double precision x1(a1:a2,b1:b2),
     * x2

c------------Calculate sum itselfs-------------------------------
      write(6,*) type1,type2
c      if(me.eq.0) write(6,*) j,'pppppp', a1,a2, b1,b2

      x2 = 0.d0 
      do i=a1,a2
      do k=b1,b2
c      if(i.eq.k.and.me.eq.0) write(6,*) i,k,'zzzzz' 
      if(i.eq.k) x2 = x2 + x1(i,k)  
      enddo
      enddo 
c      if(me.eq.0) write(6,*) 'x2=',x2

      return
      end


