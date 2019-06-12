c---------------------------------------------------------------|
c  This subroutine is a superinstruction for the SIAL code      |
c  It returns array A(a,mu,nu) = CVEC(mu,a)*CVEC(nu,a)          |
c  Added by Denis Bokhan, Moscow Lomonosov State University,    |
c              Moscow, RUSSIAN FEDERATION                       |
c---------------------------------------------------------------|
      subroutine mksqcvec(x,nindex, type, bval,eval, bdim, edim,
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

      call cvec2(x, bval(1),eval(1),type(1),
     *                   bval(2),eval(2),type(2),
     *                   bval(3),eval(3),type(3),
     *                   bval(4),eval(4),type(4),
     *            x2, bval2(1),eval2(1),type2(1),
     *                bval2(2),eval2(2),type2(2),
     *                bval2(3),eval2(3),type2(3)  )

      return
      end     

c---------------------------------------------------------------
c---------------------------------------------------------------

      subroutine cvec2(x1,a1,a2,type1,b1,b2,type2,c1,c2,type3,
     *  d1,d2,type4,x2,i1,i2,type5,j1,j2,type6,k1,k2,type7)
      implicit none
      include 'parallel_info.h'        

      integer a1,a2,b1,b2,c1,c2,d1,d2,i1,i2,j1,j2,k1,k2  
      integer i,j,k,l
      integer type1,type2,type3,type4,type5,type6,type7
      double precision x1(a1:a2,b1:b2,c1:c2,d1:d2),
     * x2(i1:i2,j1:j2,k1:k2)

      if(me.eq.0) write(6,*) a1,a2,'==========>',d1,d2  

      do i=a1,a2
      do j=b1,b2
      do k=c1,c2
      x2(i,j,k) = x1(i,j,k,i)
      enddo
      enddo
      enddo      

      return
      end  
c---------------------------------------------------------------
