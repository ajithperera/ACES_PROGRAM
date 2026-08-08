c---------------------------------------------------------------|
c  This subroutine is a superinstruction for the SIAL code      |
c  Added by Denis Bokhan, Moscow Lomonosov State University,    |
c              Moscow, RUSSIAN FEDERATION                       |
c---------------------------------------------------------------|
      subroutine copy24(x,nindex, type, bval,eval, bdim, edim,
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

      if(nindex.eq.2) then
      call do_copy24(x,bval(1),eval(1),type(1),bval(2),eval(2),type(2),
     & x2,bval2(1),eval2(1),type2(1),bval2(2),eval2(2),type2(2),
     &    bval2(3),eval2(3),type2(3),bval2(4),eval2(4),type2(4) )
      endif

      if(nindex.eq.4) then
      call do_copy46(x,bval(1),eval(1),type(1),bval(2),eval(2),type(2),
     &                 bval(3),eval(3),type(3),bval(4),eval(4),type(4),
     & x2,bval2(1),eval2(1),type2(1),bval2(2),eval2(2),type2(2),
     &    bval2(3),eval2(3),type2(3),bval2(4),eval2(4),type2(4), 
     &    bval2(5),eval2(5),type2(5),bval2(6),eval2(6),type2(6) )
      endif

c----------------------------------------------------------------    
      return
      end

C----------------------------------------------------------------
C----------------------------------------------------------------

      subroutine do_copy24(x1,a1,a2,type1,b1,b2,type2,
     * x2,i1,i2,type11,j1,j2,type22,c1,c2,type33,d1,d2,type44)
      implicit none
      include 'parallel_info.h'
      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'

      integer a1,a2,b1,b2,c1,c2,d1,d2,i1,i2,j1,j2,k1,k2
      integer i,j,k,l,a,b,l1,l2,g
      integer type1,type2,type3
      integer type11,type22,type33,type44
      double precision x1(a1:a2,b1:b2),
     * x2(i1:i2,j1:j2,c1:c2,d1:d2)

c        if(me.eq.0) write(6,*) 'i1=',i1,'j1=',j1
c        if(me.eq.0) write(6,*) 'i2=',i2,'j2=',j2

        do i=i1,i2
        do j=j1,j2
         DO k=a1,a2
         DO l=b1,b2
         x2(i,j,k,l) = x1(k,l)
         ENDDO
         ENDDO
        enddo
        enddo 

      return
      end

C----------------------------------------------------------------
C----------------------------------------------------------------

      subroutine do_copy46(x1,a1,a2,type1,b1,b2,type2,c1,c2,type3,d1,d2,
     * type4,x2,i1,i2,type11,j1,j2,type22,k1,k2,type33,l1,l2,type44,
     * n1,n2,type55,m1,m2,type66)
      implicit none
      include 'parallel_info.h'
      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
C      include 'F12.h'
      include 'F12_gr.h'

      integer a1,a2,b1,b2,c1,c2,d1,d2,i1,i2,j1,j2,k1,k2
      integer l1,l2,m1,m2,n1,n2
      integer i,j,k,l,m,n,a,b,g
      integer type1,type2,type3,type4
      integer type11,type22,type33,type44,type55,type66
      double precision x1(a1:a2,b1:b2,c1:c2,d1:d2),
     * x2(i1:i2,j1:j2,k1:k2,l1:l2,n1:n2,m1:m2)

        if(me.eq.0) write(6,*) 'i1=',i1,'j1=',j1
        if(me.eq.0) write(6,*) 'i2=',i2,'j2=',j2


        do i=i1,i2
        do j=j1,j2

         DO k=a1,a2
         DO l=b1,b2
         DO n=c1,c2
         DO m=d1,d2

         x2(i,j,k,l,n,m) = x1(k,l,n,m)            

         enddo
         enddo
         enddo
         enddo

        enddo
        enddo

      return
      end
c---------------------------------------------------------------
c---------------------------------------------------------------
