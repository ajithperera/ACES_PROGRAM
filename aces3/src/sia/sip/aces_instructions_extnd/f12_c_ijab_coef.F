c---------------------------------------------------------------|
c  This subroutine is a superinstruction for the SIAL code      |
c  It returns arrays C(i,j,a,b) = 3/8 and C(i,a,j,b) = 1/8      |
c  which are necessary for EA-EOM-CCSD(F12) calculation         |
c  Added by Denis Bokhan, Moscow Lomonosov State University,    |
c              Moscow, RUSSIAN FEDERATION                       |
c---------------------------------------------------------------|
      subroutine f12_c_ijab_coef(x,nindex, type, bval,eval, bdim, edim,
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

      call do_f12_c_ijab_coef(x, bval(1),eval(1),type(1),
     *                   bval(2),eval(2),type(2),
     *                   bval(3),eval(3),type(3),
     *                   bval(4),eval(4),type(4),
     *            x2, bval2(1),eval2(1),type2(1),
     *                bval2(2),eval2(2),type2(2),
     *                bval2(3),eval2(3),type2(3),  
     *                bval2(4),eval2(4),type2(4)  )

      return
      end     

c---------------------------------------------------------------
c---------------------------------------------------------------

      subroutine do_f12_c_ijab_coef(x1,a1,a2,type1,b1,b2,type2,c1,c2,
     * type3,d1,d2,type4,x2,i1,i2,type5,j1,j2,type6,k1,k2,type7,
     * l1,l2,type8)
      implicit none
      include 'parallel_info.h'        
      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
C      include 'F12.h'

      integer a1,a2,b1,b2,c1,c2,d1,d2,i1,i2,j1,j2,k1,k2  
      integer i,j,k,l,a,b,l1,l2
      integer type1,type2,type3,type4,type5,type6,type7,type8
      double precision x1(a1:a2,b1:b2,c1:c2,d1:d2),
     * x2(i1:i2,j1:j2,k1:k2,l1:l2)

c------------C_ijab matrix------------------------------------
      do i=a1,a2
      do j=b1,b2
      do a=c1,c2
      do b=d1,d2 
      x1(i,j,a,b) = 0.d0
      if(i.eq.j.and.a.eq.b) x1(i,j,a,b)=-(3.d0/8.d0)/F12GAMMA
      enddo
      enddo
      enddo
      enddo 
c------------C_iajb matrix------------------------------------
      do i=i1,i2
      do a=j1,j2
      do j=k1,k2
      do b=l1,l2
      x2(i,a,j,b) = 0.d0
      if(i.eq.j.and.a.eq.b) x2(i,a,j,b) = -(1.d0/8.d0)/F12GAMMA
      enddo
      enddo
      enddo
      enddo
c-------------------------------------------------------------


      return
      end  
c---------------------------------------------------------------
