c---------------------------------------------------------------|
c  This subroutine is a superinstruction for the SIAL code      |
c  It performs multiplication of array by x value               |
c  Added by Denis Bokhan, Moscow Lomonosov State University,    |
c              Moscow, RUSSIAN FEDERATION                       |
c---------------------------------------------------------------|
      subroutine efgninty(x,nindex, type, bval,eval, bdim, edim,
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

      call do_efgninty(x, bval(1),eval(1),type(1),
     *                    bval(2),eval(2),type(2),x2 )
c----------------------------------------------------------------
      return
      end

C----------------------------------------------------------------
C----------------------------------------------------------------

      subroutine do_efgninty(x1,a1,a2,type1,b1,b2,type2,
     * x2)
      implicit none
      include 'parallel_info.h'
      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
      include 'F12_gr.h'
      include 'efgn.h'

      integer a1,a2,b1,b2,c1,c2,d1,d2,i1,i2,j1,j2,k1,k2
      integer i,j,k,l,a,b,l1,l2,g,kk1,ll1
      integer type1,type2,type3,type4,iatom
      double precision x1(a1:a2,b1:b2),qq,rsq,
     * x2,gx,gy,gz,XCOORD,YCOORD,ZCOORD,RX,RY,RZ,sqroot,
     * den


C-------------Define index of requested nucleus----------------------------------
      iatom = atindex 
C-------------Get atomic coordinates---------------------------------------------
            XCOORD = acenter(iatom,1)
            YCOORD = acenter(iatom,2)
            ZCOORD = acenter(iatom,3)
C--------------------------------------------------------------------------------
      g = x2
      gx  = PTX(g)
      gy  = PTY(g)
      gz  = PTZ(g)  
C--------------------------------------------------------------------------------
      RX = gx-XCOORD
      RY = gy-YCOORD
      RZ = gz-ZCOORD

      rsq = RX*RX+RY*RY+RZ*RZ
      sqroot = DSQRT(rsq)  
      den = sqroot**5.d0

      qq = (3.d0*RY*RY - rsq)/den


c------------Multiply orbital overlap by field gradient value at nucleus-------
      DO ll1=a1,a2 
      DO kk1=b1,b2   
         x1(ll1,kk1) = x1(ll1,kk1)*qq 
      ENDDO
      ENDDO 
c--------------------------------------------------------------------------------
      return
      end


