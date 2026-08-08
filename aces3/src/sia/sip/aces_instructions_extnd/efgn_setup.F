c---------------------------------------------------------------|
c  This subroutine is a superinstruction for the SIAL code      |
c  Added by Denis Bokhan, Moscow Lomonosov State University,    |
c              Moscow, RUSSIAN FEDERATION                       |
c---------------------------------------------------------------|
      subroutine efgn_setup(x,nindex, type, bval,eval, bdim, edim) 


      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)
c----------------------------------------------------------------

      call do_efgn_setup(x)
c----------------------------------------------------------------
      return
      end

C----------------------------------------------------------------
C----------------------------------------------------------------

      subroutine do_efgn_setup(x1)
      implicit none
      include 'parallel_info.h'
      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
C      include 'F12.h'
      include 'efgn.h'

      integer a1,a2,b1,b2,c1,c2,d1,d2,i1,i2,j1,j2,k1,k2
      integer i,j,k,l,a,b,l1,l2,g,kk1,ll1
      integer type1,type2,type3,type4,iatom
      double precision x1

      ll1 = x1
      atindex = ll1 
c      if(me.eq.0) write(6,*) 'atindedx=',atindex 

c--------------------------------------------------------------------------------
      return
      end


