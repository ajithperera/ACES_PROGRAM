c--------------------------------------------------------------------|
c   This is a superinstruction for SIAL code                         |
c   It prints out EOM-CC roots to the summary.out file               |
c   Added by Denis Bokhan, Moscow Lomonosov State University,        |
c                  MOSCOW, RUSSIAN FEDERATION                        |
c--------------------------------------------------------------------|


      subroutine set_tr_sin(x,nindex, type, bval,eval,
     & bdim, edim, x2, nindex2, type2, bval2, eval2, bdim2, edim2)

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'

      double precision x(*), ntr, nsin
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)
      double precision x2(*)
      integer nindex2, type2(*), bval2(*), eval2(*)
      integer bdim2(*), edim2(*)
c----------------------------------------------------------------
      call do_set_tr_sin(x,x2)
c----------------------------------------------------------------
      return
      end
c==================================================================
c==================================================================
      subroutine do_set_tr_sin(x1,x2)

            implicit none
      include 'parallel_info.h'
      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
      include 'spin_orb.h'

      integer a1,a2,b1,b2,c1,c2,d1,d2,i1,i2,j1,j2,k1,k2
      integer i,j,k,l,a,b,l1,l2
      integer type1,type2,type3,type4,type5,type6,type7,type8
      double precision x1, x2

      n_tr = x1
      n_si = x2


      return
      end
