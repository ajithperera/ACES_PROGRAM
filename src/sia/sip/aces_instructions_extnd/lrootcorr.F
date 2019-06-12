c--------------------------------------------------------------------|
c   This is a superinstruction for SIAL code                         |
c   It makes correspondence between left and right eigenvectors      |
c                     using overlap matrix                           |
c   Added by Denis Bokhan, Moscow Lomonosov State University,        |
c                  MOSCOW, RUSSIAN FEDERATION                        |
c--------------------------------------------------------------------|


       subroutine lrootcorr(x,nindex, type, bval,eval,
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

c----------------------------------------------------------------
      call do_lrootcorr (x,bval(1),eval(1),type(1),
     &                     bval(2),eval(2),type(2),
     &                 x2, bval2(1),eval2(1),type2(1),
     &                     bval2(2),eval2(2),type2(2) )

c----------------------------------------------------------------
       return
       end

C=====================================================================
C=====================================================================
      subroutine do_lrootcorr(x1,a1,a2,type1,b1,b2,type2,
     &                                x2,i1,i2,type3,j1,j2,type4)
      implicit none
      include 'parallel_info.h'
      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'

      integer a1,a2,b1,b2,c1,c2,d1,d2,i1,i2,j1,j2,k1,k2
      integer i,j,k,l,a,b,l1,l2,indmax
      integer type1,type2,type3,type4,type5,type6,type7,type8
      double precision x1(a1:a2,b1:b2),
     & x2(i1:i2,j1:j2), qq, qq1, ff

c-----------Search for maximal element in the overlap array-------------
      do i = b1,b2
       indmax = 1
       do j = a1,a2
       qq = x1(indmax,i)
       qq1 = x1(j,i)
           if(dabs(qq1).gt.dabs(qq)) then
            indmax = j
           endif
       enddo
      if(me.eq.0) write(6,*) i,'////|||\\\',indmax 
       
      if (dabs(x1(indmax,i)).lt.1d-3) indmax = i 

      ff = indmax
      do k=j1,j2
      x2(i,k) = ff 
      enddo

      enddo 



      return
      end  
