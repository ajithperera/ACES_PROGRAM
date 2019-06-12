c---------------------------------------------------------------|
c  This subroutine is a superinstruction for the SIAL code      |
c  It calculates inverted U(+) and U(-) matrices for STEOM      |
c  calculations.                                                | 
c  Added by Denis Bokhan, Moscow Lomonosov State University,    |
c              Moscow, RUSSIAN FEDERATION                       |
c---------------------------------------------------------------|
      subroutine steom_u_inv (array_table,
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)


c-------------------------------------------------------------------------- 
      implicit none
      include 'mpif.h'
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'pst_functions.h'
      include 'daveom.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

c---------------Variables--------------------------------------------
      integer narray_table, nindex_table, nsegment_table,
     *        nblock_map_table, array
      integer op(loptable_entry),sind
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table), get_index_from_base
      integer*8 get_block_data_index
      integer find_current_block, block

      integer array_type, nindex, inda(2), indb(2)

      logical direct_flag

      integer i, j, k,qq,qq1,is_it_worker(max_procsx)
      integer nind, nsend,ind1,w_worker(max_procsx)
      integer ierr, handle,our_workers
      integer blk, blkndx, maxblk
      integer*8 iarray, farray, indblk, get_block_index
      integer*8 integral_scr,daddr,daddr1,iarray1
      integer seg, val1(mx_array_index), val2(mx_array_index)
      integer msg(len_sip_server_message)
      integer flag,ind(mx_array_index)
      integer stack,g,n1,n2,m1,m2

      double precision x(1),int_type_d,dind1
      double precision dbuf(1), temp
      integer ibuf(1)
      integer int_type
c----------------------------------------------------------------------
#ifdef ALTIX
      pointer (dptr, dbuf)
      pointer (xptr, x)
      pointer (iptr, ibuf)
#else
      common dbuf, x
      equivalence (ibuf(1), dbuf(1))
#endif

      common /d2int_com/jatom, jx, jcenter
      integer jatom, jx, jcenter

      array = op(c_result_array)

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
      xptr = dshptr
#endif

c-------------------S(+/-) input array---------------------------------
      array      = op(c_result_array)

      array_type = array_table(c_array_type, array)
      nindex     = array_table(c_nindex, array)

      daddr = address_table(array)
      iarray  = get_index_from_base(daddr, x, 2)


      do i = 1, nindex
         ind(i) = array_table(c_index_array1+i-1,array)
      enddo

      n1 = index_table(c_index_size, ind(1))  ! pick up length of index
      n2 = index_table(c_index_size, ind(2))  ! pick up length of index
c-------------------U_inv(+/-) output array------------------------------
      array = op(c_op1_array)

      array_type = array_table(c_array_type, array)
      nindex     = array_table(c_nindex, array)

      daddr1 = address_table(array)
      iarray1  = get_index_from_base(daddr1, x, 2)

      do i = 1, nindex
         ind(i) = array_table(c_index_array1+i-1,array)
      enddo

      m1 = index_table(c_index_size, ind(1))  ! pick up length of index
      m2 = index_table(c_index_size, ind(2))  ! pick up length of index
c-------------------------------------------------------------------------
      call do_u_inv(x(iarray),n1,n2,x(iarray1),m1,m2,ip_roots,ea_roots)
c-------------------------------------------------------------------------

      return
      end

c#########################################################################
c#########################################################################         

      subroutine do_u_inv(S_aux,n1,n2,U_inv,m1,m2,n_ip_roots,n_ea_roots)
      implicit none
      include 'parallel_info.h'
      include 'daveom.h'
      include 'int_gen_parms.h'

      integer i,j,n1,n2,m1,m2,izz,k,info,l
      integer n_ip_roots,n_ea_roots
      double precision S_aux(n1,n2), U_inv(m1,m2),qq
      double precision U_ip(n_ip_roots,n_ip_roots)
      double precision U_ea(n_ea_roots,n_ea_roots)
      double precision U_ip_aux_inv(n_ip_roots,n_ip_roots)
      double precision U_ea_aux_inv(n_ea_roots,n_ea_roots)
      integer ipiv(n_ip_roots),ipiv2(n_ea_roots,n_ea_roots)
      double precision test_ip(n_ip_roots,n_ip_roots),
     & test_ea(n_ea_roots,n_ea_roots)

c-----------Active occupied part---------------------------
c       write(6,*) n1,n2,'===>>',m1,m2

       n_ip_roots = ip_roots

       IF(RL.eq.1) then

       do i=1,n_ip_roots
       do j=1,n_ip_roots
       U_ip_aux_inv(i,j) = 0.d0
       if(i.eq.j) U_ip_aux_inv(i,j) = 1.d0
       enddo
       enddo 


       do j=1,n_ip_roots
       do i=1,n_ip_roots
       k = n1-n_ip_roots+i
       U_ip(i,j) = S_aux(k,j) 
       enddo
       enddo
        

       call dgesv(n_ip_roots,n_ip_roots,U_ip,n_ip_roots,ipiv,
     & U_ip_aux_inv,n_ip_roots,info)

       if(me.eq.0) write(6,*) 'info=',info

         do i=1,m2
          do j=1,m1
          k = m1-n_ip_roots+1
          l = n_ip_roots-m1+j
           if(i.gt.n_ip_roots.or.j.lt.k) then
           U_inv(j,i) = 0.d0
           else
           U_inv(j,i) = U_ip_aux_inv(i,l)
c           if (me.eq.0) write(6,*) U_inv(j,i),'>>>>>>>>',j,i 
           endif
          enddo
         enddo

      endif
c--------Active virtual part-------------------------------

       IF(RL.eq.2) then

       do i=1,n_ea_roots
       do j=1,n_ea_roots
       U_ea_aux_inv(i,j) = 0.d0
       if(i.eq.j) U_ea_aux_inv(i,j) = 1.d0
       enddo
       enddo

       do j=1,n_ea_roots
       do i=1,n_ea_roots
       U_ea(i,j) = S_aux(i,j) 
       enddo
       enddo
      

       call dgesv(n_ea_roots,n_ea_roots,U_ea,n_ea_roots,ipiv2,
     & U_ea_aux_inv, n_ea_roots,info)

       if(me.eq.0) write(6,*) 'info=',info

         do i=1,m2
          do j=1,m1

           if(i.gt.n_ea_roots.or.j.gt.n_ea_roots) then
           U_inv(j,i) = 0.d0
           else
           U_inv(j,i) = U_ea_aux_inv(i,j)
           endif
          enddo
         enddo

       endif

c         do i=1,m1
c          do j=1,m2
c          if (me.eq.0) write(6,*) U_inv(i,j),'>>>>>>>>',i,j
c          enddo
c         enddo

      

c        do i=1,n_ea_roots
c        do j=1,n_ea_roots
c        qq = 0.d0
c        do k=1,n_ea_roots
c        qq = qq + test_ea(i,k)*U_ea_aux_inv(k,j)
c        enddo
c        if(me.eq.0) write(6,*) 'qq=',qq,'^^>',i,j
c        enddo 
c        enddo



      return
      end 
