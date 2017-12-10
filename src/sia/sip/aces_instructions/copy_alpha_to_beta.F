
      subroutine Copy_Alpha_to_Beta(
     *                      array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)

c---------------------------------------------------------------|
      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'mpif.h'
      include 'parallel_info.h'
      include 'int_gen_parms.h'
#ifdef ALTIX
      include 'sheap.h'
#endif
c---------------------------------------------------------------|

      integer narray_table, nindex_table, nsegment_table,
     *        nblock_map_table
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table)
      integer*8 iarray_in, iarray_out, ievec, daddr

      integer i, j, k,n1,n2,n3,n4,m1,m2,m3
      integer array_in,array_out, index, nindex_in,nindex_out, ierr
      integer block, blkndx, seg
      integer find_current_block
      integer array_in_type, array_out_type
      integer*8 indblk1, indblk2, get_block_index
      integer stack,ind_in(mx_array_index),ind_out(mx_array_index)
      integer nsend
      integer comm
      integer type(mx_array_index)
      integer na1, na2, ni1, ni2
      integer Idim,Jdim,Kdim,Ldim,Mdim,Ndim,IJKdim,LMNdim
      integer*8 addr, get_index_from_base
      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

#ifdef ALTIX
      dptr = dshptr
#endif
      array_in = op(c_result_array)
      array_out = op(c_op1_array)

      array_in_type   = array_table(c_array_type, array_in)
      array_out_type  = array_table(c_array_type, array_out)
C
      if (array_in_type .ne. static_array .or.array_out_type 
     &                  .ne. static_array) then
         Write(6,"(a,a)") "@-COPY_ALPHA_TO_BETA: Both arrays",
     &                    " must be static."
         call abort_job()
      endif

      nindex_in  = array_table(c_nindex, array_in)
      nindex_out = array_table(c_nindex,array_out)

      do i = 1, nindex_in
         ind_in(i) = array_table(c_index_array1+i-1,array_in)
      enddo

      do j = 1, nindex_out
         ind_out(j) = array_table(c_index_array1+j-1,array_out)
      enddo

       n1 = index_table(c_index_size, ind_in(1))  
       n2 = index_table(c_index_size, ind_in(2)) 
       n3 = index_table(c_index_size, ind_in(3)) 

       m1 = index_table(c_index_size, ind_out(1))  
       m2 = index_table(c_index_size, ind_out(2)) 
       m3 = index_table(c_index_size, ind_out(3)) 

c---------------------------------------------------------------------------
c   Look up each array's address.
c---------------------------------------------------------------------------

      iarray_in  = get_index_from_base(address_table(array_in), x, 2)
      iarray_out = get_index_from_base(address_table(array_out), x, 2)

      Idim = n1
      Jdim = n2
      Kdim = n3
   
      Ldim = m1
      Mdim = m2
      Ndim = m3
     
      IJKdim = n1*n2*n3
      LMNdim = m1*m2*n3

#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a,8i4)") "@-COPY_ALPHA_TO_BETA: variables",
     +           Idim,Jdim,Kdim,Ldim,Mdim,Ndim,IJKdim,LMNdim
      Endif 
#endif 

      Call Copy_A_to_B(x(iarray_in),IJKdim,x(iarray_out),LMNdim)

      Return
      End

      Subroutine Copy_A_to_B(A,IJKdim,B,LMNdim)

      implicit none 
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      Integer IJKdim,LMNdim

      Double Precision A(IJKdim) 
      Double Precision B(LMNdim)

      If (IJKdim .NE. LMNdim) Then
         If (me. eq.0) then
         Write(6,"(a,a)") "@-COPY_ALPHA_TO_BETA: the dimensions ",
     +                    " of both arrays must be identical to copy"
         Endif 
         Call abort_job()
      Endif 

      Call Dcopy(IJKdim, A, 1, B, 1)

      Return
      End

