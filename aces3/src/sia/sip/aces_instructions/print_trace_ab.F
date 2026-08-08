
      subroutine Print_trace_ab(
     *                 array_table, narray_table,
     *                 index_table,
     *                 nindex_table, segment_table, nsegment_table,
     *                 block_map_table, nblock_map_table,
     *                 scalar_table, nscalar_table,
     *                 address_table, op)

c---------------------------------------------------------------|
      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'mpif.h'
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
      integer*8 iarray_in1,iarray_in2,ievec,daddr,iarray_ca

      integer i, j, k,n1,n2,m1,m2,l1,l2
      integer array_in,array_out,index,nindex_in,nindex_out, ierr
      integer nindex_ca,array_ca
      integer block, blkndx, seg
      integer find_current_block
      integer array_in_type, array_out_type,array_ca_type
      integer*8 indblk1, indblk2, get_block_index
      integer stack,ind_in(mx_array_index),ind_out(mx_array_index)
      integer ind_ca(mx_array_index)
      integer nsend, Ndim

      integer comm

      integer fop1(mx_array_index), fop2(mx_array_index)
      integer sop1(mx_array_index), sop2(mx_array_index)
      integer sindex(mx_array_index), findex(mx_array_index)
      integer type(mx_array_index)
      integer na1, na2, ni1, ni2
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
         Write(6,"(a)") "@-PRINT_TRACE_AB Both arrays must be static."
         call abort_job()
      endif

      nindex_in  = array_table(c_nindex,array_in)
      nindex_out = array_table(c_nindex,array_out)

      do i = 1, nindex_in
         ind_in(i) = array_table(c_index_array1+i-1,array_in)
      enddo

      do j = 1, nindex_out
         ind_out(j) = array_table(c_index_array1+j-1,array_out)
      enddo

       n1 = index_table(c_index_size, ind_in(1))  ! pick up length of index
       n2 = index_table(c_index_size, ind_in(2))  ! pick up length of index

       m1 = index_table(c_index_size, ind_out(1))  ! pick up length of index
       m2 = index_table(c_index_size, ind_out(2))  ! pick up length of index

c---------------------------------------------------------------------------
c   Look up each array's address.
c---------------------------------------------------------------------------

      iarray_in1 = get_index_from_base(address_table(array_in), x, 2)
      iarray_in2 = get_index_from_base(address_table(array_out), x, 2)

      Ndim = N1

      Call Print_trace(x(iarray_in1),1,n1,1,n2,x(iarray_in2),1,m1,
     &                    1,m2,Ndim)

      Return
      End
C
      Subroutine Print_trace(H,i1,i2,j1,j2,D,a1,a2,b1,b2,Ndim)

      implicit none 
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      Integer i1,i2,j1,j2,a1,a2,b1,b2,i,j
      Integer Ndim 

      Double precision H(i1:i2,j1:j2),D(a1:a2,b1:b2)
      Double precision Suma 

      Suma = 0.0D0 
      Do i=1, Ndim
      Do j=1, Ndim
         Suma = Suma + H(i,j) * D(j,i)
      Enddo 
      Enddo 

      If (me.eq.0) Then
         write(6,*) 
         write(6,"(a,1X,F15.9)") " Trace(h*D) =",suma
      Endif 

      Return
      End
