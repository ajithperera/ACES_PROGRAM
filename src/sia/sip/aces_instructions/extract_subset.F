
      subroutine Extract_subset(
     *
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
      integer*8 iarray_in, iarray_out, ievec, daddr

      integer i, j, k,n1,n2,m1,m2
      integer array_in,array_out, index, nindex_in,nindex_out, ierr
      integer block, blkndx, seg
      integer find_current_block
      integer array_in_type, array_out_type
      integer*8 indblk1, indblk2, get_block_index
      integer stack,ind_in(mx_array_index),ind_out(mx_array_index)
      integer nsend

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
         Write(6,"(a)") "@-EXTRACT_SUBSET: Both arrays must be static."
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

       n1 = index_table(c_index_size, ind_in(1))  ! pick up length of index
       n2 = index_table(c_index_size, ind_in(2))  ! pick up length of index

       m1 = index_table(c_index_size, ind_out(1))  ! pick up length of index
       m2 = index_table(c_index_size, ind_out(2))  ! pick up length of index

c---------------------------------------------------------------------------
c   Look up each array's address.
c---------------------------------------------------------------------------

      iarray_in  = get_index_from_base(address_table(array_in), x, 2)
      iarray_out = get_index_from_base(address_table(array_out), x, 2)

      Call Strip(x(iarray_in),1,n1,1,n2,x(iarray_out),1,m1,1,m2)

      return
      end
C
      Subroutine Strip(A_in,i1,i2,j1,j2,A_out,a1,a2,b1,b2)

      implicit none 
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      Integer i,j,k,l,i1,i2,j1,j2,a1,a2,b1,b2,Ndim,Nroots 
      Integer Istart,Maxdim, Ierror,Out_dim,In_dim
      Integer Ndim1, Ndim2
      
      Parameter(Maxdim = 100)

      Double precision A_in(i1:i2,j1:j2),A_out(a1:a2,b1:b2)

      In_dim  = j2-j1+1
      Out_dim = b2-b1+1
      Ndim1   = i2-i1+1
      Ndim2   = a2-a1+1

      If (Out_dim .gt. In_dim) Then
         If (me .eq. 0) Then
         Write(6,"(a,a)") "@-Extract_subset: The sub-matrix is bigger",
     &                    " than the super-matrix"
         call abort_job()
         Endif 
      Endif 

      If (Ndim1 .NE. Ndim1) Then
         If (me .eq. 0) Then
         Write(6,"(a,a)") "@-Extract_subset: A dimension mismatch",
     &                    " of the two incomming arrays"
         Write(6,"(a,a,2(1x,i3))") "@-Extract_subset: The number of",
     &                    " rows must be identical but they are: ",
     &                      Ndim1,Ndim2
         call abort_job()
         Endif 
      Endif 

      If (Ndim1 .EQ. Nalpha_occupied .OR.
     +    Ndim1 .EQ. Nbeta_occupied) Then
          l = Int(Ip_roots)
      Elseif (Ndim1 .EQ. Nalpha_virtual .OR.
     +        Ndim1 .EQ. Nbeta_virtual) Then
          l = Int(Ea_roots)
      ENdif

      Do j = 1, l
         Do i = a1, a2
            A_out(i,j) = A_in(i,j)
         Enddo 
      Enddo 

#ifdef _DEBUG_LVL0
      if (me.eq.0) then
      k = nalpha_occupied 
CSS      k = nalpha_virtual
      l = eom_roots 
      Write(6,"(a)") "The incoming vector"
      call output(A_in,1,k,1,j2,i2,j2,1)
      Write(6,"(a)") "The outgoing vector"
      call output(A_out,1,k,1,l,a2,b2,1)
      endif 
#endif 

      Return
      End

