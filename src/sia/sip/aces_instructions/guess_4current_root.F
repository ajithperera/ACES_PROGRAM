
      subroutine Guess_4current_root (
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

      integer i, j, k,n1,n2,n3,m1,m2,m3
      integer array_in,array_out, index, nindex_in,nindex_out, ierr
      integer block, blkndx, seg
      integer find_current_block
      integer array_in_type, array, array_type
      integer*8 indblk1, indblk2, get_block_index
      integer stack,ind_in(mx_array_index),ind_out(mx_array_index)
      integer nsend,sind,nindex 
      integer comm,type
      integer fop1(mx_array_index), fop2(mx_array_index)
      integer sop1(mx_array_index), sop2(mx_array_index)
      integer sindex(mx_array_index), findex(mx_array_index)
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
      array_in_type   = array_table(c_array_type, array_in)
C
      if (array_in_type .ne. static_array) then
         Write(6,"(a,a)") "@-GUESS_4CURRENT_ROOT: First array must",
     &                  " be static."
         call abort_job()
      endif

      nindex_in  = array_table(c_nindex, array_in)

      do i = 1, nindex_in
         ind_in(i) = array_table(c_index_array1+i-1,array_in)
      enddo

      n1 = index_table(c_index_size, ind_in(1))  ! pick up length of index
      n2 = index_table(c_index_size, ind_in(2))  ! pick up length of index
      n3 = index_table(c_index_size, ind_in(3))  ! pick up length of index

      array       = op(c_op1_array)
      array_type  = array_table(c_array_type, array)
      nindex      = array_table(c_nindex, array)

      if (array_type .ne. scalar_value) then
         Write(6,"(a,a)") "@-GUESS_4CURRENT_ROOT: Second argument",
     &                  " Must be an scalar."
         call abort_job()
      endif

      sind =  array_table(c_scalar_index, array)
      type = INT (scalar_table(sind))

C---------------------------------------------------------------------------
c   Look up each a3rray's address.
c---------------------------------------------------------------------------

      iarray_in  = get_index_from_base(address_table(array_in), x, 2)

      Call Return_guess(x(iarray_in),type,1,n1,1,n2,1,n3)
 
      If (Type .Lt. 0) scalar_table(sind) = -Type 

      return
      end

      Subroutine Return_guess(Vecs_in,Root,n1,n2,m1,m2,l1,l2)

      implicit none 
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      Integer n1,n2,m1,m2,l1,l2
      Integer Ndim
      Integer Kdim, k
      INteger Root 
      Logical Triplet 

      Double Precision Vecs_in(n1:n2,m1:m2,l1:l2)

      Ndim  = (n2-n1+1)*(m2-m1+1) 
      Kdim  = (l2-l1+1)

C#ifdef _DEBUG_LVL0
      If (me .eq.0) then
         Write(6,"(a,2I)") "Root index:", Root
         Write(6,"(a)") "@-GET_GUESS_4CURRENT; Input vectors"
         Do k =1, eom_roots 
           call output(Vecs_in(1,1,k),n1,n2,m1,m2,n2,m2,1)
         Enddo 
      Endif 
C#endif 
      If (Root.le. 0.0D0) Then
          Triplet = .True.
          Root    = -Root 
      Endif 
      If (Root .Eq. 1 .and. .not. Triplet) Return

      If (Root .Eq. 1 .and. Triplet) Then
         Call Dscal(Ndim, -1.0D0, Vecs_in(1,1,1), 1)
      Else 
         Call Dcopy(Ndim, Vecs_in(1,1,Root), 1, Vecs_in(1,1,1), 1)
         If (Triplet) Call Dscal(Ndim, -1.0D0, Vecs_in(1,1,1), 1)
      Endif 

C#ifdef _DEBUG_LVL0
      If (me .eq.0) then
         Write(6,"(a)") "@-GET_GUESS_4CURRENT; Outgoing vector"
         call output(Vecs_in(1,1,1),n1,n2,m1,m2,n2,m2,1)
      Endif 
C#endif 

      Return
      End

