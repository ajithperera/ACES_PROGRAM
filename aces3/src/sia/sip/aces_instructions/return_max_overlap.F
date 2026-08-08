
      subroutine Return_max_overlap(
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
      include 'parallel_info.h'
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
      integer array_in,index, nindex_in,nindex_out, ierr
      integer block, blkndx, seg
      integer find_current_block
      integer array_in_type, array,array_type
      integer*8 indblk1, indblk2, get_block_index
      integer stack,ind_in(mx_array_index),ind_out(mx_array_index)
      integer nsend,sind,nindex 

      integer comm

      integer fop1(mx_array_index), fop2(mx_array_index)
      integer sop1(mx_array_index), sop2(mx_array_index)
      integer sindex(mx_array_index), findex(mx_array_index)
      integer type
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
         Write(6,"(a)") "@-RETURN_MAX_OVERLAP: Array must be static."
         call abort_job()
      endif

      nindex_in  = array_table(c_nindex, array_in)

      do i = 1, nindex_in
         ind_in(i) = array_table(c_index_array1+i-1,array_in)
      enddo

       n1 = index_table(c_index_size, ind_in(1))  ! pick up length of index

c   Look up each array's address.

      iarray_in  = get_index_from_base(address_table(array_in), x, 2)

      array       = op(c_op1_array)
      array_type  = array_table(c_array_type, array)
      nindex      = array_table(c_nindex, array)

      if (array_type .ne. scalar_value) then
         Write(6,"(a)") "@-RETURN_MAX_OVERLAP: Must be an scalar."
         call abort_job()
      endif

      sind =  array_table(c_scalar_index, array)
      type = INT (scalar_table(sind))

CSSS      If (me.eq.0) Write(6,*) "@-Return_max_overlap: in = ", Type
      Call Overlap_max(x(iarray_in),1,n1,Type)

CSSS      If (me.eq.0) Write(6,*) "@-Return_max_overlap: out = ", Type

      scalar_table(sind) = Dble(Type)

      return
      end
C
      Subroutine Overlap_max(Overlap,i1,i2,Over_max)

      implicit none 
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      Integer i,j,k,l,i1,i2
      Integer Idim
      Integer Over_max 

      Double Precision Overlap(i1:i2)

      Idim = i2 -i1 + 1

#ifdef _DEBUG_LVL0 
      if (me .eq.0) then
      l=0
      do i=1,idim
         if (Dabs(Overlap(i)) .ne. 0.0d0) l=l+1
      enddo 
      write(6,"(a)") "@-RETURN_MAX_OVERLP The overlap array"
      Write(6,"(6(1x,F10.5))") (Dabs(Overlap(i)),i=1,l)
      endif
#endif 
      Over_max = 1 
      Do I = 1, Idim
        if (Dabs(Overlap(i)) .gt. Dabs(Overlap(over_max))) then
            Over_max = i
        Endif
      Enddo 
#ifdef _DEBUG_LVL0 
      if (me .eq.0) then
      write(6,"(a,i3)") "The position of the maximum overalp:",
     +                   Over_max
      endif 
#endif 

      Return
      End

