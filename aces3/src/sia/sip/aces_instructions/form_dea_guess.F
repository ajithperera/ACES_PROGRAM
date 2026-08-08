 
      subroutine Form_DEA_guess(
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

      integer i, j, k,n1,n2,m1,m2,m3
      integer array_in,index,nindex_in,nindex,sind,ierr,Root
      integer block, blkndx, seg
      integer find_current_block
      integer array_in_type,array,array_type 
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
      array_in_type   = array_table(c_array_type, array_in)
C
      if (array_in_type .ne. static_array) then
         Write(6,"(a)") "@-FORM_DEA_GUESS: Array must be static."
         call abort_job()
      endif

      nindex_in  = array_table(c_nindex, array_in)

      do i = 1, nindex_in
         ind_in(i) = array_table(c_index_array1+i-1,array_in)
      enddo

      array       = op(c_op1_array)
      array_type  = array_table(c_array_type, array)
      nindex      = array_table(c_nindex, array)

      if (array_type .ne. scalar_value) then
         Write(6,"(a,a)") "@-FORM_DEA_GUESS: Second argument",
     &                  " Must be an scalar."
         call abort_job()
      endif

      sind =  array_table(c_scalar_index, array)
      Root = INT (scalar_table(sind))

      n1 = index_table(c_index_size, ind_in(1))  ! pick up length of index
      n2 = index_table(c_index_size, ind_in(2))  ! pick up length of index

c---------------------------------------------------------------------------
c   Look up each a3rray's address.
c---------------------------------------------------------------------------

      iarray_in  = get_index_from_base(address_table(array_in), x, 2)

      Call DEA_guess(x(iarray_in),Root,1,n1,1,n2)

      If (Root .Lt. 0) Then
         scalar_table(sind) = -Root 
      Else 
         scalar_table(sind) = Root 
      Endif 

      return
      end

      Subroutine DEA_guess(Vecs,Root,n1,n2,m1,m2)

      implicit none 
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      Integer n1,n2,m1,m2
      Integer Ndim,Mdim
      Integer i,j,Nroots,Row,Col,Root,Rootp1,Rootp2
      Logical Triplet 

      Double Precision Vecs((n2-n1+1),(m2-m1+1))

      If (Root .lt. 0) Then 
         Triplet = .True.
         Root    = -Root
      Endif 

      Nroots = Ea_roots 
      Ndim   = (n2-n1+1)
      Mdim   = (m2-m1+1)
      Rootp1 = Root 
      Rootp2 = Root + 1

#ifdef _DEBUG_LVL0
      if (me .eq. 0) then
      Write(6,"(a,4i,l)") "Ndim,Mdim,Rootm1,Triplet",Ndim,Mdim,
     &                 Rootm1,Rootm2,Triplet
      endif 
#endif 

      If (Triplet) Then

         do j = m1, m2
         do i = n1, n2

            Col = Rootp2
            Row = Rootp1

#ifdef _DEBUG_LVL0
      if (me .eq.0) then
      Write(6,"(a,2i)") "Row,Col", Row,Col
      endif
#endif
            If (Col .EQ. j .AND. Row .EQ. i) Then
              Vecs(Row,COl) =  1.0D0
              Vecs(COl,Row) = -1.0D0
            endif

         enddo
         enddo


      Else 

         do j = m1, m2
         do i = n1, n2

            Col =  Rootp1
            Row =  Rootp1

#ifdef _DEBUG_LVL0
      if (me .eq.0) then
      Write(6,"(a,2i)") "Row,Col", Row,Col
      endif 
#endif 
            If (Col .EQ. j .AND. Row .EQ. i) Then
              Vecs(Row,COl) =  1.0D0
              Vecs(COl,Row) =  1.0D0
            endif

         enddo
         enddo

CSSS      If (Root .eq. 1) then
CSSS          Vecs(1,1) = 1.0D0
CSSS      Elseif (Root .eq. 2) then
CSSS          Vecs(1,2) = 1.0D0
CSSS          Vecs(2,1) = 1.0D0
CSSS      Elseif (Root .eq. 3) then
CSSS          Vecs(1,1) = 1.0D0
CSSS          Vecs(2,2) = 1.0D0
CSSS      Endif 

      Endif 

#ifdef _DEBUG_LVL0 
      if (me .eq.0) then 
      write(6,"(a)") "@-DEA-GUESS vectors"
      call output(Vecs,1,n2,1,m2,n2,m2,1)
      endif
#endif 

      Return
      End

