
      subroutine Form_RPA_guess(
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
         Write(6,"(a)") "@-FORM_RPA_GUESS: Both arrays must be static."
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
       m3 = index_table(c_index_size, ind_out(3))  ! pick up length of index

c---------------------------------------------------------------------------
c   Look up each a3rray's address.
c---------------------------------------------------------------------------

      iarray_in  = get_index_from_base(address_table(array_in), x, 2)
      iarray_out = get_index_from_base(address_table(array_out), x, 2)

      Call RPA_guess(x(iarray_in),x(iarray_out),1,n1,1,n2,1,m1,1,m2,
     &               1,m3)
      return
      end

      Subroutine RPA_guess(Diags,Vecs,n1,n2,m1,m2,k1,k2,l1,l2,p1,p2)

      implicit none 
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      Integer n1,n2,m1,m2,k1,k2,l1,l2,p1,p2
      Integer Ndim,Itmp,Jtmp,Mdim,Ibgn
      Integer i,j,k,Nroots,Max_roots,ondx,vndx
      Parameter(Max_roots=30)

      Double Precision Diags((n2-n1+1),(m2-m1+1))
      Double Precision Vecs(k1:k2,l1:l2,p1:p2)
      Double Precision Dtmp((n2-n1+1),(m2-m1+1))
      Double Precision Xmin,Tol
      Integer O_Label(Max_roots)
      Integer V_Label(Max_roots)

      Nroots = Eom_roots 
      Ndim   = (n2-n1+1)
      Mdim   = (m2-m1+1)
      Tol    = 1.0D-12

      Call Dcopy(Ndim*Mdim, Diags, 1, Dtmp, 1)

#ifdef _DEBUG_LVL0 
      if (me .eq. 0) then
      write(6,"(a,3i3)") "Input dimensions:", k2,l2,p2
      write(6,"(a)") "@-RPA-GUESS @-entry, the diagonals"
      call output(Diags,1,n2,1,m2,n2,m2,1)
      endif 
#endif 

      Do k = 1, Nroots 

          Xmin = Diags(1,1)

          Do j= 1, Mdim
             If (j .Eq. 1) Then
                 Ibgn = 2 
             Else
                 Ibgn = 1
             Endif 
             Do i=Ibgn, Ndim 
                If (Dtmp(i,j) .Gt. Tol) Then
                If (Dtmp(I,J) .LT. Xmin) Then
                    Xmin  = Dtmp(I,J)
                    Itmp  = I
                    Jtmp  = J
                Endif 
                Endif 
             Enddo 
          Enddo

          O_label(k) = Jtmp
          V_label(k) = Itmp 
          Dtmp(Itmp,Jtmp) = 1.0D10
      Enddo 

#ifdef _DEBUG_LVL0 
      if (me .eq.0) then
      write(6,"(a)") "@-RPA-GUESS The position of smallest digonals"
      write(6,"(6(1x,i3))") (O_label(i),i=1,Nroots)
      write(6,"(6(1x,i3))") (V_label(i),i=1,Nroots)
      endif
#endif 
    
      Do k = 1, Nroots 
         Ondx = O_label(k)
         Vndx = V_label(k)
         Vecs(Vndx,Ondx,k) = 1.0D0 
      Enddo

#ifdef _DEBUG_LVL0 
      if (me .eq.0) then 
      write(6,"(a)") "@-RPA-GUESS vectors"
      Do k = 1, Nroots
          call output(Vecs(1,1,k),k1,k2,l1,l2,k2,l2,1)
      Enddo 
      endif
#endif 

      Return
      End

