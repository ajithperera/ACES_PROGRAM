
      subroutine Form_norm_with_bmat(
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
      If (me .eq. 0) Then
         Write(6,"(a,a)") "@-FORM_NORM_WITH_BMAT: Both arrays must",
     &                  " be static."
      endif 
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

      Call Get_norms(x(iarray_in),1,n1,1,n2,x(iarray_out),1,m1,1,m2)

      return
      end
C
      Subroutine Get_norms(Evecs,i1,i2,j1,j2,Pmat,a1,a2,b1,b2)

      implicit none 
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'ex_mask_vars.h'

      Integer i,l,i1,i2,j1,j2,a1,a2,b1,b2,Ndim,Ndim1,Ndim2,Nroots 
      Integer Istart,Ierror,k

      Double precision Tmp(Max_pdim,Max_pdim)
      Double precision Evecs(i1:i2,j1:j2),Pmat(a1:a2,b1:b2)
      Double precision Ddot 

      Ndim1 = (i2 -i1) + 1
      Ndim2 = (b2 -b1) + 1

      If (Ndim1 .NE. Ndim2) Then
         If (Me .eq. 0) Then
         Write(6,"(2a)")" The two incomming matrices are not",
     +                  " not allowed to have two different dimensions"
         Endif 
         Call Abort_job()
      Endif 

      Ndim = Ndim1

      If (Ndim .Gt. Max_pdim) Then
         If (Me .eq. 0) Then
         Write(6,"(2a,i3)") "@-FORM_NORM_WITH_BAMT: Exceeded the max.",
     +                      " dimension allowed, the maximum is ", 
     +                        Max_pdim
         Endif 
         Call abort_job()
      Endif 

#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a)") "@-Entry to GET_NORMS the two input matrices"
      call output(Evecs,1,Ndim,1,Ndim,Ndim,Ndim,1)
      Write(6,*)
      call output(Pmat,1,Ndim,1,Ndim,Ndim,Ndim,1)
      call zero
      Endif
#endif
      Call Xgemm("N","N",Ndim,Ndim,Ndim,1.0D0,Pmat,Ndim,Evecs,Ndim,
     +            0.0D0,Tmp,Max_pdim)
      Do I = 1, Ndim 
         P_ovlp(I) = Ddot(Ndim*Ndim,Tmp(1,I),1,Evecs(1,I),1)
      Enddo

#ifdef _DEBUG_LVL0
      if (me.eq.0) then
      k =0
      do i=1,Ndim
      If (P_ovlp(i) .ne. 0.0D0) k=k+1   
      enddo 
      Write(6,"(a)") "The C^(t)PC array"
      Write(6,"(5(1x,F15.10))") (P_ovlp(i),i=1,k)
      endif 
#endif 

      Return
      End

