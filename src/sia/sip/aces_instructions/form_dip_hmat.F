
      subroutine Form_DIP_Hmat(
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

      integer fop1(mx_array_index), fop2(mx_array_index)
      integer sop1(mx_array_index), sop2(mx_array_index)
      integer sindex(mx_array_index), findex(mx_array_index)
      integer type(mx_array_index)
      integer na1, na2, ni1, ni2
      integer Idim,Jdim,Kdim,Ldim,IJDim,KLDim,MMDim
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
         Write(6,"(a)") "@-FORM_DIP_HMAT: Both arrays must be static."
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
       n4 = index_table(c_index_size, ind_in(4)) 

       m1 = index_table(c_index_size, ind_out(1))  
       m2 = index_table(c_index_size, ind_out(2)) 

c---------------------------------------------------------------------------
c   Look up each array's address.
c---------------------------------------------------------------------------

      iarray_in  = get_index_from_base(address_table(array_in), x, 2)
      iarray_out = get_index_from_base(address_table(array_out), x, 2)

      Idim = n1
      Jdim = n2
      Kdim = n3
      Ldim = n4
     
#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a,9i4)") "@-Entry to pack4d_to_2D variables",
     +                   Idim,jdim,kdim,ldim,m1,m2
      Endif 
#endif 
      IJdim= n1*n2
      KLdim= n3*n4
      MMdim= m1*m2
  
      If (IJdim .Ne. KLdim .And. IJdim .Ne. MMdim .And.
     +    KLdim .Ne. MMdim) Then
          Write(6,"(a)") "The dimensions of the matirces do not match"
          Call abort_job()
      Endif 

      Call Form_Hmat(x(iarray_in),Idim,Jdim,Kdim,Ldim,
     +               x(iarray_out)) 

      Return
      End

      Subroutine Form_Hmat(Hbar_ijij,Idim,Jdim,Kdim,Ldim,Hbar_ij)

      implicit none 
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      Integer I,J,K,L,P,Idim,Jdim,Kdim,Ldim

      Double Precision Hbar_ijij(Idim,Jdim,kdim,Ldim) 
      Double Precision Hbar_ij(Idim,Jdim)

#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a,3i4)") "@-pack4d_to_2D variables;before packing",
     +                    Idim,jdim,kdim,ldim
      Endif 
#endif 
      
      DO J =1, Jdim
      DO I =1, Idim
         DO L =1, Ldim
         DO K =1, Kdim
            IF (J .Eq. L)
     +       Hbar_IJIJ(I,J,K,L) = Hbar_IJIJ(I,J,K,L) + Hbar_ij(I,K)
         ENDDO
         ENDDO
      ENDDO
      ENDDO

      DO J =1, Jdim
      DO I =1, Idim
         DO L =1, Ldim
         DO K =1, Kdim
            IF (I .Eq. K)
     +       Hbar_IJIJ(I,J,K,L) = Hbar_IJIJ(I,J,K,L) + Hbar_ij(J,L)
         ENDDO
         ENDDO
      ENDDO
      ENDDO

      Return
      End
