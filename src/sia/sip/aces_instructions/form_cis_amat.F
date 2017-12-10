
      subroutine Form_CIS_Amat(
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
      integer stack,ind_in_bgn(mx_array_index)
      integer ind_in_end(mx_array_index),ind(mx_array_index)
      Integer ind_out(mx_array_index),n(mx_array_index)
      integer nsend

      integer comm

      integer fop1(mx_array_index), fop2(mx_array_index)
      integer sop1(mx_array_index), sop2(mx_array_index)
      integer sindex(mx_array_index), findex(mx_array_index)
      integer type(mx_array_index)
      integer na1, na2, ni1, ni2
      integer Idim,Jdim,Kdim,Ldim,Mdim,Ndim
      integer Ibgn,Jbgn,Kbgn,Lbgn,Iend,Jend,Kend,Lend
      integer IJdim,KLdim,MNdim
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
      if (array_in_type .ne. temp_array .or.array_out_type 
     &                  .ne. static_array) then
         Write(6,"(a,a)") "@-FORM_CIS_AMAT: Arrays must be temp and",
     &                  "  static in that order"
         call abort_job()
      endif

      nindex_in  = array_table(c_nindex, array_in)
      nindex_out = array_table(c_nindex, array_out)

      If (nindex_in .NE. 4 .OR. nindex_out .NE. 2) Then
          Write(6,"(a,a)") "@-FORM_CIS_AMAT: The arrays must be",
     &                     " rank 4 and 2 respectively"
          call abort_job()
      Endif 

      block  = find_current_block(array_in, array_table(1,array_in),
     *                            index_table, nindex_table,
     *                            segment_table, nsegment_table,
     *                            block_map_table, blkndx)
      stack = array_table(c_array_stack,array_in)
      iarray_in = get_block_index(array_in, block, stack, blkndx, 
     *                            x, .true.)

      do i = 1, nindex_in
         ind(i) = array_table(c_index_array1+i-1,array_in)
         n(i)   = index_table(c_index_size, ind(i))  ! pick up length of index
         seg    = index_table(c_current_seg,ind(i))

         call get_index_segment(ind(i), seg, segment_table,
     *                          nsegment_table, index_table,
     *                          nindex_table, ind_in_bgn(i),
     *                          ind_in_end(i))
      enddo

      do j = 1, nindex_out
         ind_out(j) = array_table(c_index_array1+j-1,array_out)
      enddo
      iarray_out = get_index_from_base(address_table(array_out), x, 2)

       m1 = index_table(c_index_size, ind_out(1))  
       m2 = index_table(c_index_size, ind_out(2)) 

c---------------------------------------------------------------------------
c   Look up each array's address.
c---------------------------------------------------------------------------
      Kbgn = ind_in_bgn(1)
      Lbgn = ind_in_bgn(2)
      Ibgn = ind_in_bgn(3) 
      Jbgn = ind_in_bgn(4)
      Kend = ind_in_end(1)
      Lend = ind_in_end(2)
      Iend = ind_in_end(3)
      Jend = ind_in_end(4)

      Idim = Iend-Ibgn+1
      Jdim = Jend-Jbgn+1
      Kdim = Kend-Kbgn+1
      Ldim = Lend-Lbgn+1

      Mdim = m1
      Ndim = m2
     
#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,*)
      Write(6,"(a,6(1x,i4))") "@-Entry to Form_CIS_Amat variables",
     +                   Idim,jdim,kdim,ldim,m1,m2
      Endif 
#endif 
      IJdim= Idim * Jdim
      KLdim= Kdim * Ldim
      MNdim= m1*m2

      If (MNdim .LT. IJDim*KLdim) Then
         Write(6,"(3a,i8,a,i8)")
     +                   "The space allocated for the CIS A matrix",
     +                   " from the SIAL code is smaller than what is",
     +                   " required, Allocated :", MNDim, "Requested:",
     +                   IJDim*KLdim
         Call abort_job()
      Endif 
#ifdef _ASSERT 
      If (IJdim .Ne. KLdim) Then 
          Write(6,"(a,a)") "The row and columns dimensions of", 
     +                     " the first matrix do not match"
          Call abort_job()
      Endif 
#endif 

      Call Form_Amat(x(iarray_in),Idim,Jdim,Kdim,Ldim,Mdim,Ndim,
     +               Ibgn,Jbgn,Kbgn,Lbgn,Iend,Jend,Kend,Lend,
     +               x(iarray_out)) 

      Return
      End

      Subroutine Form_Amat(A_iaia,Idim,Jdim,Kdim,Ldim,Mdim,Ndim,
     +                     Ibgn,Jbgn,Kbgn,Lbgn,Iend,Jend,Kend,Lend,
     +                     A_ai)

      implicit none 
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'aces3_vars.h'

      Integer I,J,K,L,P,Idim,Jdim,Kdim,Ldim,Mdim,Ndim,Nai
      Integer IJ,KL,Ioff,Joff,n
      Integer nocc_alpha,nvrt_alpha,nocc_beta,nvrt_beta
      Integer nocc_nvrt_a,nocc_nvrt_b
      integer Ibgn,Jbgn,Kbgn,Lbgn,Iend,Jend,Kend,Lend
      integer Vrbgn,Vlbgn 

      Double Precision A_iaia(Kbgn:Kend,Lbgn:Lend,Ibgn:Iend,Jbgn:Jend) 
      Double Precision A_ai(Mdim,Ndim)

#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,*)
      Write(6,"(2a,6i4)") "@-Form_Amat variables before forming", 
     +                    " the CIS matrix:", Idim,jdim,kdim,ldim,
     +                      mdim,ndim
      Endif 
#endif 
      
      nocc_alpha = nalpha_occupied
      nvrt_alpha = nalpha_virtual 
      nocc_beta  = nbeta_occupied
      nvrt_beta  = nbeta_virtual 

      nocc_nvrt_a = nocc_alpha * nvrt_alpha 
      nocc_nvrt_b = nocc_beta  * nvrt_beta 
      
      Nai = nocc_nvrt_a + nocc_nvrt_b 

#ifdef _DEBUG_LVL0
      If (me .eq. 0) then
         Write(6,"(a,a)") "@-Entry to Form_amat internal variables:", 
     +                      " nocc_alpha,nvrt_alpha"
         Write(6,"(a)")     "nocc_beta,nvrt_beta: "
         Write(6,"(4i3)") nocc_alpha,nvrt_alpha,nocc_beta,nvrt_beta
         Write(6,"(a,3i2)") "The spin-types:",AAAA,BBBB,ABAB
         Write(6,"(a,8i3)") "Ibgn,Jbgn....Iend,Jend", 
     +                       Ibgn,Jbgn,Kbgn,Lbgn,Iend,Jend,Kend,Lend
      Endif
#endif

      If (AAAA .EQ. 1) Then 

C Alpha,Alpha block

         DO J =Jbgn, Jend
         DO I =Ibgn, Iend
            Vrbgn = I - nocc_alpha
            IJ    = (J-1) * nvrt_alpha  + Vrbgn
            KL = 0
            DO L =Lbgn, Lend
            DO K =Kbgn, Kend
               Vlbgn = K - nocc_alpha 
               KL    = (L-1) * nvrt_alpha + Vlbgn 
               A_ai(KL,IJ) = A_iaia(K,L,I,J) + A_ai(KL,IJ)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      Endif 

      If (BBBB .EQ. 1) Then 

C Beta/Beta block

         Ioff = nocc_nvrt_a 
         Joff = nocc_nvrt_a 

         DO J =Jbgn, Jend
         DO I =Ibgn, Iend
            Vrbgn = I - nocc_beta
            IJ    = (j-1) *  nvrt_beta + Vrbgn
            DO L =Lbgn, Lend
            DO K =Kbgn, Kend
               Vlbgn = K - nocc_beta 
               KL    =  (L-1) * nvrt_beta + Vlbgn
               A_ai(KL+Ioff,IJ+Joff) = A_iaia(K,L,I,J) + 
     +                                 A_ai(KL+Ioff,IJ+Joff)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      Endif 

      If (ABAB .EQ. 1) Then 

C Alpha/Beta and Beta/Alpha blocks

         Joff = nocc_nvrt_b
         Ioff = nocc_nvrt_a

         DO J =Jbgn, Jend
         DO I =Ibgn, Iend
            Vrbgn = I - nocc_beta 
            IJ = (J-1) * nvrt_beta + Vrbgn
            DO L =Lbgn, Lend
            DO K =Kbgn, Kend
               Vlbgn = K - nocc_alpha 
               KL    = (L-1) * Nvrt_alpha + Vlbgn

               A_ai(KL,IJ+Joff) = A_iaia(K,L,I,J) + A_ai(KL,IJ+Joff)
               A_ai(KL+Ioff,IJ) = A_iaia(K,L,I,J) + A_ai(KL+Ioff,IJ)
            ENDDO
            ENDDO
         ENDDO
         ENDDO

      Endif 
#ifdef _DEBUG_LVL0
      If (me .eq. 0) then
         Write(6,"(a)") " The CIS matrix"
         n= nocc_nvrt_a
         If (AAAA .EQ. 1) Then 
         Write(6,"(a)")  "The CIS matrix: AA block"
         call output(A_ai,1,N,1,N,Ndim,Ndim,1)
         ElseIf (BBBB .EQ. 1) then
         Write(6,"(a)")  "The CIS matrix: BB block"
         call output(A_ai(ioff+1,joff+1),1,N,1,N,Ndim,Ndim,1)
         ElseIf (ABAB .EQ. 1) then
         Write(6,"(a)")  "The CIS matrix: AB block"
         call output(A_ai(1,1+Joff),1,N,1,N,Ndim,Ndim,1)
         Write(6,"(a)")  "The CIS matrix: BA block"
         call output(A_ai(1+Ioff,1),1,N,1,N,Ndim,Ndim,1)
         Endif 
      Endif
#endif
      Return
      End
