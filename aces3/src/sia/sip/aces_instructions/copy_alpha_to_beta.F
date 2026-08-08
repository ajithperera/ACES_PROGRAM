C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
      subroutine Copy_alpha_to_beta(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

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
      integer*8 indblk, get_block_index
      integer block, blkndx, seg
      integer find_current_block

      integer i, j, ii, source, target, company, comm, ierr
      integer source_type, target_type
      integer nindex_source, nindex_target, stack
      integer*8 isource, itarget 
      integer n(mx_array_index),ind(mx_array_index) 
      integer ind_in0_bgn(mx_array_index),ind_in1_bgn(mx_array_index)
      integer ind_in0_end(mx_array_index),ind_in1_end(mx_array_index)
      integer Idim,Jdim,Kdim,IJKdim
      integer Mdim,Ndim,Pdim,MNPdim
      integer Ibgn,Jbgn,Kbgn,Mbgn,Nbgn,Pbgn
      integer Iend,Jend,Kend,Mend,Nend,Pend
      integer msg(len_sip_server_message)
      integer status(MPI_STATUS_SIZE) 

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

#ifdef ALTIX
      dptr = dshptr
#endif
c--------------------------------------------------------------------------
c   Locate the data for both arrays.  
c--------------------------------------------------------------------------

      target = op(c_op1_array)
      source = op(c_result_array)

c---------------------------------------------------------------------------
c   Look up source's address.
c---------------------------------------------------------------------------

      block  = find_current_block(source, array_table(1,source),
     *                            index_table, nindex_table,
     *                            segment_table, nsegment_table,
     *                            block_map_table, blkndx)
      stack = array_table(c_array_stack,source)
      isource = get_block_index(source, block, stack, blkndx, x, .true.) 

c---------------------------------------------------------------------------
c   Look up target's address.
c---------------------------------------------------------------------------

      block  = find_current_block(target, array_table(1,target),
     *                            index_table, nindex_table,
     *                            segment_table, nsegment_table,
     *                            block_map_table, blkndx)
      stack = array_table(c_array_stack,target)
      itarget = get_block_index(target, block, stack, blkndx, x, .true.) 

c--------------------------------------------------------------------------
c   Check the dimensions of both arrays both arrays.  
c--------------------------------------------------------------------------
      
      source_type = array_table(c_array_type, source)
      target_type = array_table(c_array_type, target)

C      if (source_type .ne. distributed_array .or. target_type
C     &                   .ne. distributed_array) then
      if (source_type .ne. static_array  .or. target_type
     &                   .ne. static_array) then
         Write(6,"(a,a)") "@-COPY-ALPHA_TO_BETA: Both arrays must be",
     &                    " distributed arrays."
         call abort_job()
      endif

      nindex_source = array_table(c_nindex, source)
      nindex_target = array_table(c_nindex, target)

#ifdef _DEBUG_LVL0
      if (me.eq.0) then
      Write(6,*)
      Write(6,"(2a,2i4)") "@-Entry:Copy_alpha_to_beta The array 1 and 2",
     +                    " dimensions: ",
     +                    nindex_source,  nindex_target 
      endif 
#endif

      if (nindex_source .NE. nindex_target ) then
          Write(6,"(a,a,a)") "@-COPY_ALPHA_TO_BETA: The number of",
     &                     " indices of the two arrays must be",
     &                     " indentical and they are not!."  
          call abort_job()
      endif

c-------------------------------------------------------------------------
c   Get segment ranges of the array source.
c-------------------------------------------------------------------------

      do i = 1, nindex_source 
         ind(i) = array_table(c_index_array1+i-1,source)
         n(i)   = index_table(c_index_size, ind(i))  ! pick up length of index 
         seg    = index_table(c_current_seg,ind(i)) 

         call get_index_segment(ind(i), seg, segment_table,
     *                          nsegment_table, index_table,
     *                          nindex_table, ind_in0_bgn(i), 
     *                          ind_in0_end(i))
      enddo

c-------------------------------------------------------------------------
c   Get segment ranges of the array target.
c-------------------------------------------------------------------------

      do i = 1, nindex_target 
         ind(i) = array_table(c_index_array1+i-1,target)
         n(i)   = index_table(c_index_size, ind(i))  ! pick up length of index 
         seg    = index_table(c_current_seg,ind(i)) 

         call get_index_segment(ind(i), seg, segment_table,
     *                          nsegment_table, index_table,
     *                          nindex_table, ind_in1_bgn(i),
     *                          ind_in1_end(i))
      enddo

      If (nindex_source .EQ. 3 .AND. nindex_target .EQ. 3) Then

          Ibgn = ind_in0_bgn(1)
          Jbgn = ind_in0_bgn(2)
          Kbgn = ind_in0_bgn(3)
          Iend = ind_in0_end(1)
          Jend = ind_in0_end(2)
          Kend = ind_in0_end(3)

          Idim = Iend-Ibgn+1
          Jdim = Jend-Jbgn+1
          Kdim = Kend-Kbgn+1
          IJKdim = Idim * Jdim * Kdim

          Mbgn = ind_in1_bgn(1)
          Nbgn = ind_in1_bgn(2)
          Pbgn = ind_in1_bgn(3)
          Mend = ind_in1_end(1)
          Nend = ind_in1_end(2)
          Pend = ind_in1_end(3)

          Mdim = Mend-Mbgn+1
          Ndim = Nend-Nbgn+1
          Pdim = Pend-Pbgn+1
          MNPdim = Mdim * Ndim * Pdim

#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a,a,4i4)") "@-Entry to Copy_alpha_to_beta:",
     +                     " array2 variabels: ",
     +                      Mdim,Ndim,Pdim,Qdim
      Endif
#endif


      Call Copy_A_to_B_33(x(isource),Idim,Jdim,Kdim,
     +                    Mdim,Ndim,Pdim,Ibgn,Jbgn,Kbgn,
     +                    Iend,Jend,Kend,Mbgn,Nbgn,Pbgn,
     +                    Mend,Nend,Pend,x(itarget))

      Endif

      Return
      End

      Subroutine Copy_A_to_B_33(A,Idim,Jdim,Kdim,Mdim,Ndim,Pdim,
     +                          Ibgn,Jbgn,Kbgn,Iend,Jend,Kend,
     +                          Mbgn,Nbgn,Pbgn,Mend,Nend,Pend,B)

      implicit none
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      Integer I,J,K,Idim,Jdim,Kdim,Mdim,Ndim,Pdim
      Integer Ibgn,Jbgn,Kbgn,Iend,Jend,Kend
      Integer Mbgn,Nbgn,Pbgn,Mend,Nend,Pend

      Double Precision A(Ibgn:Iend,Jbgn:Jend,Kbgn:Kend)
      Double Precision B(Mbgn:Mend,Nbgn:Nend,Pbgn:Pend)

#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,*)
      Write(6,"(a,6i4)") "@-Copy_A_to_B_33 variables: ",
     +                    Idim,jdim,kdim,Mdim,Ndim,Pdim
      Endif
#endif

      DO I =Ibgn, Iend
      DO J =Jbgn, Jend

         DO K =Kbgn, Kend

#ifdef _DEBUG_LVL0
        if (me.eq.0) write(6,"(5I3,2F15.6)") i,j,k,
     +                                       A(I,J,K,L),
#endif
            B(I,J,K) = A(I,J,K)

         ENDDO

      ENDDO
      ENDDO

      Return
      End



