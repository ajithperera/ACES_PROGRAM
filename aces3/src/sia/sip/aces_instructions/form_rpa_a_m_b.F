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
      subroutine Form_hphp_matrix(array_table, narray_table, 
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
      integer Idim,Jdim,Kdim,Ldim,IJdim,KLdim
      integer Mdim,Ndim,Pdim,Qdim,MNdim,PQdim
      integer Ibgn,Jbgn,Kbgn,Lbgn,Mbgn,Nbgn,Pbgn,Qbgn
      integer Iend,Jend,Kend,Lend,Mend,Nend,Pend,Qend
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

      if (source_type .ne. temp_array .or. target_type
     &                   .ne. temp_array) then
         Write(6,"(a,a)") "@-FORM_HPHP_MATRIX: Both arrays must be",
     &                    " temp. arrays."
         call abort_job()
      endif

      nindex_source = array_table(c_nindex, source)
      nindex_target = array_table(c_nindex, target)
#ifdef _DEBUG_LVL0
      if (me.eq.0) then
      Write(6,*)
      Write(6,"(a,2i4)") "@-Entry: The array 1 and 2 dimensions: ",
     +                   nindex_source,  nindex_target 
      endif 
#endif

      if (nindex_source .ne. 4 .AND. (nindex_target .NE. 2 .or.
     &    nindex_target .NE. 4)) then
          Write(6,"(a,a)") "@-FORM_HPHP_MATRIX: The first array",
     &                     " must always be a 4 index"
          Write(6,"(a,a)") " and the second array must be 4 or 2 index",
     &                     " array."  
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

#ifdef _DEBUG_LVL0
      if (me.eq.0) then
      Write(6,*)
      Write(6,"(a,2i4)") "@-Entry: The array 1 and 2 dimensions: ",
     +                   nindex_source,  nindex_target 
      endif
#endif 

      If (nindex_source .EQ. 2) Then
          Ibgn = ind_in0_bgn(1)
          Jbgn = ind_in0_bgn(2)
          Iend = ind_in0_end(1)
          Jend = ind_in0_end(2)
          Idim = Iend-Ibgn+1
          Jdim = Iend-Jbgn+1

#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a,a,2i4)") "@-Entry to Form_hphp_matirx:",
     +                     " array1 variabels: ", Idim,Jdim
      Endif
#endif

      Else if (nindex_source .EQ. 4) Then
          Ibgn = ind_in0_bgn(1)
          Jbgn = ind_in0_bgn(2)
          Kbgn = ind_in0_bgn(3)
          Lbgn = ind_in0_bgn(4)
          Iend = ind_in0_end(1)
          Jend = ind_in0_end(2)
          Kend = ind_in0_end(3)
          Lend = ind_in0_end(4)
          Idim = Iend-Ibgn+1
          Jdim = Jend-Jbgn+1
          Kdim = Kend-Kbgn+1
          Ldim = Lend-Lbgn+1
          IJdim= Idim * Jdim 
          KLdim= Kdim * Ldim 

#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a,a,6i4)") "@-Entry to Form_hphp_matirx:",
     +                   " array1 variabels: ", Idim,Jdim,Kdim,Ldim,  
     +                     Ijdim,KLdim
      Endif
#endif 
      Endif 

      If (nindex_target .EQ. 2) Then
          Mbgn = ind_in1_bgn(1)
          Nbgn = ind_in1_bgn(2)
          Mend = ind_in1_end(1)
          Nend = ind_in1_end(2)
          Mdim = Mend-Mbgn+1
          Ndim = Nend-Nbgn+1
          MNdim= Mdim * Ndim 

#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a,a,2i4)") "@-Entry to Form_hphp_matirx:",
     +                    " array2 variabels: ",Mdim,Ndim
      Endif
#endif
      Else if (nindex_target .EQ. 4) Then
          Mbgn = ind_in1_bgn(1)
          Nbgn = ind_in1_bgn(2)
          Pbgn = ind_in1_bgn(3)
          Qbgn = ind_in1_bgn(4)
          Mend = ind_in1_end(1)
          Nend = ind_in1_end(2)
          Pend = ind_in1_end(3)
          Qend = ind_in1_end(4)
          Mdim = Mend-Mbgn+1
          Ndim = Nend-Nbgn+1
          Pdim = Pend-Pbgn+1
          Qdim = Qend-Qbgn+1

          MNdim= Mdim * Ndim 
          PQdim= Pdim * Qdim 

#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a,a,4i4)") "@-Entry to Form_hphp_matirx:",
     +                     " array2 variabels: ",
     +                      Mdim,Ndim,Pdim,Qdim
      Endif
#endif
      Endif
#ifdef _ASSERT 
      If (nindex_source .EQ. 4) Then
          If (me.eq.0) write(6,*)  "IDim and KLdim", IJdim, KLdim 
         If (IJdim .Ne. KLdim) Then
             If (me.eq.0) write(6,*)  "IDim and KLdim", IJdim, KLdim 
             Write(6,"(a,a)") "The raw and column dimensions of the",
     +                        " first matrix do not match."
             Call abort_job()
         Endif
      Endif 

      If (nindex_source .EQ. 4 .and. nindex_target .EQ.2) Then
          If (IJdim .Ne. MNdim) Then
             Write(6,"(a,a)") "The dimensions of source and target",
     +                        " matrices do not match."
             Call abort_job()
          Endif
      Elseif (nindex_target .EQ. 4) Then
          If (MNdim .Ne. PQdim) Then
             Write(6,"(a,a)") "The raw and column dimensions of the",
     +                        " second matrix do not match."
             Call abort_job()
          Endif
          If (IJdim .Ne. MNdim .AND. IJdim .Ne. PQdim) Then
             Write(6,"(a,a)") "The raw and column dimensions of the",
     +                        " second matrix do not match."
              Call abort_job()
          Endif
      Endif
#endif 

C At the moment I consider two possible options. Other two possibilities
C can be added as needed.

      If (Nindex_source .EQ. 4 .AND. Nindex_target .EQ. 2) Then

          Call Form_Amat_42(x(isource),Idim,Jdim,Kdim,Ldim,
     +                      Mdim,Ndim,Ibgn,Jbgn,Kbgn,Lbgn,
     +                      Iend,Jend,Kend,Lend,Mbgn,Nbgn,
     +                      Mend,Nend,x(itarget))
      Endif

      If (Nindex_source .EQ. 4 .AND. Nindex_target .EQ. 4) Then

          Call Form_Amat_44(x(isource),Idim,Jdim,Kdim,Ldim,
     +                      Mdim,Ndim,Pdim,Ibgn,Jbgn,Kbgn,Lbgn,
     +                      Iend,Jend,Kend,Lend,Mbgn,Nbgn,Pbgn,
     +                      Qbgn,Mend,Nend,Pend,Qend,x(itarget))
      Endif

      Return
      End

      Subroutine Form_Amat_42(A_hphp,Idim,Jdim,Kdim,Ldim,Mdim,Ndim,
     +                        Ibgn,Jbgn,Kbgn,Lbgn,Iend,Jend,Kend,
     +                        Lend,Mbgn,Nbgn,Mend,Nend,B_hp)

      implicit none
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      Integer I,J,K,L,Idim,Jdim,Kdim,Ldim,IA
      Integer Ibgn,Jbgn,Kbgn,Lbgn,Iend,Jend,Kend,Lend
      Integer Mbgn,Nbgn,Mend,Nend
      Integer M,N,Mdim,Ndim,MNdim

      Double Precision A_hphp(Ibgn:Iend,Jbgn:Jend,Kbgn:Kend,Lbgn:Lend)
      Double Precision B_hp(Mbgn:Mend,Nbgn:Nend)

#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,*)
      Write(6,"(a,6i4)") "@-Form_Amat_42 @-Entry: variables: ",
     +                    Idim,jdim,kdim,Ldim,Mdim,Ndim
      Endif
      Write(6,"(a))") "(eps_i - eps_a)"
      call output(B_hp,Mbgn,Mend,Nbgn,Nend,Mdim,Ndim,1)
#endif
      DO L =Lbgn, Lend
      DO K =Kbgn, Kend
         DO J =Jbgn, Jend
         DO I =Ibgn, Iend
            IF (L .EQ. J .AND. K .EQ. I )  Then
#ifdef _DEBUG_LVL0
        if (me.eq.0) write(6,"(5I3,2F15.6)") i,j,k,l,ia,
     +                                       A_hphp(I,J,K,L),
     +                                       B_hp(L,K)
#endif
            A_hphp(I,J,K,L) = A_hphp(I,J,K,L) + B_hp(L,K)

            Endif
         ENDDO
         ENDDO
      eNDDO
      ENDDO

      Return
      End

      Subroutine Form_Amat_44(A_hphp,Idim,Jdim,Kdim,Ldim,Mdim,Ndim,
     +                        Pdim,Qdim,Ibgn,Jbgn,Kbgn,Lbgn,Iend,
     +                        Jend,Kend,Lend,Mbgn,Nbgn,Pbgn,Qbgn,
     +                        Mend,Nend,Pend,Qend,B_hphp)

      implicit none
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      Integer I,J,K,L,Idim,Jdim,Kdim,Ldim
      Integer Ibgn,Jbgn,Kbgn,Lbgn,Iend,Jend,Kend,Lend
      Integer Mbgn,Nbgn,Pbgn,Qbgn,Mend,Nend,Pend,Qend
      Integer M,N,P,Q,Mdim,Ndim,Pdim,Qdim

      Double Precision A_hphp(Ibgn:Iend,Jbgn:Jend,
     +                        Kbgn:Kend,Lbgn:Lend)
      Double Precision B_hphp(Mbgn:Mend,Nbgn:Nend,
     +                        Pbgn:Pend,Qbgn:Qend)

#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a,3i4)") "@-Form_Amat_44 @-Entry: variables: ",
     +                    Idim,jdim,kdim,Ldim
      Endif
#endif

      DO L =Lbgn, Lend
      DO K =Kbgn, Kend
         DO J =Jbgn, Jend
         DO I =Ibgn, Iend
            A_hphp(I,J,K,L) = A_hphp(I,J,K,L) + B_hphp(I,J,K,L)
         ENDDO
         ENDDO
      ENDDO
      ENDDO

      Return
      End


