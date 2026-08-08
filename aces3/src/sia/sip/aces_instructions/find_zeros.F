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
      subroutine Find_zeros(array_table, narray_table, 
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
      integer*8 addr, get_index_from_base
      integer n(mx_array_index),ind(mx_array_index)
      integer ind_in(mx_array_index),ind_out(mx_array_index)
      integer ind_in0_bgn(mx_array_index),ind_in1_bgn(mx_array_index)
      integer ind_in0_end(mx_array_index),ind_in1_end(mx_array_index)
      integer Ibgn,Jbgn,Kbgn,Lbgn
      integer Iend,Jend,Kend,Lend
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

      source = op(c_result_array)
      target = op(c_op1_array)

      source_type = array_table(c_array_type, source)
      target_type = array_table(c_array_type, target)

      if (.NOT.(source_type .eq. static_array  .or. source_type 
     &                   .eq. static_array)) then
         Write(6,"(a,a)") "@-Strip_zero_vecs: Array must be",
     &                    " static arrays"
         call abort_job()
      endif

      nindex_source = array_table(c_nindex, source)
      nindex_target = array_table(c_nindex, target)

#ifdef _DEBUG_LVL0
      if (me.eq.0) then
      Write(6,*)
      Write(6,"(2a,2i4)") "@-Entry:FIND_ZEROS The array 1 and 2",
     +                    " dimensions: ",
     +                      nindex_source,  nindex_target 
      endif 
#endif
      Do i = 1, nindex_source
         ind_in(i) = array_table(c_index_array1+i-1,source)
      Enddo 
      do j = 1, nindex_target
         ind_out(j) = array_table(c_index_array1+j-1,target)
      enddo

      Ibgn = 1
      Jbgn = 1
      Kbgn = 1
      Lbgn = 1
      Iend = index_table(c_index_size, ind_in(1))
      Jend = index_table(c_index_size, ind_in(2))
      Kend = index_table(c_index_size, ind_out(1))
      Lend = index_table(c_index_size, ind_out(2))

      isource  = get_index_from_base(address_table(source), x, 2)
      itarget  = get_index_from_base(address_table(target), x, 2)

#ifdef _DEBUG_LVL0
      if (me.eq.0) Write(6,"(8i4)") Ibgn,Jbgn,Kbgn,Lbgn,Iend,
     +             Jend,Kend,Lend
#endif

      Call Flag_zeros(x(isource),Ibgn,Jbgn,Iend,Jend,x(itarget),
     +                Kbgn,Lbgn,Kend,Lend)

      Return
      End

      Subroutine Flag_zeros(A,Ibgn,Jbgn,Iend,Jend,B,Kbgn,Lbgn,
     +                      Kend,Lend)

      implicit none
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      Integer I,J,Ndim,Mdim
      Integer Ibgn,Jbgn,Iend,Jend,Kbgn,Lbgn,Kend,Lend

      Double Precision A(Ibgn:Iend,Jbgn:Jend)
      Double Precision B(Kbgn:Kend,Lbgn:Lend)
      Double Precision Sum,Threshold
      Integer Nroots

      Nroots     = Eom_roots 
      Threshold  = eom_tol 
      Mdim = (Iend-Ibgn+1) * (Jend-Jbgn+1)
      Ndim = (Kend-Kbgn+1) * (Lend-Lbgn+1)

      If (Mdim .NE. Ndim) Then
         If (me. eq. 0) Then 

          Write(6,"(8i4)") Ibgn,Jbgn,Kbgn,Lbgn,Iend,
     +             Jend,Kend,Lend
         Write(6,"(a,a)") "@-FLAG_ZEORS; Both source and target arrays",
     +                    " must be of the same dimensions"
         Call abort_job()
         Endif 
      Endif 

      Call Zero(B,Ndim)

      If (Nroots .Gt. 12) Then
         Write(6,"(3a)") "@-FLAG_zeros: The current maximum",
     +                  " number of roots that", 
     +                  " can be proccessed is 12."
         Call abort_job()
      Endif 

#ifdef _DEBUG_LVL0
      If (me .eq. 0) then
      Write(6,"(a,i5,F15.6)") "@-Entry Nroots, Threshold: ",
     +                         Nroots,Threshold
      Write(6,"(a)") "At entry FLAG_ZEROS"
      Write(6,"(5(1x,F15.6))") ((A(i,j),i=ibgn,iend), j=jbgn,jend)
      Endif 
#endif 

      DO J =Jbgn, Jend

         DO I =Ibgn, Iend
            Sum  = Dabs(A(I,J))
            If (Sum .Lt. Threshold) B(I,J) = 1.0
         ENDDO
      ENDDO 

#ifdef _DEBUG_LVL0
      If (me .eq. 0) then
      Write(6,"(a)") "At return array elements from FLAG_ZEORS"
      Write(6,"(5(1x,F15.6))") ((B(i,j),i=ibgn,iend), j=jbgn,jend)
      Endif 
#endif 

      Return
      End
