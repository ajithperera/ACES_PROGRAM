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
      subroutine Drive_core_states(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'mpif.h'
      include 'parallel_info.h'
      include 'int_gen_parms.h'
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

      integer i, j, k
      integer array, index, nindex, ierr, sind, arg
      integer block, blkndx, seg, array_type 
      integer find_current_block
      integer*8 indblk, get_block_index
      integer stack
      integer comm,Iprj 
      Logical Sta

      integer val1(mx_array_index), val2(mx_array_index)
      integer type(mx_array_index)
      integer na1, na2, ni1, ni2,nk1,nk2
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

c---------------------------------------------------------------------------
c   Find the indices of the array block.
c---------------------------------------------------------------------------
       
      array = op(c_result_array)
      nindex = array_table(c_nindex, array)
      do i = 1, nindex
         index = array_table(c_index_array1+i-1,array)
         type(i) = index_table(c_index_type, index)
         seg = index_table(c_current_seg,index)

c-------------------------------------------------------------------------
c   Get segment ranges.
c-------------------------------------------------------------------------

         call get_index_segment(index, seg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1(i), val2(i))
      enddo

c---------------------------------------------------------------------------
c   Get array data address.
c---------------------------------------------------------------------------

      if (array_table(c_array_type,array) .eq. static_array) then
         addr = address_table(array)
         indblk = get_index_from_base(addr, x, 2)
      else
         block = find_current_block(array, array_table(1,array),
     *                             index_table, nindex_table,
     *                             segment_table, nsegment_table,
     *                             block_map_table, blkndx)

         stack = array_table(c_array_stack,array)
         indblk = get_block_index(array, block, stack,
     *                            blkndx, x, .true.)
      endif 

C The second argument indicate which side get scalled.

      array      = op(c_op1_array)
      array_type = array_table(c_array_type, array)

      if (array_type .ne. scalar_value) then
         write(6,"(a)") " Error: The second argument in",
     *                  " DRIVE_CORE_STATES must be a scalar"
         call abort_job()
      endif

      sind = array_table(c_scalar_index, array)
      Iprj = Scalar_table(sind)

C#ifdef _DEBUG_LVL0
      If (me.eq.0) write(6,"(2a,i2,1x,i2)") " Incoming value of", 
     &                                " Iprj and Nindex: ", Nindex,
     &                                  Iprj
C#endif 

      If (nindex .eq. 4) then
         If (Iprj .Eq. 1) Then
         call Project_occ41(x(indblk), val1(1),val2(1),type(1),
     *                 val1(2),val2(2),type(2),
     *                 val1(3),val2(3),type(3),
     *                 val1(4),val2(4),type(4))
         Else if (Iprj .Eq. 2) Then 
         call Project_occ42(x(indblk), val1(1),val2(1),type(1),
     *                 val1(2),val2(2),type(2),
     *                 val1(3),val2(3),type(3),
     *                 val1(4),val2(4),type(4))
         Endif 
      Else if (nindex .eq. 2) then
         array = op(c_result_array)
         if (array_table(c_array_type,array) .eq. static_array) then
            na1 = array_table(c_index_range1, array)
            na2 = array_table(c_index_range2,array)
            ni1 = array_table(c_index_range1+1,array)
            ni2 = array_table(c_index_range2+1,array)
            sta = .True.
         else
            na1 = val1(1)
            na2 = val2(1)
            ni1 = val1(2)
            ni2 = val2(2)
         endif 

         If (Iprj .Eq. 1) Then
            call Project_occ21(x(indblk), val1(1),val2(1),type(1),
     *                 na1,na2,val1(2),val2(2),type(2),ni1,ni2)
         Else 

            call Set_guess_2(x(indblk),val1(1),val2(1),type(1),
     *                     na1,na2,val1(2),val2(2),type(2),ni1,ni2,
     *                     Sta,Iprj)
         Endif 

      Else if (nindex .eq. 3) then
         array = op(c_result_array)
         if (array_table(c_array_type,array) .eq. static_array) then
            na1 = array_table(c_index_range1, array)
            na2 = array_table(c_index_range2,array)
            ni1 = array_table(c_index_range1+1,array)
            ni2 = array_table(c_index_range2+1,array)
            nk1 = array_table(c_index_range1+1,array)
            nk2 = array_table(c_index_range2+1,array)
            sta = .True.
         else
            na1 = val1(1)
            na2 = val2(1)
            ni1 = val1(2)
            ni2 = val2(2)
            nk1 = val1(3)
            nk2 = val2(3)
         endif

         If (Iprj .Eq. -1) Then
            call Set_guess_3(x(indblk),val1(1),val2(1),type(1),
     *                     na1,na2,val1(2),val2(2),type(2),ni1,ni2,
     *                     Nk1,Nk2,Sta)
         Endif 

      Else
         write(6,"(a,a)") " Error in DRIVE_CORE_STATES: target array",
     *                    " can not have more than 4 indices" 
         Write(6,"(a)")   " or the second argument can not be zero."
         Write(6,"(a,i2,a,i2,a)") " Currently it has", nindex,
     *                       " indices and ", Iprj, "second argument."
         call abort_job()

      endif

      return
      end

      subroutine Project_occ42(x,a1,a2,atype,i1,i2,itype,b1,b2,
     *                         btype,j1,j2,jtype) 

      implicit none
      include 'interpreter.h'
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'ex_mask_vars.h'

      integer a1,a2,b1,b2,i1,i2,j1,j2
      integer Itype, Jtype, Btype, Atype 
      double precision x(a1:a2,i1:i2,b1:b2,j1:j2)

      integer a,b,i,j

      do j = j1,j2
      do b = b1,b2
      do i = i1,i2
      do a = a1,a2
         If (Itype .EQ. Moaindex .and. Jtype .EQ. Moaindex) Then 
            If (Ibelong_a(I) .NE. 0 .OR. Ibelong_a(J) .NE. 0) Then
               x(a,i,b,j) = x(a,i,b,j)
            Else 
               x(a,i,b,j) = 0.0D0
            Endif 
         Elseif (Itype .EQ. Moaindex .and. Jtype .EQ. Mobindex) Then
            If (Ibelong_a(I) .NE. 0 .OR. Ibelong_b(J) .NE. 0) Then
              x(a,i,b,j) =  x(a,i,b,j)
            Else 
               x(a,i,b,j) = 0.0D0 
            Endif 
         Elseif (Itype .EQ. Mobindex .and. Jtype .EQ. Moaindex) Then 
            If (Ibelong_b(I) .NE. 0 .OR. Ibelong_a(J) .NE. 0) Then
              x(a,i,b,j) =  x(a,i,b,j)
            Else 
              x(a,i,b,j) =  0.0D0
            Endif 
         Elseif (Itype .EQ. Mobindex .and. Jtype .EQ. Mobindex) Then 
            If (Ibelong_b(I) .NE. 0 .OR. Ibelong_b(J) .NE. 0) Then
              x(a,i,b,j) =  x(a,i,b,j)
            Else 
              x(a,i,b,j) =  0.0D0
            Endif 
         Endif 
      enddo
      enddo
      enddo
      enddo

      return
      end

      subroutine Project_occ41(x,a1,a2,atype,i1,i2,itype,b1,b2,
     *                         btype,j1,j2,jtype)

      implicit none
      include 'interpreter.h'
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'ex_mask_vars.h'
      include 'aces3_vars.h'

      integer a1,a2,b1,b2,i1,i2,j1,j2
      integer Itype, Jtype, Btype, Atype
      double precision x(a1:a2,i1:i2,b1:b2,j1:j2)

      integer a,b,i,j

      If (Ip) Then 
      do j = j1,j2
      do b = b1,b2
      do i = i1,i2
      do a = a1,a2
         If (Jtype .EQ. Moaindex) Then
            If (Ibelong_a(J) .NE. 0) Then
               x(a,i,b,j) = x(a,i,b,j)
            Else
               x(a,i,b,j) = 0.0D0
            Endif
         Elseif (Jtype .EQ. Mobindex) Then
            If (Ibelong_b(J) .NE. 0) Then
              x(a,i,b,j) =  x(a,i,b,j)
            Else
               x(a,i,b,j) = 0.0D0
            Endif
         Elseif (Jtype .EQ. Moaindex) Then
            If (Ibelong_a(J) .NE. 0) Then
              x(a,i,b,j) =  x(a,i,b,j)
            Else
              x(a,i,b,j) =  0.0D0
            Endif
         Elseif (Jtype .EQ. Mobindex) Then
            If (Ibelong_b(J) .NE. 0) Then
              x(a,i,b,j) =  x(a,i,b,j)
            Else
              x(a,i,b,j) =  0.0D0
            Endif
         Endif
      enddo
      enddo
      enddo
      enddo
      Endif 

      Return
      End 

      subroutine Project_occ21(x,a1,a2,atype,na1,na2,i1,i2,itype,
     *                         ni1,ni2)
      implicit none
      include 'interpreter.h'
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'ex_mask_vars.h'
      include 'aces3_vars.h'

      integer a1,a2,i1,i2,Itype,Atype
      integer na1, na2, ni1, ni2
      double precision x(na1:na2,ni1:ni2)

      integer a,i

      If (Ee) Then

      do i = i1, i2
      do a = a1, a2
         If (Itype .EQ. Moaindex) Then
            If (Ibelong_a(I) .NE. 0) Then
               x(a,i) = x(a,i)
            Else
               x(a,i) = 0.0D0
            Endif 
         Else
            If (Ibelong_b(I) .NE. 0) Then
               x(a,i) = x(a,i)
            Else
               x(a,i) = 0.0D0
            Endif 
         Endif 
      enddo
      enddo

      Elseif (Ip) then

      do a = a1, a2
      do i = i1, i2
         If (Itype .EQ. Moaindex) Then
            If (Ibelong_a(I) .NE. 0) Then
               x(i,a) = x(i,a)
            Else
               x(i,a) = 0.0D0
            Endif
         Else
            If (Ibelong_b(I) .NE. 0) Then
               x(i,a) = x(i,a)
            Else
               x(i,a) = 0.0D0
            Endif
         Endif
      enddo
      enddo

      endif 


      return
      end

      subroutine Set_guess_3(x,a1,a2,atype,na1,na2,i1,i2,itype,
     *                       ni1,ni2,Nk1,Nk2,Sta)
      implicit none
      include 'interpreter.h'
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'ex_mask_vars.h'

      integer a1,a2,i1,i2,Itype,Atype
      logical Lda,Ldb,Ldc
      integer na1,na2,ni1,ni2,nk1,nk2
      integer Iroots,Nroots,Isites,Ioff 
      double precision x(na1:na2,ni1:ni2,nk1:nk2)
      Logical Sta

      integer a,i

      If (.Not. Sta) Then
         Write(6,"(2a)") " @-DRIVE_CORE_STATES: The incomming array is",
     +                   " not static"
         call abort_job()
      Endif 

      Lda = Na2-Na1+1
      Ldb = Ni2-Ni1+1
      Ldc = Nk2-Nk1+1

      Nroots = Eom_roots 
      Nroots = Nroots * Nsites 

      If (Nroots .GT. 12 .OR. Ldc .GT. 12) Then
         Write(6,"(2a)") " The requested number of roots exceeded",
     +                   " the maximum allowed 12 roots." 
         Write(6,"(a,i3,a,i3)") " The number of sites is ", Nsites,
     +                   " and the number of roots requested per",
     +                   "  site is ", Eom_roots
         Write(6,"(2a)") " Change these limits in the sial codes"
     +                   " resubmit"
         Call abort_job()
      Endif 
C#ifdef _DEBUG_LVL0
      If (me .eq. 0) Then 
      Write(6,*) 
      Write(6,"(a,2i3)") " The index ranges-row: ", na1,na2
      Write(6,"(a,2i3)") " The index ranges-col: ", ni1,ni2
      Write(6,"(a)") " The markers for the excitation mask"
      Write(6,"(6I4)") (Ibelong_a(i),i=1, Nalpha_occupied)
      Write(6,"(6I4)") (Ibelong_b(i),i=1, Nbeta_occupied)
      endif 
C#endif 

      Do Isites = 1, Nsites 
         Ioff = 0
      Do Iroots = 1, Eom_roots 
         ioff =  Ioff + 1 
      Do i = ni1, ni2
      Do a = na1, na2
         If (Itype .EQ. Moaindex) Then
            If (Ibelong_a(I) .NE. 0 .AND. A .EQ. Nalpha_occupied+Ioff)
     +      Then
               x(a,i,Iroots) = 1.0D0
            Else
               x(a,i,Iroots) = 0.0D0
            Endif
         Else if (Itype .EQ. Mobindex) Then 
            If (Ibelong_b(I) .NE. 0 .AND. A .EQ. Nalpha_occupied+Ioff)
     +      Then
               x(a,i,Iroots) = 1.0D0
            Else
               x(a,i,Iroots) = 0.0D0
            Endif
         Endif
      Enddo
      Enddo
      Enddo
      Enddo 

C#ifdef _DEBUG_LVL0
      If (me.eq.0) then 
      do iroots=1,Eom_roots
      Write(6,"(a,a,I2)") " The initial guess @-DRIVE_CORE_STATES for",
     +                    " root: ", iroots
      Write(6,"(8(1x,ES8.2E2))") ((x(a,i,iroots),a=na1,na2),
     +                                           i=ni1,ni2)
      enddo 
      endif 
C#endif 

      Return
      End 
    
      subroutine Set_guess_2(x,a1,a2,atype,na1,na2,i1,i2,itype,
     *                       ni1,ni2,Sta,Iwant)
      implicit none
      include 'interpreter.h'
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'ex_mask_vars.h'

      integer a1,a2,i1,i2,Itype,Atype
      logical Lda,Ldb,Ldc,Iwant
      integer na1,na2,ni1,ni2
      integer Iroots,Nroots,Isites,Ipick 
      double precision x(na1:na2,ni1:ni2)
      Integer Imap(100)
      Logical Sta,Allowed 

      integer a,i

      If (.Not. Sta) Then
         If (me .eq. 0) then
         Write(6,"(2a)") " @-DRIVE_CORE_STATES: The incomming array is",
     +                   " not static"
         call abort_job()
         endif 
      Endif 

      Lda = Na2-Na1+1
      Ldb = Ni2-Ni1+1

      Nroots  = Eom_roots 
      Nroots  = Nroots * Nsites 
      Allowed = Nroots .Le. 12 .AND. (.Not. (Nroots .Gt. 100))

      If (.Not. Allowed) Then
         If (me .eq. 0) then
         Write(6,"(2a)") " The requested number of roots exceeded",
     +                   " the maximum allowed 12 roots." 
         Write(6,"(a,i2,2a,i2)") " The number of sites is ", Nsites,
     +                   " and the number of roots requested per",
     +                   " site is ", Nroots 
         Write(6,"(2a)") " Change these limits in the sial codes",
     +                   " and resubmit."
         Write(6,*) 
         Call abort_job()
         endif 
      Endif 

#ifdef _DEBUG_LVL0
      If (me .eq. 0) Then 
      write(6,*)
      Write(6,"(a,2i3)") " The index ranges-row: ", na1,na2
      Write(6,"(a,2i3)") " The index ranges-col: ", ni1,ni2
      Write(6,"(a)") " The markers for the excitation mask"
      Write(6,"(6I4)") (Ibelong_a(i),i=1, Nalpha_occupied)
      Write(6,"(6I4)") (Ibelong_b(i),i=1, Nbeta_occupied)
      endif 
#endif 

      Do I = 1, Nroots
         Imap(9+I) = I
      Enddo 

      Do Isites = 1, Nsites 
         Ipick = Imap(Iwant) 
      Do i = ni1, ni2
      Do a = na1, na2
         If (Itype .EQ. Moaindex) Then
            If (Ibelong_a(I) .NE. 0 .AND. A .EQ. Nalpha_occupied+Ipick)
     +      Then
               x(a,i) = 1.0D0
            Else
               x(a,i) = 0.0D0
            Endif
         Else if (Itype .EQ. Mobindex) Then 
            If (Ibelong_b(I) .NE. 0 .AND. A .EQ. Nalpha_occupied+Ipick)
     +      Then
               x(a,i) = 1.0D0
            Else
               x(a,i) = 0.0D0
            Endif
         Endif
      Enddo
      Enddo
      Enddo 

#ifdef _DEBUG_LVL0
      If (me.eq.0) then 
      Write(6,"(a,a,I2)") " The initial guess @-DRIVE_CORE_STATES",
     +                    " for root ", Imap(Iwant)
      Write(6,"(8(1x,ES8.2E2))") ((x(a,i),a=na1,na2),
     +                                    i=ni1,ni2)
      endif 
#endif 

      Return
      End 



