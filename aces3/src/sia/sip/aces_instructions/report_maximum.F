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
      subroutine Report_maximum(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)

      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'proto_events.h'
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
      integer*8 address_table(narray_table),get_index_from_base
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      double precision Max_2,Max_4
      double precision x(1), input_constant 

      integer eps,i,n,index4_max
      integer na1,na2,ni1,ni2,ierr
      integer array, index, nindex,arg
      integer array_type, inda(2)
      integer sind 
      integer atype,btype
      integer*8 addr
      integer*8 indblk, get_block_index
      integer block, blkndx, seg,stack
      integer master, pst_get_master
      integer status(mpi_status_size)
      integer find_current_block
      integer ind(mx_array_index)
      integer val1(mx_array_index), val2(mx_array_index)
      integer type(mx_array_index)
      integer company_comm, pst_get_company_comm
      integer ibuf(1)
      logical msg_present

      include 'int_gen_parms.h'
      include 'epsilon.h'

#ifdef ALTIX
      pointer (iptr, ibuf)
      pointer (dptr, x)
#else
      common x
      equivalence (ibuf, x)
#endif

c-------------------------------------------------------------------------
c   Locate array address in array_table.
c-------------------------------------------------------------------------

#ifdef ALTIX
      iptr = ishptr
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

c-------------------------------------------------------------------------
c   Get segment ranges of the array target.
c-------------------------------------------------------------------------

      array      = op(c_op1_array)
      array_type = array_table(c_array_type, array)

      if (array_type .ne. scalar_value) then
         write(6,"(a)") " Error: The second argument in",
     *                  " REPORT_MAXIMUM must be a scalar"
         call abort_job()
      endif

      sind = array_table(c_scalar_index, array)

      If (Nindex .Eq. 4) Then
         Call Find_maximum4(x(indblk),val1(1),val2(1),type(1),
     *                      val1(2),val2(2),type(2),
     *                      val1(3),val2(3),type(3),
     *                      val1(4),val2(4),type(4),Max_4)
C#ifdef _DEBUG_LVL0
        if (me .eq.0) then
        write(6,"(a,E8.2E2)") " The largest AAAA,BBBB or ABAB: ",
     *                          Max_4
        endif 
C#endif 
          Scalar_table(sind) = Max_4
           
      Elseif (Nindex .Eq. 2) Then
         na1 = val1(1)
         na2 = val2(1)
         ni1 = val1(2)
         ni2 = val2(2)
         Call Find_maximum2(x(indblk),val1(1),val2(1),type(1),
     *                      na1,na2,val1(2),val2(2),type(2),ni1,
     *                      ni2,Max_2)

C#ifdef _DEBUG_LVL0
        if (me .eq.0) then
        write(6,"(a,E8.2E2)") " The largest AA or BB         : ",
     *                          Max_2
        endif 
C#endif 
          Scalar_table(sind) = Max_2
      Endif 

      return
      end
C-------------------------------------------------------------------------

      Subroutine find_maximum4(xxx,a1,a2,atype,i1,i2,itype,b1,b2,
     *                         btype,j1,j2,jtype,Max_4)
      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'mpif.h'
      include 'proto_events.h'
      include 'parallel_info.h'

      integer a1,a2,b1,b2,i1,i2,j1,j2,i,j,k,l
      integer Itype, Jtype, Btype, Atype
      integer Lda,Ldb,Len
      double precision xxx(a1:a2,i1:i2,b1:b2,j1:j2)
      double precision zq(6),Max_4

      lda = (a2-a1+1)*(i2-i1+1)
      ldb = (b2-b1+1)*(j2-j1+1) 
      Len =  Lda*Ldb
C#ifdef _DEBUG_LVL0
        if (me .eq.0) then
        write(6,"(a,E8.2E2,1x,4I4)") " In find_maximum4, array length ", 
     +                         Dble(Len),a2,a1,j2,j1
        endif
        if (me .eq. 0) then
        Call Vstat_extension(xxx,zq,Len)
        endif 
C#ifdef _DEBUG_LVL0
        if (me .eq.0) then
        write(6,"(a,a,E8.2E2)") " In find_maximum4, The largest", 
     +                          " AAAA,BBBB or ABAB: ",Zq(1)
        endif
C#endif
      Max_4 = Zq(1)

      Return
      End

      Subroutine find_maximum2(xxx,a1,a2,atype,na1,na2,i1,i2,itype,
     *                         ni1,ni2,Max_2)

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'mpif.h'
      include 'proto_events.h'
      include 'parallel_info.h'

      integer a1,a2,i1,i2,Itype,Atype,i,j
      integer na1, na2, ni1, ni2,Len,lda,ldb
      double precision xxx(na1:na2,ni1:ni2)
      double precision zq(6),Max_2

      lda = (na2-na1+1)
      ldb = (ni2-ni1+1)
      Len =  Lda*Ldb
C#ifdef _DEBUG_LVL0
        if (me .eq.0) then
        write(6,"(a,E8.2E2)") " In find_maximum2, array length ", 
     +        Dble(Len)
        endif
C#endif

      Call Vstat_extension(xxx,zq,Len)
C#ifdef _DEBUG_LVL0
        if (me .eq.0) then
        write(6,"(a,a,E8.2E2)") " In find_maximum2, The largest",
     +                          " AA or BB: ",Zq(1)
        endif
C#endif
      Max_2 = Zq(1)

      Return
      End 

      SUBROUTINE VSTAT_EXTENSION(V,ZQ,LENGTH)
C
C RETURNS STATISTICAL INFO ABOUT VECTOR V in ZQ
C     ZQ(1)  Largest absolute magnitude
C     ZQ(2)  Smallest absolute magnitude
C     ZQ(3)  Largest value
C     ZQ(4)  Smallest value
C     ZQ(5)  2-norm
C     ZQ(6)  Dynamic range of the vector (abs. min. - abs. max.)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(LENGTH),ZQ(6)
      U=0.D0
      CALL ZERO(ZQ,6)
      ZQ(2)=DABS(V(1))
      ZQ(4)=V(1)
      DO 20 I=1,LENGTH
      ZQ(1)=MAX(ZQ(1),DABS(V(I)))
      ZQ(2)=MIN(ZQ(2),DABS(V(I)))
      ZQ(3)=MAX(ZQ(3),V(I))
      ZQ(4)=MIN(ZQ(4),V(I))
20    U=U+V(I)*V(I)
      If (Length .ne. 0) ZQ(5)=DSQRT(U/LENGTH)
   
      RETURN
      END
