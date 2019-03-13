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
      subroutine apply_frac_occs(array_table, narray_table, 
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
      include 'epsilon.h'
      include 'frac_occ_vars.h'
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
      integer array, index, nindex, ierr, sind
      integer block, blkndx, seg, array_type 
      integer find_current_block
      integer*8 indblk, get_block_index
      integer stack
      
      integer comm

      integer val1(mx_array_index), val2(mx_array_index)
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
     *                  " apply_frac_occs must be a scalar"
         call abort_job()
      endif

      sind =  array_table(c_scalar_index, array)

      if (nindex .eq. 4 .and. contract .eq. 1) then
         call frac_occ41(x(indblk), val1(1),val2(1),type(1),
     *                 val1(2),val2(2),type(2),
     *                 val1(3),val2(3),type(3),
     *                 val1(4),val2(4),type(4), 
     *                 epsilon,epsilonb,Sind)
      elseif (nindex .eq. 4 .and. contract .eq. 2) then
        call frac_occ42(x(indblk), val1(1),val2(1),type(1),
     *                 val1(2),val2(2),type(2),
     *                 val1(3),val2(3),type(3),
     *                 val1(4),val2(4),type(4),
     *                 epsilon,epsilonb,Sind)
      elseif (nindex .eq. 4 .and. contract .eq. 3) then
        call frac_occ43(x(indblk), val1(1),val2(1),type(1),
     *                 val1(2),val2(2),type(2),
     *                 val1(3),val2(3),type(3),
     *                 val1(4),val2(4),type(4),
     *                 epsilon,epsilonb,Sind)
      elseif (nindex .eq. 4 .and. contract .eq. 4) then
        call frac_occ44(x(indblk), val1(1),val2(1),type(1),
     *                 val1(2),val2(2),type(2),
     *                 val1(3),val2(3),type(3),
     *                 val1(4),val2(4),type(4),
     *                 epsilon,epsilonb,Sind)
      else if (nindex .eq. 2) then
         if (array_table(c_array_type,array) .eq. static_array) then
            na1 = array_table(c_index_range1, array)
            na2 = array_table(c_index_range2,array)
            ni1 = array_table(c_index_range1+1,array)
            ni2 = array_table(c_index_range2+1,array)
         else
            na1 = val1(1)
            na2 = val2(1)
            ni1 = val1(2)
            ni2 = val2(2)
         endif 

         If (Contract .Eq. 1) Then
         call frac_occ21(x(indblk), val1(1),val2(1),type(1),
     *                 na1, na2,
     *                 val1(2),val2(2),type(2),ni1,ni2,
     *                 epsilon,epsilonb,sind)    
         Elseif (Contract .Eq. 2) then
         call frac_occ22(x(indblk), val1(1),val2(1),type(1),
     *                 na1, na2,
     *                 val1(2),val2(2),type(2),ni1,ni2,
     *                 epsilon,epsilonb,sind)    
         Endif 
      
      else
        
         write(6,"(a,a)") " Error in apply_frac_occs: target array",
     *                    " must have either 4 or 2 indices." 
         Write(6,"(a,i2,a)") " array has", nindex, "indices"
         call abort_job()
      endif

      return
      end

      subroutine frac_occ41(x, a1,a2,atype,
     *                     i1,i2,itype,b1,b2,btype,
     *                     j1,j2, jtype, epsilon, epsilonb,sind)
      implicit none
      include 'interpreter.h'
      include 'frac_occ_vars.h'

      integer a1,a2,b1,b2,i1,i2,j1,j2
      integer atype, itype, btype, jtype
      double precision x(a1:a2,i1:i2,b1:b2,j1:j2)

      double precision epsilon(*), epsilonb(*)

      integer a,b,i,j,sind
      double precision eps, epsa, epsb, epsi, epsj
      double precision val,Done 
      Logical Lft_one,Rit_one
      Logical O_type,V_type

      Data Done /1.0D0/
  
      If (sind .eq. 1) Then
         Lft_one = .True.
      Elseif (sind .eq. 2) then
         Rit_one = .True.
      Endif 

      If (O .eq. 1) Then 
         O_type = .True.
      ElseIf (V .eq. 1) Then 
         V_type = .True.
      Endif

      do j = j1,j2
      do b = b1,b2
      do i = i1,i2
      do a = a1,a2
         If (Lft_one) Then
            If (O_type) Then
               val =  x(a,i,b,j)*Occ_nums(a)
            Elseif (V_type) Then 
               val =  x(a,i,b,j)*(Done-Occ_nums(a))
            Endif 
         Elseif (Rit_one) Then 
            If (O_type) Then
               val =  x(a,i,b,j)**Occ_nums(j)
            Elseif (V_type) Then 
               val =  x(a,i,b,j)*(Done-Occ_nums(j))
            Endif 
         Endif 
      enddo
      enddo
      enddo
      enddo

      return
      end

      subroutine frac_occ42(x, a1,a2,atype,
     *                     i1,i2,itype,b1,b2,btype,
     *                     j1,j2, jtype, epsilon, epsilonb,sind)
      implicit none
      include 'interpreter.h'
      include 'frac_occ_vars.h'

      integer a1,a2,b1,b2,i1,i2,j1,j2
      integer atype, itype, btype, jtype
      double precision x(a1:a2,i1:i2,b1:b2,j1:j2)

      double precision epsilon(*), epsilonb(*)

      integer a,b,i,j,sind
      double precision eps, epsa, epsb, epsi, epsj
      double precision val,Done 
      Logical Lft_one,Rit_one
      Logical OO_type,OV_type,VO_type,VV_type

      Data Done /1.0D0/
  
      Lft_one = .False.
      Rit_one = .False.

      If (sind .eq. 1) Then
         Lft_one = .True.
      Elseif (sind .eq. 2) then
         Lft_one = .True.
      Endif 

      OO_type = .False.
      VV_type = .False.
      OV_type = .False.
      VO_type = .False.

      If (OO(1) .eq. 1 .and. OO(2) .eq. 1) Then 
         OO_type = .True.
      Elseif (VV(1) .eq. 1 .and. VV(2) .eq. 1) Then 
         VV_type = .True.
      Elseif (OO(1) .eq. 1 .and. VV(2) .eq. 1) Then 
         OV_type = .True.
      Elseif (VV(1) .eq. 1 .and. OO(2) .eq. 1) Then 
         VO_type = .True.
      Endif

      do j = j1,j2
      do b = b1,b2
      do i = i1,i2
      do a = a1,a2
         If (Lft_one) Then
            If (OO_type) Then
               val =  x(a,i,b,j)*Occ_nums(a)*Occ_nums(i)
            Elseif (OV_type) Then 
               val =  x(a,i,b,j)*Occ_nums(a)*(Done-Occ_nums(i))
            Elseif (VO_type) Then 
               val =  x(a,i,b,j)*(Done-Occ_nums(a))*Occ_nums(i)
            Elseif (VV_type) Then 
               val =  x(a,i,b,j)*(Done-Occ_nums(a))*(Done-Occ_nums(i))
            Endif 
         Elseif (Rit_one) Then 
            If (OO_type) Then
               val =  x(a,i,b,j)*Occ_nums(b)*Occ_nums(j)
            Elseif (OV_type) Then 
               val =  x(a,i,b,j)*Occ_nums(b)*(Done-Occ_nums(j))
            Elseif (VO_type) Then 
               val =  x(a,i,b,j)*(Done-Occ_nums(b))*Occ_nums(j)
            Elseif (VV_type) Then 
               val =  x(a,i,b,j)*(Done-Occ_nums(b))*(Done-Occ_nums(j))
            Endif 
         Endif 
      enddo
      enddo
      enddo
      enddo

      return
      end

      subroutine frac_occ43(x, a1,a2,atype,
     *                     i1,i2,itype,b1,b2,btype,
     *                     j1,j2, jtype, epsilon, epsilonb,sind)
      implicit none
      include 'interpreter.h'
      include 'frac_occ_vars.h'

      integer a1,a2,b1,b2,i1,i2,j1,j2
      integer atype, itype, btype, jtype
      double precision x(a1:a2,i1:i2,b1:b2,j1:j2)

      double precision epsilon(*), epsilonb(*)

      integer a,b,i,j,sind
      double precision eps, epsa, epsb, epsi, epsj
      double precision val,Done 
      Logical OOO_type,OVV_type,VVO_type,OOV_type,VOV_type,
     +        OVO_type,VVV_type,VOO_type
      Logical Lft_three,Rit_three 

      Data Done /1.0D0/
  
      Lft_three = .False.
      Rit_three = .False.

      If (sind .eq. 1) Then
         Lft_three = .True.
      Elseif (sind .eq. 2) then
         Rit_three = .True.
      Endif 
      OOO_type = .False.
      VVV_type = .False.
      OVO_type = .False.
      VOV_type = .False.
      VVO_type = .False.
      OVV_type = .False.
      OOV_type = .False.
      VOO_type = .False.

      If (OOO(1) .eq. 1 .and. OOO(2) .eq. 1 .and. OOO(3) .eq. 1) Then
         OOO_type = .True.
      ElseIf (VV(1) .eq. 1 .and. VV(2) .eq. 1 .and. VV(3) .eq. 1) Then
         VVV_type = .True.
      ElseIf (OO(1) .eq. 1 .and. OO(2) .eq. 1 .and. VV(3) .Eq. 1) Then 
         OOV_type = .True.
      ElseIf (OO(1) .eq. 1 .and. OO(2) .eq. 1 .and. VV(3) .Eq. 1) Then 
         VOO_type = .True.
      ElseIf (VV(1) .eq. 1 .and. VV(2) .eq. 1 .and. OO(3) .Eq. 1) Then 
         VVO_type = .True.
      ElseIf (OO(1) .eq. 1 .and. VV(2) .eq. 1 .and. OO(3) .Eq. 1) Then 
         OVO_type = .True.
      ElseIf (OO(1) .eq. 1 .and. VV(1) .eq. 1 .and. VV(3) .Eq. 1) Then 
         VOV_type = .True.
      ElseIf (OO(1) .eq. 1 .and. VV(2) .eq. 1 .and. VV(3) .Eq. 1) Then 
         OVV_type = .True.
      Endif 

      do j = j1,j2
      do b = b1,b2
      do i = i1,i2
      do a = a1,a2
         If (Lft_three) Then
            If (OOO_type) Then
               val =  x(a,i,b,j)*Occ_nums(a)*Occ_nums(i)*Occ_nums(b)
            Elseif (OVO_type) Then 
               val =  x(a,i,b,j)*Occ_nums(a)*(Done-Occ_nums(i))*
     +                           Occ_nums(b)
            Elseif (VOV_type) Then 
               val =  x(a,i,b,j)*(Done-Occ_nums(a))*Occ_nums(i)*
     +                           (Done-Occ_nums(b))
            Elseif (OOV_type) Then 
               val =  x(a,i,b,j)*Occ_nums(a)*Occ_nums(i)*
     +                           (Done-Occ_nums(b))
            Elseif (VOO_type) Then 
               val =  x(a,i,b,j)*(Done-Occ_nums(a))*Occ_nums(i)*
     +                            Occ_nums(b)
            Elseif (VVO_type) Then 
               val =  x(a,i,b,j)*(Done-Occ_nums(a))*(Done-Occ_nums(i))*
     +                           Occ_nums(b)
            Elseif (OVV_type) Then 
               val =  x(a,i,b,j)*Occ_nums(a)*(Done-Occ_nums(i))*
     +                           (Done-Occ_nums(b))
            Elseif (VVV_type) Then 
               val =  x(a,i,b,j)*(Done-Occ_nums(a))*(Done-Occ_nums(i))*
     +                           (Done-Occ_nums(b))
            Endif 
         Elseif (Rit_three) Then 
           If (OOO_type) Then
               val =  x(a,i,b,j)*Occ_nums(i)*Occ_nums(b)*Occ_nums(j)
            Elseif (OVO_type) Then
               val =  x(a,i,b,j)*Occ_nums(i)*(Done-Occ_nums(b))*
     +                           Occ_nums(j)
            Elseif (VOV_type) Then
               val =  x(a,i,b,j)*(Done-Occ_nums(i))*Occ_nums(b)*
     +                           (Done-Occ_nums(j))
            Elseif (OOV_type) Then
               val =  x(a,i,b,j)*Occ_nums(i)*Occ_nums(b)* 
     +                           (Done-Occ_nums(j))
            Elseif (VOO_type) Then
               val =  x(a,i,b,j)*(Done-Occ_nums(i))*Occ_nums(b)* 
     +                           Occ_nums(j)
            Elseif (VVO_type) Then
               val =  x(a,i,b,j)*(Done-Occ_nums(i))*(Done-Occ_nums(b))*
     +                           Occ_nums(j)
            Elseif (OVV_type) Then
               val =  x(a,i,b,j)*Occ_nums(i)*(Done-Occ_nums(b))*
     +                           (Done-Occ_nums(j))
            Elseif (VVV_type) Then
               val =  x(a,i,b,j)*(Done-Occ_nums(i))*(Done-Occ_nums(b))*
     +                           (Done-Occ_nums(j))
            Endif
         Endif 
      enddo
      enddo
      enddo
      enddo

      return
      end

      subroutine frac_occ44(x, a1,a2,atype,
     *                     i1,i2,itype,b1,b2,btype,
     *                     j1,j2, jtype, epsilon, epsilonb,sind)
      implicit none
      include 'interpreter.h'
      include 'frac_occ_vars.h'

      integer a1,a2,b1,b2,i1,i2,j1,j2
      integer atype, itype, btype, jtype
      double precision x(a1:a2,i1:i2,b1:b2,j1:j2)

      double precision epsilon(*), epsilonb(*)

      integer a,b,i,j,sind
      double precision eps, epsa, epsb, epsi, epsj
      double precision val,Done
      Logical  OVOV_type,VOVO_type

      Data Done /1.0D0/

      OVOV_type = .False.
      VOVO_type = .False.

      If (OOOO(1) .eq. 1 .and. VVVV(2) .eq. 1 .and. OOOO(3) .eq. 1
     +    .and. VVVV(4) .eq. 1) Then
          OVOV_type = .True.
      Elseif (VVVV(1) .eq. 1 .and. OOOO(2) .eq. 1 .and. VVVV(3) .eq. 1
     +    .and. OOOO(4) .eq. 1) then
         VOVO_type = .True.
      Endif

      do j = j1,j2
      do b = b1,b2
      do i = i1,i2
      do a = a1,a2
         If (OVOV_type) Then
            val =  x(a,i,b,j)*Occ_nums(a)*(Done-Occ_nums(i))*
     +                        Occ_nums(b)*(Done-Occ_nums(j))
         Elseif (VOVO_type) Then 
            val =  x(a,i,b,j)*(Done-Occ_nums(a))*Occ_nums(i)*
     +                        (Done-Occ_nums(b))*Occ_nums(j)
         Endif
      enddo
      enddo
      enddo
      enddo

      return
      end

      subroutine frac_occ21(x, a1,a2,atype,na1,na2,i1,i2,itype,
     *                     ni1,ni2,epsilon,epsilonb,sind)
      implicit none
      include 'interpreter.h'
      include 'frac_occ_vars.h'

      integer a1,a2,i1,i2
      integer atype, itype
      integer na1, na2, ni1, ni2,sind
      double precision x(na1:na2,ni1:ni2)
      double precision Done,epsilon,epsilonb

      integer a,i
      double precision eps, epsa, epsi
      double precision val
      Logical Lft_one,Rit_one
 
      Data Done /1.0D0/

      Lft_one = .False.
      Rit_one = .False.

      If (sind .eq. 1) Lft_one = .True.
      If (sind .eq. 2) Rit_one = .True.

      do i = i1, i2
      do a = a1,a2
         If (Lft_one) Then
            If (O .eq. 1) Then
               val =  x(a,i)*Occ_nums(a)
            ElseIf (V .eq. 1) then
               val =  x(a,i)*(Done-Occ_nums(a))
            Endif 
         Elseif (Rit_one) Then
            If (O .eq. 1) Then
               val =  x(a,i)*Occ_nums(i)
            ElseIf (V .eq. 1) then
               val =  x(a,i)*(Done-Occ_nums(i))
            Endif 
         Endif 
      enddo
      enddo

      return
      end

      subroutine frac_occ22(x, a1,a2,atype,na1,na2,i1,i2,itype,
     *                     ni1,ni2,epsilon,epsilonb,sind)
      implicit none
      include 'interpreter.h'
      include 'frac_occ_vars.h'

      integer a1,a2,i1,i2
      integer atype, itype
      integer na1, na2, ni1, ni2,sind
      double precision x(na1:na2,ni1:ni2)
      double precision Done,epsilon,epsilonb

      integer a,i
      double precision eps, epsa, epsi
      double precision val
      Logical OO_type,VV_type,OV_type,VO_type

      Data Done /1.0D0/

      OO_type = .False.
      VV_type = .False.
      VO_type = .False.
      OV_type = .False.

      If (OO(1) .eq. 1 .and. OO(2) .eq. 1) Then
         OO_type = .True.
      Elseif (VV(1) .eq. 1 .and. VV(2) .eq. 1) Then
         VV_type = .True.
      Elseif (OO(1) .eq. 1 .and. VV(2) .eq. 1) Then
         OV_type = .True.
      Elseif (VV(1) .eq. 1 .and. OO(2) .eq. 1) Then
         VO_type = .True.
      Endif

      do i = i1, i2
      do a = a1,a2
         If (OO_type) Then
            val =  x(a,i)*Occ_nums(a)*Occ_nums(i)
         Elseif (VV_type) Then
            val =  x(a,i)*(Done-Occ_nums(a))*(Done-Occ_nums(i))
         Elseif (VO_type) Then
            val =  x(a,i)*(Done-Occ_nums(a))*Occ_nums(i)
         Elseif (OV_type) Then
            val =  x(a,i)*Occ_nums(a)*(Done-Occ_nums(i))
         Endif 
      enddo
      enddo

      return
      end



