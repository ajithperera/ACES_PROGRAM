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

      sind = array_table(c_scalar_index, array)
      Arg  = Scalar_table(sind)

#ifdef _DEBUG_LVL0
      If (me.eq.0) write(6,"(a,i2)") "Arg: ", Arg
#endif 

      if (nindex .eq. 4 .and. contract .eq. 1) then
         call frac_occ41(x(indblk), val1(1),val2(1),type(1),
     *                 val1(2),val2(2),type(2),
     *                 val1(3),val2(3),type(3),
     *                 val1(4),val2(4),type(4), 
     *                 epsilon,epsilonb,Arg)
      elseif (nindex .eq. 4 .and. contract .eq. 2) then
        call frac_occ42(x(indblk), val1(1),val2(1),type(1),
     *                 val1(2),val2(2),type(2),
     *                 val1(3),val2(3),type(3),
     *                 val1(4),val2(4),type(4),
     *                 epsilon,epsilonb,Arg)
      elseif (nindex .eq. 4 .and. contract .eq. 3) then
        call frac_occ43(x(indblk), val1(1),val2(1),type(1),
     *                 val1(2),val2(2),type(2),
     *                 val1(3),val2(3),type(3),
     *                 val1(4),val2(4),type(4),
     *                 epsilon,epsilonb,Arg)
      elseif (nindex .eq. 4 .and. contract .eq. 4) then
        call frac_occ44(x(indblk), val1(1),val2(1),type(1),
     *                 val1(2),val2(2),type(2),
     *                 val1(3),val2(3),type(3),
     *                 val1(4),val2(4),type(4),
     *                 epsilon,epsilonb,Arg)
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
     *                 epsilon,epsilonb,Arg)    
         Elseif (Contract .Eq. 2) then
         call frac_occ22(x(indblk), val1(1),val2(1),type(1),
     *                 na1, na2,
     *                 val1(2),val2(2),type(2),ni1,ni2,
     *                 epsilon,epsilonb,Arg)    
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
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      integer a1,a2,b1,b2,i1,i2,j1,j2
      integer atype, itype, btype, jtype
      double precision x(a1:a2,i1:i2,b1:b2,j1:j2)

      double precision epsilon(*), epsilonb(*)

      integer a,b,i,j,sind
      double precision eps, epsa, epsb, epsi, epsj
      double precision val,Done 
      Logical Lft_one,Rit_one
      Logical O_type,V_type
      Logical AL,BL,AR,BR

      Data Done /1.0D0/
  
      If (sind .eq. 1) Then
         Lft_one = .True.
      Elseif (sind .eq. 2) then
         Rit_one = .True.
      Endif 

      O_type = .False.
      V_type = .False.

      AL     = .False.
      BL     = .False.
      AR     = .False.
      BR     = .False.

      If (O .eq. 1) Then 
         O_type = .True.
      ElseIf (V .eq. 1) Then 
         V_type = .True.
      Endif
      
      If (Lft_one) Then
      AL = Atype .eq. Moaindex 
      BL = Atype .eq. Mobindex 
      Endif 

      If (Rit_one) Then
      AR = Jtype .eq. Moaindex 
      BR = Jtype .eq. Mobindex 
      Endif

#ifdef _DEBUG_LVL0
      if (me .eq. 0) then
         write(6,"(4l)")AL,BL,AR,BR
         write(6,"(4l)")o_type,v_type,o,v
         Write(6,"(a)") "The Occ_num arrays"
         Write(6,"(6F4.1)") (Occ_numsa(i),i=1, Nbasis)
         Write(6,"(6F4.1)") (Occ_numsb(i),i=1, Nbasis)
      endif
#endif
      do j = j1,j2
      do b = b1,b2
      do i = i1,i2
      do a = a1,a2
         If (Lft_one) Then
            If (O_type) Then
               If (AL) x(a,i,b,j) =  x(a,i,b,j)*Occ_numsa(a)
               If (BL) x(a,i,b,j) =  x(a,i,b,j)*Occ_numsb(a)
            Elseif (V_type) Then 
               If (AL) x(a,i,b,j) =  x(a,i,b,j)*(Done-Occ_numsa(a))
               If (BL) x(a,i,b,j) =  x(a,i,b,j)*(Done-Occ_numsb(a))
            Endif 
         Elseif (Rit_one) Then 
            If (O_type) Then
               If (AR) x(a,i,b,j) =  x(a,i,b,j)*Occ_numsa(j)
               If (BR) x(a,i,b,j) =  x(a,i,b,j)*Occ_numsb(j)
            Elseif (V_type) Then 
               If (AR) x(a,i,b,j) =  x(a,i,b,j)*(Done-Occ_numsa(j))
               If (BR) x(a,i,b,j) =  x(a,i,b,j)*(Done-Occ_numsb(j))
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
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      integer a1,a2,b1,b2,i1,i2,j1,j2
      integer atype, itype, btype, jtype
      double precision x(a1:a2,i1:i2,b1:b2,j1:j2)

      double precision epsilon(*), epsilonb(*)

      integer a,b,i,j,sind
      double precision eps, epsa, epsb, epsi, epsj
      double precision val,Done 
      Logical Lft_two,Rit_two
      Logical OO_type,OV_type,VO_type,VV_type
      Logical AAL,BBL,ABL,BAL,AAR,BBR,ABR,BAR

      Data Done /1.0D0/
  
      Lft_two = .False.
      Rit_two = .False.

      If (sind .eq. 1) Then
         Lft_two = .True.
      Elseif (sind .eq. 2) then
         Rit_two = .True.
      Endif 

      OO_type = .False.
      VV_type = .False.
      OV_type = .False.
      VO_type = .False.
      AAL     = .False.
      BBL     = .False. 
      ABL     = .False.
      BAL     = .False.
      AAR     = .False.
      BBR     = .False.
      ABR     = .False.
      BAR     = .False.

#ifdef _DEBUG_LVL0
      If(me .eq.0) then
       write(6,"(a,4I2)") "OO and VV :",OO(1),OO(2),VV(1),VV(2)
      Endif
#endif

      If (OO(1) .eq. 1 .and. OO(2) .eq. 1) Then 
         OO_type = .True.
      Elseif (VV(1) .eq. 1 .and. VV(2) .eq. 1) Then 
         VV_type = .True.
      Elseif (OO(1) .eq. 1 .and. VV(2) .eq. 1) Then 
         OV_type = .True.
      Elseif (VV(1) .eq. 1 .and. OO(2) .eq. 1) Then 
         VO_type = .True.
      Endif

      If (Lft_two) Then
      AAL = Atype .eq. Moaindex .and. Itype .eq. Moaindex
      BBL = Atype .eq. Mobindex .and. Itype .eq. Mobindex
      ABL = Atype .eq. Moaindex .and. Itype .eq. Mobindex
      BAL = Atype .eq. Mobindex .and. Itype .eq. Moaindex
      Endif 

      If (Rit_two) Then
      AAR = Btype .eq. Moaindex .and. Jtype .eq. Moaindex
      BBR = Btype .eq. Mobindex .and. Jtype .eq. Mobindex
      ABR = Btype .eq. Moaindex .and. Jtype .eq. Mobindex
      BAR = Btype .eq. Mobindex .and. Jtype .eq. Moaindex
      Endif 

#ifdef _DEBUG_LVL0
      if (me .eq. 0) then
C      write(6,"(a,3l)") "The AO-Ladder spin: ",AO_AA,AO_BB,AO_AB
      write(6,"(a,a,6l)")"Lft_two,Rit_two,OO_type,VV_type,OV_type,",
     +                   "VO_type,AAL-BBR",Lft_two,Rit_two,OO_type,
     +                    VV_type,OV_type,VO_type
      write(6,"(a,8l)") "type:",aal,bbl,abl,bal,aar,bbr,abr,bar
      Write(6,"(a)") "The Occ_num arrays"
      Write(6,"(6F4.1)") (Occ_numsa(i),i=1, Nbasis)
      Write(6,"(6F4.1)") (Occ_numsb(i),i=1, Nbasis)
      endif
#endif 
    
      do j = j1,j2
      do b = b1,b2
      do i = i1,i2
      do a = a1,a2
         If (Lft_two) Then
            If (OO_type) Then
               If(AAL) x(a,i,b,j) =  x(a,i,b,j)*Occ_numsa(a)*
     +                                          Occ_numsa(i)
               If(BBL) x(a,i,b,j) =  x(a,i,b,j)*Occ_numsb(a)*
     +                                          Occ_numsb(i)
               If(ABL) x(a,i,b,j) =  x(a,i,b,j)*Occ_numsa(a)*
     +                                          Occ_numsb(i)
               If(BAL) x(a,i,b,j) =  x(a,i,b,j)*Occ_numsb(a)*
     +                                          Occ_numsa(i)
            Elseif (OV_type) Then 
               If (AAL) x(a,i,b,j) =  x(a,i,b,j)*Occ_numsa(a)*
     +                                          (Done-Occ_numsb(i))
               If (BBL) x(a,i,b,j) =  x(a,i,b,j)*Occ_numsb(a)*
     +                                          (Done-Occ_numsb(i))
               If (ABL) x(a,i,b,j) =  x(a,i,b,j)*Occ_numsa(a)*
     +                                          (Done-Occ_numsb(i))
               If (BAL) x(a,i,b,j) =  x(a,i,b,j)*Occ_numsb(a)*
     +                                          (Done-Occ_numsa(i))
            Elseif (VO_type) Then 
               If (AAL) x(a,i,b,j) =  x(a,i,b,j)*(Done-Occ_numsa(a))*
     +                                            Occ_numsa(i)
               If (BBL) x(a,i,b,j) =  x(a,i,b,j)*(Done-Occ_numsb(a))*
     +                                            Occ_numsb(i)
               If (ABL) x(a,i,b,j) =  x(a,i,b,j)*(Done-Occ_numsa(a))*
     +                                            Occ_numsb(i)
               If (BAL) x(a,i,b,j) =  x(a,i,b,j)*(Done-Occ_numsb(a))*
     +                                            Occ_numsa(i)
            Elseif (VV_type) Then  
               If (AAL) x(a,i,b,j) =  x(a,i,b,j)*(Done-Occ_numsa(a))*
     +                                           (Done-Occ_numsa(i))
               If (BBL) x(a,i,b,j) =  x(a,i,b,j)*(Done-Occ_numsb(a))*
     +                                           (Done-Occ_numsb(i))
               If (ABL) x(a,i,b,j) =  x(a,i,b,j)*(Done-Occ_numsa(a))*
     +                                           (Done-Occ_numsb(i))
               If (BAL) x(a,i,b,j) =  x(a,i,b,j)*(Done-Occ_numsb(a))*
     +                                           (Done-Occ_numsa(i))
            Endif 
         Elseif (Rit_two) Then 
            If (OO_type) Then
               If(AAR) x(a,i,b,j) =  x(a,i,b,j)*Occ_numsa(b)*
     +                                          Occ_numsa(j)
               If(BBR) x(a,i,b,j) =  x(a,i,b,j)*Occ_numsb(b)*
     +                                          Occ_numsb(j)
               If(ABR) x(a,i,b,j) =  x(a,i,b,j)*Occ_numsa(b)*
     +                                          Occ_numsb(j)
               If(BAR) x(a,i,b,j) =  x(a,i,b,j)*Occ_numsb(b)*
     +                                          Occ_numsa(j)
            Elseif (OV_type) Then 
               If (AAR) x(a,i,b,j) =  x(a,i,b,j)*Occ_numsa(b)*
     +                                          (Done-Occ_numsb(j))
               If (BBR) x(a,i,b,j) =  x(a,i,b,j)*Occ_numsb(b)*
     +                                          (Done-Occ_numsb(j))
               If (ABR) x(a,i,b,j) =  x(a,i,b,j)*Occ_numsa(b)*
     +                                          (Done-Occ_numsb(j))
               If (BAR) x(a,i,b,j) =  x(a,i,b,j)*Occ_numsb(b)*
     +                                          (Done-Occ_numsa(j))
            Elseif (VO_type) Then 
               If (AAR) x(a,i,b,j) =  x(a,i,b,j)*(Done-Occ_numsa(b))*
     +                                            Occ_numsa(j)
               If (BBR) x(a,i,b,j) =  x(a,i,b,j)*(Done-Occ_numsb(b))*
     +                                            Occ_numsb(j)
               If (ABR) x(a,i,b,j) =  x(a,i,b,j)*(Done-Occ_numsa(b))*
     +                                            Occ_numsb(j)
               If (BAR) x(a,i,b,j) =  x(a,i,b,j)*(Done-Occ_numsb(b))*
     +                                            Occ_numsa(j)
            Elseif (VV_type) Then 
               If (AAR) x(a,i,b,j) =  x(a,i,b,j)*(Done-Occ_numsa(b))*
     +                                           (Done-Occ_numsa(j))
               If (BBR) x(a,i,b,j) =  x(a,i,b,j)*(Done-Occ_numsb(b))*
     +                                           (Done-Occ_numsb(j))
               If (ABR) x(a,i,b,j) =  x(a,i,b,j)*(Done-Occ_numsa(b))*
     +                                           (Done-Occ_numsb(j))
               If (BAR) x(a,i,b,j) =  x(a,i,b,j)*(Done-Occ_numsb(b))*
     +                                           (Done-Occ_numsa(j))
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
      include 'parallel_info.h'
      include 'int_gen_parms.h'

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
      Logical AAAL,BBBL,ABAL,BABL,AAAR,BBBR,ABAR,BABR,ABBL,ABBR
      Logical AABR,AABL

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
      AAAL     = .False.
      BBBL     = .False.
      ABAL     = .False.
      BABL     = .False.
      ABBL     = .False.
      AABL     = .False.
      AAAR     = .False.
      BBBR     = .False.
      ABAR     = .False.
      BABR     = .False.
      ABBR     = .False.
      AABR     = .False.

      If (OOO(1) .eq. 1 .and. OOO(2) .eq. 1 .and. OOO(3) .eq. 1) Then
         OOO_type = .True.
      ElseIf (VVV(1) .eq. 1 .and. VVV(2) .eq. 1 .and. VVV(3) 
     +        .eq. 1) Then
         VVV_type = .True.
      ElseIf (OOO(1) .eq. 1 .and. OOO(2) .eq. 1 .and. VVV(3)
     +        .Eq. 1) Then 
         OOV_type = .True.
      ElseIf (OOO(2) .eq. 1 .and. OOO(3) .eq. 1 .and. VVV(1)
     +        .Eq. 1) Then 
         VOO_type = .True.
      ElseIf (VVV(1) .eq. 1 .and. VVV(2) .eq. 1 .and. OOO(3) 
     +        .Eq. 1) Then 
         VVO_type = .True.
      ElseIf (OOO(1) .eq. 1 .and. VVV(2) .eq. 1 .and. OOO(3)
     +        .Eq. 1) Then 
         OVO_type = .True.
      ElseIf (OOO(2) .eq. 1 .and. VVV(1) .eq. 1 .and. VVV(3)
     +        .Eq. 1) Then 
         VOV_type = .True.
      ElseIf (OOO(1) .eq. 1 .and. VVV(2) .eq. 1 .and. VVV(3)
     +        .Eq. 1) Then 
         OVV_type = .True.
      Endif 

      If (Lft_three) Then
      AAAL = Atype .eq. Moaindex .and. Itype .eq. Moaindex .and. 
     +       Btype .eq. Moaindex 
      BBBL = Atype .eq. Mobindex .and. Itype .eq. Mobindex .and. 
     +       Btype .eq. Mobindex 
      ABAL = Atype .eq. Moaindex .and. Itype .eq. Mobindex .and. 
     +       Btype .eq. Moaindex 
      BABL = Atype .eq. Mobindex .and. Itype .eq. Moaindex .and. 
     +       Btype .eq. Mobindex 
      ABBL = Atype .eq. Moaindex .and. Itype .eq. Mobindex .and. 
     +       Btype .eq. Mobindex 
      AABL = Atype .eq. Moaindex .and. Itype .eq. Moaindex .and.
     +       Btype .eq. Mobindex
      Endif 

      If (Rit_three) Then
      AAAR = Itype .eq. Moaindex .and. Btype .eq. Moaindex .and.
     +       Jtype .eq. Moaindex
      BBBR = Itype .eq. Mobindex .and. Btype .eq. Mobindex .and.
     +       Jtype .eq. Mobindex
      ABAR = Itype .eq. Moaindex .and. Btype .eq. Mobindex .and.
     +       Jtype .eq. Moaindex
      BABR = Itype .eq. Mobindex .and. Btype .eq. Moaindex .and.
     +       Jtype .eq. Mobindex
      ABBR = Itype .eq. Moaindex .and. Btype .eq. Mobindex .and. 
     +       Jtype .eq. Mobindex 
      AABR = Itype .eq. Moaindex .and. Btype .eq. Moaindex .and.
     +       Jtype .eq. Mobindex
      Endif 

#ifdef _DEBUG_LVL0
      if (me .eq. 0) then
          write(6,"(a,4(1x,i4))") "orbital counts", nalpha_occupied,
     +    nbeta_occupied,nalpha_virtual,nbeta_virtual 
          write(6,"(12l)")AAAL,BBBL,ABAL,BABL,ABBL,AABL,AAAR,BBBR,
     +                    ABAR,BABR,ABBR,AABR
          write(6,"(8l)")ooo_type,vvv_type,oov_type,voo_type,
     +                   vvo_type,ovo_type,vov_type,ovv_type

      endif 
#endif
      Return 
      do j = j1,j2
      do b = b1,b2
      do i = i1,i2
      do a = a1,a2
         If (Lft_three) Then
            If (OOO_type) Then
               If (AAAL) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(a)*
     +                                           Occ_numsa(i)*
     +                                           Occ_numsa(b)
               If (BBBL) x(a,i,b,j) = x(a,i,b,j)*Occ_numsb(a)*
     +                                           Occ_numsb(i)*
     +                                           Occ_numsb(b)
               If (ABAL) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(a)*
     +                                           Occ_numsb(i)*
     +                                           Occ_numsa(b)
               If (BABL) x(a,i,b,j) = x(a,i,b,j)*Occ_numsb(a)*
     +                                           Occ_numsa(i)*
     +                                           Occ_numsb(b)
               If (ABBL) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(a)*
     +                                           Occ_numsb(i)*
     +                                           Occ_numsb(b)
               If (AABL) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(a)*
     +                                           Occ_numsa(i)*
     +                                           Occ_numsb(b)
            Elseif (OVO_type) Then 
               If (AAAL) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(a)*
     +                                           (Done-Occ_numsa(i))*
     +                                           Occ_numsa(b)
               If (BBBL) x(a,i,b,j) = x(a,i,b,j)*Occ_numsb(a)*
     +                                           (Done-Occ_numsb(i))*
     +                                           Occ_numsb(b)
               If (ABAL) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(a)*
     +                                           (Done-Occ_numsb(i))*
     +                                           Occ_numsa(b)
               If (BABL) x(a,i,b,j) = x(a,i,b,j)*Occ_numsb(a)*
     +                                           (Done-Occ_numsa(i))*
     +                                           Occ_numsb(b)
               If (ABBL) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(a)*
     +                                           (Done-Occ_numsb(i))*
     +                                           Occ_numsb(b)
               If (AABL) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(a)*
     +                                           (Done-Occ_numsa(i))*
     +                                           Occ_numsb(b)
            Elseif (VOV_type) Then 
               If (AAAL) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(a))*
     +                                Occ_numsa(i)*(Done-Occ_numsa(b))
               If (BBBL) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsb(a))*
     +                                Occ_numsb(i)*(Done-Occ_numsb(b))
               If (ABAL) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(a))*
     +                                Occ_numsb(i)*(Done-Occ_numsa(b))
               If (BABL) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsb(a))*
     +                                Occ_numsa(i)*(Done-Occ_numsb(b))
               If (ABBL) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(a))*
     +                                Occ_numsb(i)*(Done-Occ_numsb(b))
               If (AABL) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(a))*
     +                                Occ_numsa(i)*(Done-Occ_numsb(b))
            Elseif (OOV_type) Then 
               If (AAAL) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(a)*
     +                                Occ_numsa(i)*(Done-Occ_numsa(b))
               If (AAAL) x(a,i,b,j) = x(a,i,b,j)*Occ_numsb(a)*
     +                                Occ_numsb(i)*(Done-Occ_numsb(b))
               If (ABAL) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(a)*
     +                                Occ_numsb(i)*(Done-Occ_numsa(b))
               If (BABL) x(a,i,b,j) = x(a,i,b,j)*Occ_numsb(a)*
     +                                Occ_numsa(i)*(Done-Occ_numsb(b))
               If (ABBL) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(a)*
     +                                Occ_numsb(i)*(Done-Occ_numsb(b))
               If (AABL) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(a)*
     +                                Occ_numsa(i)*(Done-Occ_numsb(b))
            Elseif (VOO_type) Then 
               If (AAAL) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(a))*
     +                                Occ_numsa(i)*Occ_numsa(b)
               If (BBBL) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsb(a))*
     +                                Occ_numsb(i)*Occ_numsb(b)
               If (ABAL) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(a))*
     +                                Occ_numsb(i)*Occ_numsa(b)
               If (BABL) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsb(a))*
     +                                Occ_numsa(i)*Occ_numsb(b)
               If (BABL) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(a))*
     +                                Occ_numsb(i)*Occ_numsb(b)
               If (AABL) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(a))*
     +                                Occ_numsa(i)*Occ_numsb(b)
            Elseif (VVO_type) Then 
               If (AAAL) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(a))*
     +                                (Done-Occ_numsa(i))*Occ_numsa(b)
               If (BBBL) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsb(a))*
     +                                (Done-Occ_numsb(i))*Occ_numsb(b)
               If (ABAL) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(a))*
     +                                (Done-Occ_numsb(i))*Occ_numsa(b)
               If (BABL) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsb(a))*
     +                                (Done-Occ_numsa(i))*Occ_numsb(b)
               If (BABL) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(a))*
     +                                (Done-Occ_numsb(i))*Occ_numsb(b)
               If (AABL) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(a))*
     +                                (Done-Occ_numsa(i))*Occ_numsb(b)
            Elseif (OVV_type) Then 
               If (AAAL) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(a)*
     +                                (Done-Occ_numsa(i))*
     +                                (Done-Occ_numsa(b))
               If (BBBL) x(a,i,b,j) = x(a,i,b,j)*Occ_numsb(a)*
     +                                (Done-Occ_numsb(i))*
     +                                (Done-Occ_numsb(b))
               If (ABAL) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(a)*
     +                                (Done-Occ_numsb(i))*
     +                                (Done-Occ_numsa(b))
               If (BABL) x(a,i,b,j) = x(a,i,b,j)*Occ_numsb(a)*
     +                                (Done-Occ_numsa(i))*
     +                                (Done-Occ_numsb(b))
               If (ABBL) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(a)*
     +                                (Done-Occ_numsb(i))*
     +                                (Done-Occ_numsb(b))
               If (AABL) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(a)*
     +                                (Done-Occ_numsa(i))*
     +                                (Done-Occ_numsb(b))
            Elseif (VVV_type) Then 
               If (AAAL) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(a))*
     +                                (Done-Occ_numsa(i))*
     +                                (Done-Occ_numsa(b))
               If (BBBL) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsb(a))*
     +                                (Done-Occ_numsb(i))*
     +                                (Done-Occ_numsb(b))
               If (ABAL) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(a))*
     +                                (Done-Occ_numsb(i))*
     +                                (Done-Occ_numsa(b))
               If (BABL) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsb(a))*
     +                                (Done-Occ_numsa(i))*
     +                                (Done-Occ_numsb(b))
               If (ABBL) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(a))*
     +                                (Done-Occ_numsb(i))*
     +                                (Done-Occ_numsb(b))
               If (AABL) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(a))*
     +                                (Done-Occ_numsa(i))*
     +                                (Done-Occ_numsb(b))
            Endif 
         Elseif (Rit_three) Then 
           If (OOO_type) Then
              If (AAAR) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(i)*
     +                                          Occ_numsa(b)*
     +                                          Occ_numsa(j)
              If (BBBR) x(a,i,b,j) = x(a,i,b,j)*Occ_numsb(i)*
     +                                          Occ_numsb(b)*
     +                                          Occ_numsb(j)
              If (ABAR) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(i)*
     +                                          Occ_numsb(b)*
     +                                          Occ_numsa(j)
              If (BABR) x(a,i,b,j) = x(a,i,b,j)*Occ_numsb(i)*
     +                                          Occ_numsa(b)*
     +                                          Occ_numsb(j)
              If (ABBR) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(i)*
     +                                          Occ_numsb(b)*
     +                                          Occ_numsb(j)
              If (ABBR) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(i)*
     +                                          Occ_numsa(b)*
     +                                          Occ_numsb(j)
            Elseif (OVO_type) Then
              If (AAAR) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(i)*
     +                                          (Done-Occ_numsa(b))*
     +                                          Occ_numsa(j)
              If (BBBR) x(a,i,b,j) = x(a,i,b,j)*Occ_numsb(i)*
     +                                          (Done-Occ_numsb(b))*
     +                                          Occ_numsb(j)
              If (ABAR) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(i)*
     +                                          (Done-Occ_numsb(b))*
     +                                          Occ_numsa(j)
              If (BABR) x(a,i,b,j) = x(a,i,b,j)*Occ_numsb(i)*
     +                                          (Done-Occ_numsa(b))*
     +                                          Occ_numsb(j)
              If (ABBR) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(i)*
     +                                          (Done-Occ_numsb(b))*
     +                                          Occ_numsb(j)
              If (AABR) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(i)*
     +                                          (Done-Occ_numsa(b))*
     +                                          Occ_numsb(j)
            Elseif (VOV_type) Then
               If (AAAR) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(i))*
     +                                Occ_numsa(b)*(Done-Occ_numsa(j))
               If (BBBR) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsb(i))*
     +                                Occ_numsb(b)*(Done-Occ_numsb(j))
               If (ABAR) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(i))*
     +                                Occ_numsb(b)*(Done-Occ_numsa(j))
               If (BABR) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsb(i))*
     +                                Occ_numsa(b)*(Done-Occ_numsb(j))
               If (ABBR) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(i))*
     +                                Occ_numsb(b)*(Done-Occ_numsb(j))
               If (AABR) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(i))*
     +                                Occ_numsa(b)*(Done-Occ_numsb(j))
            Elseif (OOV_type) Then
               If (AAAR) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(i)*
     +                                Occ_numsa(b)*(Done-Occ_numsa(j))
               If (AAAR) x(a,i,b,j) = x(a,i,b,j)*Occ_numsb(i)*
     +                                Occ_numsb(b)*(Done-Occ_numsb(j))
               If (ABAR) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(i)*
     +                                Occ_numsb(b)*(Done-Occ_numsa(j))
               If (BABR) x(a,i,b,j) = x(a,i,b,j)*Occ_numsb(i)*
     +                                Occ_numsa(b)*(Done-Occ_numsb(j))
               If (ABBR) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(i)*
     +                                Occ_numsb(b)*(Done-Occ_numsb(j))
               If (AABR) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(i)*
     +                                Occ_numsa(b)*(Done-Occ_numsb(j))
            Elseif (VOO_type) Then
               If (AAAR) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(i))*
     +                                Occ_numsa(b)*Occ_numsa(j)
               If (BBBR) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsb(i))*
     +                                Occ_numsb(b)*Occ_numsb(j)
               If (ABAR) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsb(i))*
     +                                Occ_numsa(b)*Occ_numsa(j)
               If (BABR) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsb(i))*
     +                                Occ_numsa(b)*Occ_numsb(j)
               If (ABBR) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(i))*
     +                                Occ_numsb(b)*Occ_numsb(j)
               If (AABR) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(i))*
     +                                Occ_numsa(b)*Occ_numsb(j)
            Elseif (VVO_type) Then
               If (AAAR) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(i))*
     +                                (Done-Occ_numsa(b))*Occ_numsa(j)
               If (BBBR) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsb(i))*
     +                                (Done-Occ_numsb(b))*Occ_numsb(j)
               If (ABAR) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(i))*
     +                                (Done-Occ_numsb(b))*Occ_numsa(j)
               If (BABR) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsb(i))*
     +                                (Done-Occ_numsa(b))*Occ_numsb(j)
               If (ABBR) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(i))*
     +                                (Done-Occ_numsb(b))*Occ_numsb(j)
               If (AABR) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(i))*
     +                                (Done-Occ_numsa(b))*Occ_numsb(j)
            Elseif (OVV_type) Then
               If (AAAR) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(i)*
     +                                (Done-Occ_numsa(b))*
     +                                (Done-Occ_numsa(j))
               If (BBBR) x(a,i,b,j) = x(a,i,b,j)*Occ_numsb(i)*
     +                                (Done-Occ_numsb(b))*
     +                                (Done-Occ_numsb(j))
               If (ABAR) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(i)*
     +                                (Done-Occ_numsb(b))*
     +                                (Done-Occ_numsa(j))
               If (BABR) x(a,i,b,j) = x(a,i,b,j)*Occ_numsb(i)*
     +                                (Done-Occ_numsa(b))*
     +                                (Done-Occ_numsb(j))
               If (ABBR) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(i)*
     +                                (Done-Occ_numsb(b))*
     +                                (Done-Occ_numsb(j))
               If (AABR) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(i)*
     +                                (Done-Occ_numsa(b))*
     +                                (Done-Occ_numsb(j))
            Elseif (VVV_type) Then
               If (AAAR) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(i))*
     +                                (Done-Occ_numsa(b))*
     +                                (Done-Occ_numsa(j))
               If (BBBR) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsb(i))*
     +                                (Done-Occ_numsb(b))*
     +                                (Done-Occ_numsb(j))
               If (ABAR) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(i))*
     +                                (Done-Occ_numsb(b))*
     +                                (Done-Occ_numsa(j))
               If (BABR) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsb(i))*
     +                                (Done-Occ_numsa(b))*
     +                                (Done-Occ_numsb(j))
               If (ABBR) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(i))*
     +                                (Done-Occ_numsb(b))*
     +                                (Done-Occ_numsb(j))
               If (AABR) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(i))*
     +                                (Done-Occ_numsa(b))*
     +                                (Done-Occ_numsb(j))
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
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      integer a1,a2,b1,b2,i1,i2,j1,j2
      integer atype, itype, btype, jtype
      double precision x(a1:a2,i1:i2,b1:b2,j1:j2)

      double precision epsilon(*), epsilonb(*)

      integer a,b,i,j,sind
      double precision eps, epsa, epsb, epsi, epsj
      double precision val,Done
      Logical OVOV_type,VOVO_type
      Logical AAAA,BBBB,AABB,BBAA,ABAB,BABA

      Data Done /1.0D0/

      OVOV_type = .False.
      VOVO_type = .False.
      AAAA      = .FALSE.
      BBBB      = .FALSE.
      AABB      = .FALSE.
      BBAA      = .FALSE.
      ABAB      = .FALSE.
      BABA      = .FALSE.

      If (OOOO(1) .eq. 1 .and. VVVV(2) .eq. 1 .and. OOOO(3) .eq. 1
     +    .and. VVVV(4) .eq. 1) Then
          OVOV_type = .True.
      Elseif (VVVV(1) .eq. 1 .and. OOOO(2) .eq. 1 .and. VVVV(3) .eq. 1
     +    .and. OOOO(4) .eq. 1) then
         VOVO_type = .True.
      Endif

      AAAA = Atype .eq. Moaindex .and. Itype .eq. Moaindex .and.
     +       Btype .eq. Moaindex .and. Jtype .eq. Moaindex 

      BBBB = Atype .eq. Mobindex .and. Itype .eq. Mobindex .and.
     +       Btype .eq. Mobindex .and. Jtype .eq. Mobindex 

      AABB = Atype .eq. Moaindex .and. Itype .eq. Moaindex .and.
     +       Btype .eq. Mobindex .and. Jtype .eq. Mobindex 

      BBAA = Atype .eq. Mobindex .and. Itype .eq. Mobindex .and.
     +       Btype .eq. Moaindex .and. Jtype .eq. Moaindex 

      ABAB = Atype .eq. Moaindex .and. Itype .eq. Mobindex .and.
     +       Btype .eq. Moaindex .and. Jtype .eq. Mobindex 

      BABA = Atype .eq. Mobindex .and. Itype .eq. Moaindex .and.
     +       Btype .eq. Mobindex .and. Jtype .eq. Moaindex 

#ifdef _DEBUG_LVL0
      If (me .eq. 0) Then 
         write(6,"(a,2l)") "The type:",OVOV_type,VOVO_type
         write(6,"(a,6l)") "AAAA,BBBB,ABAB,BABA,ABAB,BABA:",
     +                      AAAA,BBBB,AABB,BBAA,ABAB,BABA
         Write(6,"(a)") "The Occ_num arrays"
         Write(6,"(6F4.1)") (Occ_numsa(i),i=1, Nbasis)
         Write(6,"(6F4.1)") (Occ_numsb(i),i=1, Nbasis)
      Endif 
#endif 

      do j = j1,j2
      do b = b1,b2
      do i = i1,i2
      do a = a1,a2
         If (OVOV_type) Then
            If (AAAA) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(a)*
     +                             (Done-Occ_numsa(i))*
     +                             Occ_numsa(b)*(Done-Occ_numsa(j))
            If (BBBB) x(a,i,b,j) = x(a,i,b,j)*Occ_numsb(a)*
     +                             (Done-Occ_numsb(i))*
     +                             Occ_numsb(b)*(Done-Occ_numsb(j))
            If (AABB) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(a)*
     +                             (Done-Occ_numsa(i))*
     +                             Occ_numsb(b)*(Done-Occ_numsb(j))
            If (BBAA) x(a,i,b,j) = x(a,i,b,j)*Occ_numsb(a)*
     +                             (Done-Occ_numsb(i))*
     +                             Occ_numsa(b)*(Done-Occ_numsa(j))
            If (ABAB) x(a,i,b,j) = x(a,i,b,j)*Occ_numsa(a)*
     +                             (Done-Occ_numsb(i))*
     +                             Occ_numsa(b)*(Done-Occ_numsb(j))
            If (BABA) x(a,i,b,j) = x(a,i,b,j)*Occ_numsb(a)*
     +                             (Done-Occ_numsa(i))*
     +                             Occ_numsb(b)*(Done-Occ_numsa(j))
         Elseif (VOVO_type) Then 
            If (AAAA) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(a))*
     +                             Occ_numsa(i)*(Done-Occ_numsa(b))*
     +                             Occ_numsa(j)
            If (BBBB) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsb(a))*
     +                             Occ_numsb(i)*(Done-Occ_numsb(b))*
     +                             Occ_numsb(j)
            If (AABB) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(a))*
     +                             Occ_numsa(i)*(Done-Occ_numsb(b))*
     +                             Occ_numsb(j)
            If (BBAA) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsb(a))*
     +                             Occ_numsb(i)*(Done-Occ_numsa(b))*
     +                             Occ_numsa(j)
            If (ABAB) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsa(a))*
     +                             Occ_numsb(i)*(Done-Occ_numsa(b))*
     +                             Occ_numsb(j)
            If (BABA) x(a,i,b,j) = x(a,i,b,j)*(Done-Occ_numsb(a))*
     +                             Occ_numsa(i)*(Done-Occ_numsb(b))*
     +                             Occ_numsa(j)
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
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      integer a1,a2,i1,i2
      integer atype, itype
      integer na1, na2, ni1, ni2,sind
      double precision x(na1:na2,ni1:ni2)
      double precision Done,epsilon,epsilonb

      integer a,i
      double precision eps, epsa, epsi
      double precision val
      Logical Lft_one,Rit_one
      Logical AL,BL,AR,BR
 
      Data Done /1.0D0/

      Lft_one = .False.
      Rit_one = .False.
      AL      = .False.
      BL      = .False.
      AR      = .False.
      BR      = .False.
    
      If (sind .eq. 1) Lft_one = .True.
      If (sind .eq. 2) Rit_one = .True.

      If (Lft_one) Then
      AL = Atype .eq. Moaindex 
      BL = Atype .eq. Mobindex
      Endif 

      If (Rit_one) Then
      AR = Itype .eq. Moaindex 
      BR = Itype .eq. Mobindex 
      Endif 
#ifdef _DEBUG_LVL0
      if (me .eq. 0) then
         write(6,"(4l)")AL,BL,AR,BR
         write(6,"(2l)")o,v
         Write(6,"(a)") "The Occ_num arrays"
         Write(6,"(6F4.1)") (Occ_numsa(i),i=1, Nbasis)
         Write(6,"(6F4.1)") (Occ_numsb(i),i=1, Nbasis)
      endif
#endif

      do i = i1, i2
      do a = a1,a2
         If (Lft_one) Then
            If (O .eq. 1) Then
               If (AL) x(a,i) = x(a,i)*Occ_numsa(a)
               If (BL) x(a,i) = x(a,i)*Occ_numsb(a)
            ElseIf (V .eq. 1) then
               If (AL) x(a,i) =  x(a,i)*(Done-Occ_numsa(a))
               If (BL) x(a,i) =  x(a,i)*(Done-Occ_numsb(a))
            Endif 
         Elseif (Rit_one) Then
            If (O .eq. 1) Then
               If (AR) x(a,i) = x(a,i)*Occ_numsa(i)
               If (BR) x(a,i) = x(a,i)*Occ_numsb(i)
            ElseIf (V .eq. 1) then
               If (AR) x(a,i) =  x(a,i)*(Done-Occ_numsa(i))
               If (BR) x(a,i) =  x(a,i)*(Done-Occ_numsb(i))
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
      Logical AA,BB,AB,BA

      Data Done /1.0D0/

      OO_type = .False.
      VV_type = .False.
      VO_type = .False.
      OV_type = .False.
      AA      = .False.
      BB      = .False.
      AB      = .False.
      BA      = .False.

      If (OO(1) .eq. 1 .and. OO(2) .eq. 1) Then
         OO_type = .True.
      Elseif (VV(1) .eq. 1 .and. VV(2) .eq. 1) Then
         VV_type = .True.
      Elseif (OO(1) .eq. 1 .and. VV(2) .eq. 1) Then
         OV_type = .True.
      Elseif (VV(1) .eq. 1 .and. OO(2) .eq. 1) Then
         VO_type = .True.
      Endif

      AA = Itype .eq. Moaindex .and. Atype .eq. Moaindex 
      BB = Itype .eq. Mobindex .and. Atype .eq. Mobindex 
      AB = Itype .eq. Moaindex .and. Atype .eq. Mobindex 
      BA = Itype .eq. Mobindex .and. Atype .eq. Moaindex 

      do i = i1, i2
      do a = a1,a2
         If (OO_type) Then
            If (AA) x(a,i) =  x(a,i)*Occ_numsa(a)*Occ_numsa(i)
            If (BB) x(a,i) =  x(a,i)*Occ_numsb(a)*Occ_numsb(i)
            If (AB) x(a,i) =  x(a,i)*Occ_numsa(a)*Occ_numsb(i)
            If (BA) x(a,i) =  x(a,i)*Occ_numsb(a)*Occ_numsa(i)
         Elseif (VV_type) Then
            If  (AA) x(a,i) =  x(a,i)*(Done-Occ_numsa(a))*
     +                                (Done-Occ_numsa(i))
            If  (BB) x(a,i) =  x(a,i)*(Done-Occ_numsb(a))*
     +                                (Done-Occ_numsb(i))
            If  (AB) x(a,i) =  x(a,i)*(Done-Occ_numsa(a))*
     +                                (Done-Occ_numsb(i))
            If  (BA) x(a,i) =  x(a,i)*(Done-Occ_numsb(a))*
     +                                (Done-Occ_numsa(i))
         Elseif (VO_type) Then
            If (AA) x(a,i) =  x(a,i)*(Done-Occ_numsa(a))*Occ_numsa(i)
            If (BB) x(a,i) =  x(a,i)*(Done-Occ_numsb(a))*Occ_numsb(i)
            If (AB) x(a,i) =  x(a,i)*(Done-Occ_numsa(a))*Occ_numsb(i)
            If (BA) x(a,i) =  x(a,i)*(Done-Occ_numsb(a))*Occ_numsa(i)
         Elseif (OV_type) Then
            If (AA) x(a,i) =  x(a,i)*Occ_numsa(a)*(Done-Occ_numsa(i))
            If (BB) x(a,i) =  x(a,i)*Occ_numsb(a)*(Done-Occ_numsb(i))
            If (AB) x(a,i) =  x(a,i)*Occ_numsa(a)*(Done-Occ_numsb(i))
            If (BA) x(a,i) =  x(a,i)*Occ_numsb(a)*(Done-Occ_numsa(i))
         Endif 
      enddo
      enddo

      return
      end



