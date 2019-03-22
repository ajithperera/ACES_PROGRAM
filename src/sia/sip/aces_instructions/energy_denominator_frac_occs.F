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
      subroutine energy_denominator_frac_occs(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c   Divides each sample of the block given by the array argument by the 
c   MP2 denominator
c       eps = epsilon(i)+epsilon(j)-epsilon(a)-epsilon(b)
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'mpif.h'
      include 'epsilon.h'
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
      integer array, index, nindex, ierr
      integer block, blkndx, seg
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

      if (nindex .eq. 4) then
         call eps_divide4_frac_occ(x(indblk), val1(1),val2(1),type(1),
     *                             val1(2),val2(2),type(2),
     *                             val1(3),val2(3),type(3),
     *                             val1(4),val2(4),type(4), 
     *                             epsilon,epsilonb)
      else if (nindex .eq. 6) then
         call eps_divide6_frac_occ(x(indblk), val1(1),val2(1),type(1),
     *                             val1(2),val2(2),type(2),
     *                             val1(3),val2(3),type(3),
     *                             val1(4),val2(4),type(4),
     *                             val1(5),val2(5),type(5),
     *                             val1(6),val2(6),type(6),
     *                             epsilon,epsilonb)
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

         call eps_divide2_frac_occ(x(indblk), val1(1),val2(1),type(1),
     *                             na1, na2,val1(2),val2(2),type(2),
     *                             ni1,ni2,epsilon,epsilonb)    
      else
         write(6,"(a,a)") " Error in apply_frac_occs: target array",
     *                    " must have either 4 or 2 indices."
         Write(6,"(a,i2,a)") " array has", nindex, "indices"

         call abort_job()
      endif
      
      return
      end

      subroutine eps_divide4_frac_occ(x, a1,a2,atype,
     *                                i1,i2,itype,b1,b2,btype,
     *                                j1,j2, jtype, epsilon, epsilonb)
      implicit none
      include 'interpreter.h'
      include 'frac_occ_vars.h'
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      integer a1,a2,b1,b2,i1,i2,j1,j2
      integer atype, itype, btype, jtype
      double precision x(a1:a2,i1:i2,b1:b2,j1:j2)
      integer voff_a,voff_b,vrt_aa,vrt_ab,vrt_ba,vrt_bb
      double precision epsilon(*), epsilonb(*)

      integer a,b,i,j
      double precision eps, epsa, epsb, epsi, epsj
      double precision val, Done 

      Data Done /1.0D0/

#ifdef _DEBUG_LVL0
      if (me .eq.0) Then
          Write(6,"(a)") "@-Entry to eps_divide4_frac_occ"
          Write(6,"(a)") "The Occ_num arrays"
          Write(6,"(6F4.1)") (Occ_numsa(i),i=1, Nbasis)
          Write(6,"(6F4.1)") (Occ_numsb(i),i=1, Nbasis)
      endif
#endif 

      voff_a = Nalpha_occupied
      voff_b = Nbeta_occupied

      do j = j1,j2
         if (jtype .eq. mobindex) then
            epsj = epsilonb(j)*OCC_numsb(j)
         else
            epsj = epsilon(j)*Occ_numsa(j)
         endif
      do b = b1,b2
         if (btype .eq. mobindex) then
            vrt_bb = b
            epsb = epsilonb(b)*(Done-Occ_numsb(vrt_bb))
         else
            vrt_ba = b
            epsb = epsilon(b)*(Done-Occ_numsa(vrt_ba))
         endif
      do i = i1, i2
         if (itype .eq. mobindex) then
            epsi = epsilonb(i)*Occ_numsb(i)
         else
            epsi = epsilon(i)*Occ_numsa(i)
         endif
      do a = a1,a2
         if (atype .eq. mobindex) then
            vrt_ab = a
            epsa = epsilonb(a)*(Done-Occ_numsb(vrt_ab))
         else
            vrt_aa = a
            epsa = epsilon(a)*(Done-Occ_numsa(vrt_aa))
         endif
         val =  x(a,i,b,j)
         eps = epsi + epsj - epsa - epsb
         if (eps .ne. 0.0) then
             x(a,i,b,j) = val/eps
         else 
             x(a,i,b,j) = Done
         endif 
      enddo
      enddo
      enddo
      enddo

      return
      end

      subroutine eps_divide6_frac_occ(x, a1,a2,atype,
     *                                i1,i2,itype,b1,b2,btype,
     *                                j1,j2, jtype, c1, c2, ctype, 
     *                                k1, k2, ktype, epsilon, epsilonb)
      implicit none
      include 'interpreter.h'
      include 'frac_occ_vars.h'
      include 'int_gen_parms.h'

      integer a1,a2,b1,b2,c1,c2,i1,i2,j1,j2,k1,k2
      integer atype, itype, btype, jtype, ctype, ktype
      integer voff_a,voff_b,vrt_aa,vrt_ab,vrt_ba,vrt_bb,vrt_ca
      integer vrt_cb
      double precision x(a1:a2,i1:i2,b1:b2,j1:j2,c1:c2,k1:k2)

      double precision epsilon(*), epsilonb(*)

      integer a,b,c,i,j,k
      double precision eps, epsa, epsb, epsc, epsi, epsj, epsk
      double precision val,Done

      Data Done /1.0D0/

c--------------------------------------------------------------------
c Set index types if necassary and perform simple check. 
c--------------------------------------------------------------------

      if (atype .eq. simple_index) atype = itype 
      if (itype .eq. simple_index) itype = atype 

      if (btype .eq. simple_index) btype = jtype 
      if (jtype .eq. simple_index) jtype = btype 

      if (ctype .eq. simple_index) ctype = ktype 
      if (ktype .eq. simple_index) ktype = ctype 

      voff_a = Nalpha_occupied
      voff_b = Nbeta_occupied

      do k = k1,k2
         if (ktype .eq. mobindex) then
            epsk = epsilonb(k)*Occ_numsb(k)
         else
            epsk = epsilon(k)*Occ_numsa(k)
         endif
      do c = c1,c2
         if (ctype .eq. mobindex) then
            vrt_cb = c
            epsc = epsilonb(c)*(Done-Occ_numsb(vrt_cb))
         else
            vrt_ca = c
            epsc = epsilon(c)*(Done-Occ_numsa(vrt_ca))
         endif
      do j = j1,j2
         if (jtype .eq. mobindex) then
            epsj = epsilonb(j)*Occ_numsb(j)
         else
            epsj = epsilon(j)*Occ_numsa(j)
         endif
      do b = b1,b2
         if (btype .eq. mobindex) then
            vrt_bb = b
            epsb = epsilonb(b)*(Done-Occ_numsb(vrt_bb))
         else
            vrt_ba = b
            epsb = epsilon(b)*(Done-Occ_numsa(vrt_ba))
         endif
      do i = i1, i2
         if (itype .eq. mobindex) then
            epsi = epsilonb(i)*Occ_numsb(i)
         else
            epsi = epsilon(i)*Occ_numsa(i)
         endif
      do a = a1,a2
         if (atype .eq. mobindex) then
            vrt_ab = a
            epsa = epsilonb(a)*(Done-Occ_numsb(vrt_ab))
         else
            vrt_aa = a
            epsa = epsilon(a)*(Done-Occ_numsa(vrt_aa))
         endif
         val =  x(a,i,b,j,c,k)
         eps = epsi + epsj + epsk - epsa - epsb - epsc

         if (eps .NE. 0.0D0) then
             x(a,i,b,j,c,k) = val/eps
         else 
             x(a,i,b,j,c,k) = Done
         endif 
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      return
      end

      subroutine eps_divide2_frac_occ(x, a1,a2,atype,na1,na2,i1,i2,
     *                                itype,ni1,ni2,epsilon,epsilonb)
      implicit none
      include 'interpreter.h'
      include 'frac_occ_vars.h'
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      integer a1,a2,i1,i2
      integer atype, itype
      integer na1, na2, ni1, ni2,voff_a, Voff_b
      integer vrt_ab, vrt_aa
      double precision x(na1:na2,ni1:ni2)

      double precision epsilon(*), epsilonb(*)

      integer a,i
      double precision eps, epsa, epsi
      double precision val, Done
      Data Done /1.0D0/

#ifdef _DEBUG_LVL0
      if (me .eq.0) then
       write(6,"(a,3(1x,i4))") "itype,atype,moabindex: ",itype,
     +                          atype, mobindex 
      endif 
#endif 

      voff_a = Nalpha_occupied 
      voff_b = Nbeta_occupied 

      do i = i1, i2
         if (itype .eq. mobindex) then
            epsi = epsilonb(i)*Occ_numsb(i)
         else
            epsi = epsilon(i)*Occ_numsa(i)
         endif
      do a = a1,a2
         if (atype .eq. mobindex) then
            vrt_ab = a
            epsa = epsilonb(a)*(Done-Occ_numsb(vrt_ab))
         else
             vrt_aa =  a
            epsa = epsilon(a)*(Done-Occ_numsa(vrt_aa))
         endif
         val =  x(a,i)
         eps = epsi - epsa

#ifdef _DEBUG_LVL0
         write(6,"(a,1x,i2,1x,i2,1x,i2,1x,i2,3F10.5)") 
     +                  "a,i,epsa,epsi,epsi-epsa",a,i,vrt_ab,
     +                   vrt_aa, epsa,epsi,eps
#endif 
         if (eps .NE. 0.0D0) then 
            x(a,i) = val/eps
         else
            x(a,i) = Done
         endif

C         if (dabs(eps) .lt. 1.0d-10) then
C
C            if (a .ne. i) write(6,*) ' Small denominator being
C     *         eliminated in energy_denominator',
C     *         'a, i', a,i, 'den' ,eps, 'value', x(a,i)
C         endif

      enddo
      enddo

      return
      end
