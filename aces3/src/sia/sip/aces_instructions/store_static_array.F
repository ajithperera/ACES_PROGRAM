
      subroutine Store_static_array(
     *                 array_table, narray_table,
     *                 index_table,
     *                 nindex_table, segment_table, nsegment_table,
     *                 block_map_table, nblock_map_table,
     *                 scalar_table, nscalar_table,
     *                 address_table, op)

c---------------------------------------------------------------|
      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'mpif.h'
      include 'parallel_info.h'
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

      integer i, j, k,n1,n2,m1,m2
      integer array_in,index, nindex_in,nindex_out, ierr
      integer block, blkndx, seg
      integer find_current_block
      integer array_in_type, array,array_type
      integer*8 indblk1, indblk2, get_block_index
      integer stack,ind_in(mx_array_index),ind_out(mx_array_index)
      integer nsend,sind,nindex 

      integer comm

      integer fop1(mx_array_index), fop2(mx_array_index)
      integer sop1(mx_array_index), sop2(mx_array_index)
      integer sindex(mx_array_index), findex(mx_array_index)
      integer type
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
      array_in = op(c_result_array)
      array_in_type   = array_table(c_array_type, array_in)
C
      if (array_in_type .ne. static_array) then
         Write(6,"(a)") "@-STORE_STATIC: Array must be static."
         call abort_job()
      endif

      nindex_in  = array_table(c_nindex, array_in)

      do i = 1, nindex_in
         ind_in(i) = array_table(c_index_array1+i-1,array_in)
      enddo

       n1 = index_table(c_index_size, ind_in(1))  ! pick up length of index
       n2 = index_table(c_index_size, ind_in(2))  ! pick up length of index

c   Look up each array's address.

      iarray_in  = get_index_from_base(address_table(array_in), x, 2)

      array       = op(c_op1_array)
      array_type  = array_table(c_array_type, array)
      nindex      = array_table(c_nindex, array)

      if (array_type .ne. scalar_value) then
         Write(6,"(a)") "@-STORE_STATIC: Must be an scalar."
         call abort_job()
      endif

      sind =  array_table(c_scalar_index, array)
      type = INT (scalar_table(sind))

      Call Store(x(iarray_in),1,n1,1,n2,Type)

      return
      end
C
      Subroutine Store(Data,i1,i2,j1,j2,Type)

      implicit none 
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      Integer i,j,k,l,i1,i2,j1,j2
      Integer Type,Iroot

      Double Precision Data(i1:i2,j1:j2)


      If (me .eq. 0) Then

         If (Type .EQ. 1) Then
            Do Iroot = 1, eom_roots 
               EE_CCSD(Iroot) = Data(Iroot,1)
            Enddo
         Elseif (Type .EQ. -1) Then
            Do Iroot = 1, eom_roots 
               EE_CCSD(30+Iroot) = Data(Iroot,1)
            Enddo
         Elseif (Type .EQ. 2) Then
            Do Iroot = 1, eom_roots 
               E4_AAA(Iroot) = Data(Iroot,1)
            Enddo
         Elseif (Type .EQ. -2) Then
            Do Iroot = 1, eom_roots 
               E4_AAA(30+Iroot) = Data(Iroot,1)
            Enddo
         Elseif (Type .EQ. 3) Then
            Do Iroot = 1, eom_roots 
               E5_AAA(Iroot) = Data(Iroot,1)
            Enddo
         Elseif (Type .EQ. -3) Then
            Do Iroot = 1, eom_roots 
               E5_AAA(30+Iroot) = Data(Iroot,1)
            Enddo
         Elseif (Type .EQ. 4) Then
            Do Iroot = 1, eom_roots 
               E4_AAB(Iroot) = Data(Iroot,1)
            Enddo
         Elseif (Type .EQ. -4) Then
            Do Iroot = 1, eom_roots 
               E4_AAB(30+Iroot) = Data(Iroot,1)
            Enddo
         Elseif (Type .EQ. 5) Then
            Do Iroot = 1, eom_roots 
               E5_AAB(Iroot) = Data(Iroot,1)
            Enddo
         Elseif (Type .EQ. -5) Then
            Do Iroot = 1, eom_roots 
               E5_AAB(30+Iroot) = Data(Iroot,1)
            Enddo
         Elseif (Type .EQ. -6) Then
            Do Iroot = 1, eom_roots 
               E4_BBB(30+Iroot) = Data(Iroot,1)
            Enddo
         Elseif (Type .EQ. -7) Then
            Do Iroot = 1, eom_roots 
               E5_BBB(30+Iroot) = Data(Iroot,1)
            Enddo
         Elseif (Type .EQ. -8) Then
            Do Iroot = 1, eom_roots 
               E4_BBA(30+Iroot) = Data(Iroot,1)
            Enddo
         Elseif (Type .EQ. -9) Then
            Do Iroot = 1, eom_roots 
               E5_BBA(30+Iroot) = Data(Iroot,1)
            Enddo
         Endif
      Endif 

CSSS#ifdef _DEBUG_LVL0
      if (me. eq. 0) then
         if (type .eq. 1) then
         Write(6,"(a,6F12.6)") "S-EOM-CCSD= ",
     *                         (ee_ccsd(i),i=1,eom_roots)
         endif 
         if (type .eq. -1) then
         Write(6,"(a,6F12.6)") "T-EOM-CCSD= ",
     *                         (ee_ccsd(30+i),i=1,eom_roots)
         endif 
         If (type .eq. 2) then 
         Write(6,"(a,6F12.8)") "EOM-CCSD(T)-E4-AAAS= ",
     *                          (e4_aaa(i),i=1,eom_roots)
         endif 
         If (type .eq. -2) then 
         Write(6,"(a,6F12.8)") "EOM-CCSD(T)-E4-AAAT= ",
     *                          (e4_aaa(30+i),i=1,eom_roots)
         endif 
         If (type .eq. 3) then 
         Write(6,"(a,6F12.8)") "EOM-CCSD(T)-E5-AAAS= ",
     *                          (e5_aaa(i),i=1,eom_roots)
         endif 
         If (type .eq. -3) then 
         Write(6,"(a,6F12.8)") "EOM-CCSD(T)-E5-AAAT= ",
     *                          (e5_aaa(30+i),i=1,eom_roots)
         endif 
         If (type .eq. 4) then 
         Write(6,"(a,6F12.8)") "EOM-CCSD(T)-E4-AABS= ",
     *                          (e4_aab(i),i=1,eom_roots)
         endif 
         If (type .eq. -4) then 
         Write(6,"(a,6F12.8)") "EOM-CCSD(T)-E4-AABT= ",
     *                          (e4_aab(30+i),i=1,eom_roots)
         endif 
         If (type .eq. 5) then 
         Write(6,"(a,6F12.8)") "EOM-CCSD(T)-E5-AABS= ",
     *                          (e5_aab(i),i=1,eom_roots)
         endif 
         If (type .eq. -5) then 
         Write(6,"(a,6F12.8)") "EOM-CCSD(T)-E5-AABT= ",
     *                          (e5_aab(30+i),i=1,eom_roots)
         endif 
         if (type .eq. -6) Then
         Write(6,"(a,6F12.8)") "EOM-CCSD(T)-E4-BBBT= ",
     *                          (e4_bbb(30+i),i=1,eom_roots)
         endif 
         if (type .eq. -7) Then
         Write(6,"(a,6F12.8)") "EOM-CCSD(T)-E5-BBBT= ",
     *                          (e5_bbb(30+i),i=1,eom_roots)
         endif 
         if (type .eq. -8) Then
         Write(6,"(a,6F12.8)") "EOM-CCSD(T)-E4-BBAT= ",
     *                          (e4_bba(30+i),i=1,eom_roots)
         endif 
         if (type .eq. -9) Then
         Write(6,"(a,6F12.8)") "EOM-CCSD(T)-E5-BBAT= ",
     *                          (e5_bba(30+i),i=1,eom_roots)
         endif 
      endif 
CSSS#endif 

      Return
      End
