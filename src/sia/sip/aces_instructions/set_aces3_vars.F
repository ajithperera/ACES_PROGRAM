c---------------------------------------------------------------|
      subroutine Set_aces3_vars(   array_table,
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)
C-------------------------------------------------------------------------
C Can be uses to setup variables in the common block aces3_vars from
C a SIAL code
c-------------------------------------------------------------------------- 
      implicit none
      include 'mpif.h'
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'pst_functions.h'
      include 'aces3_vars.h'
      include 'eom_truncate.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

c---------------Variables--------------------------------------------

      integer narray_table, nindex_table, nsegment_table,
     *        nblock_map_table, array
      integer op(loptable_entry),sind1,Sind2
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      double precision qq
      integer*8 address_table(narray_table), get_index_from_base
      integer*8 get_block_data_index
      integer find_current_block, block

      integer array_type, nindex, inda(2), indb(2)

      logical direct_flag

      integer i, j, k,qq1,is_it_worker(max_procsx)
      integer nind, nsend,ind1,w_worker(max_procsx)
      integer ierr, handle,our_workers
      integer blk, blkndx, maxblk
      integer*8 iarray, farray, indblk, get_block_index
      integer*8 integral_scr
      integer flag
      integer stack,g

      double precision x(1)
      double precision dbuf(1), temp
      integer ibuf(1)
      integer int_type
c----------------------------------------------------------------------
#ifdef ALTIX
      pointer (dptr, dbuf)
      pointer (xptr, x)
      pointer (iptr, ibuf)
#else
      common dbuf, x
      equivalence (ibuf(1), dbuf(1))
#endif

      array = op(c_result_array)

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
      xptr = dshptr
#endif
c-------------------------------------------------------------------

      array = op(c_result_array)
      array_type = array_table(c_array_type, array)
      sind1 =  array_table(c_scalar_index, array)

      qq=scalar_table(sind1)
      If (Int(qq) .EQ. 0) Then
         RHS = 0
         LHS = 0
      Elseif (Int(qq) .EQ. 1) Then
         RHS = qq
         LHS = 0
      Elseif (Int(qq) .EQ. 2) Then
         RHS = 0
         LHS = (Int(qq)-1)
      Endif
#ifdef _DEBUG_LVL0
      If (me .eq.0) Write(6,"(a,2i2)") "The right and left sides: ",
     +              RHS, LHS
#endif 

c------------Information about scalar--------------------------------------

      array = op(c_op1_array)
      array_type = array_table(c_array_type, array)
      sind2 =  array_table(c_scalar_index, array)

      g=scalar_table(sind2)

      If (Int(g) .EQ. 1) Then 
         AAAA = 1
         BBBB = 0
         ABAB = 0
      Elseif (Int(g) .EQ. 2) Then 
         AAAA = 0
         BBBB = 1
         ABAB = 0
      Elseif (Int(g) .EQ. 3) Then 
         AAAA = 0
         BBBB = 0
         ABAB = 1
      Endif 
#ifdef _DEBUG_LVL0
      If (me .eq.0) Write(6,"(a,3i2)") "The spin_types: " ,
     +              AAAA,BBBB,ABAB
#endif 
C--------------updating eom_truncate variables header file----

      If (Int(g) .EQ. 10) Then
         Ndim_new = qq
      Elseif (Int(g) .EQ. 11) Then
         Ndim_old = qq
      Elseif (Int(g) .EQ. 12) Then
         Omega = Dble(qq)
      Endif 
C--------------Returning vars in eomdata header file----

      If (Int(g) .EQ. -1) Then
         scalar_table(sind1) = Iselect 
      Endif 

      return
      end

