c---------------------------------------------------------------|
      subroutine sip_test(  array_table,
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
      logical test1,test2,test3,test4,test5,test6
      character*8 BLANK

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
      if (me .eq.0) Then
         Write(6,"(a)") "------Entering SIP_TEST SIP--------"
      endif 

      array = op(c_result_array)
      array_type = array_table(c_array_type, array)
      sind1 =  array_table(c_scalar_index, array)

      array = op(c_op1_array)
      array_type = array_table(c_array_type, array)
      sind2 =  array_table(c_scalar_index, array)

      qq=scalar_table(sind1)
      g=scalar_table(sind2)
 
      if (me .eq.0) Then
         Write(6,"(a,i4,1x,i4)") "The input scalars:", int(qq), 
     *                            int(g)
      endif 
       
      scalar_table(sind2) = 5

      return
      end

