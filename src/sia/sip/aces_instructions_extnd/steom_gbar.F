c---------------------------------------------------------------|
c  This subroutine is a superinstruction for the SIAL code      |
c  It calculates Gbar matrix elements, needed for STEOM         |
c  calculations.                                                | 
c  Added by Denis Bokhan, Moscow Lomonosov State University,    |
c              Moscow, RUSSIAN FEDERATION                       |
c---------------------------------------------------------------|
      subroutine steom_gbar (array_table,
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)
c-------------------------------------------------------------------------- 
      implicit none
      include 'mpif.h'
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'pst_functions.h'
      include 'F12_gr.h'
      include 'daveom.h'
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
      integer sval1(mx_array_index),sval2(mx_array_index)
      integer tval1(mx_array_index),tval2(mx_array_index)
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

      nindex_source = array_table(c_nindex, source)
      nindex_target = array_table(c_nindex, target)

      if (nindex_source .ne. nindex_target) then
         print *,'Error: Both arrays in steom_gbar must use ',
     *           'the same number of indices.'
         print *,'First array has ',nindex_source,' indices.'
         print *,'Second array has ',nindex_target,' indices.'
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
     *                          nindex_table, sval1(i), sval2(i))

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
     *                          nindex_table, tval1(i), tval2(i))

      enddo

c-----------------Two-index arrays------------------------------------------      
       IF(nindex_source.eq.2) then
       call do_gbar_2ind(x(isource),x(itarget),
     *                      sval1(1),sval2(1),sval1(2),sval2(2),
     *                      tval1(1),tval2(1),tval1(2),tval2(2),
     *                      RL)
       ENDIF 
c-----------------Four-index arrays-----------------------------------------
       IF(nindex_source.eq.4) then
       call do_gbar_4ind(x(isource),x(itarget),
     *                      sval1(1),sval2(1),sval1(2),sval2(2),
     *                      sval1(3),sval2(3),sval1(4),sval2(4),
     *                      tval1(1),tval2(1),tval1(2),tval2(2),
     *                      tval1(3),tval2(3),tval1(4),tval2(4),
     *                      RL)
       ENDIF
c--------------------------------------------------------------------------- 
      return
      end

c################################################################################
c################################################################################

      subroutine do_gbar_2ind(W_pq,gbar,
     *                   a1,a2,b1,b2,
     *                   i1,i2,j1,j2,RL)
      implicit none
      include 'interpreter.h'
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      integer a1,a2,b1,b2,c1,c2,d1,d2
      integer i1,i2,j1,j2,k1,k2,l1,l2
      integer nvirta,nocca
      double precision W_pq(a1:a2,b1:b2)
      double precision gbar(i1:i2,j1:j2)
      integer i,j,k,l,RL,kk
      double precision val, div
      
c      if(me.eq.0) write(6,*) a1,a2,'--->>>',b1,b2
c      if(me.eq.0) write(6,*) i1,i2,'<<<>>>',j1,j2

      nocc a = Nalpha_occupied 
      nvirta = Nalpha_virtual 

c-----------------G_ik-------------------------------------------       
      IF(RL.eq.1) then

      DO i=i1,i2
      DO j=j1,j2
      DO k=a1,a2
      IF(j.eq.k) then
      gbar(i,j) = gbar(i,j) + W_pq(j,i)
      ENDIF
      ENDDO
      ENDDO
      ENDDO      

      ENDIF       
c-----------------G_ac------------------------------------------- 

      IF(RL.eq.2) then

      DO i=i1,i2
      DO j=j1,j2
      DO k=b1,b2
      kk = k + nocca
      IF(j.eq.kk) then
      gbar(i,j) = gbar(i,j) + W_pq(i,k)
      ENDIF
      ENDDO
      ENDDO
      ENDDO

      ENDIF

c-----------------G_ia-------------------------------------------
      IF(RL.eq.3) then

      DO j=j1,j2
      DO i=i1,i2
      DO k=a1,a2
      IF(i.eq.k) then
      gbar(i,j) = gbar(i,j) + W_pq(i,j)
      ENDIF
      ENDDO
      ENDDO
      ENDDO

      ENDIF


      IF(RL.eq.4) then

      DO i=i1,i2
      DO j=j1,j2
      DO k=b1,b2
      kk = k + nocca
      IF(j.eq.kk) then
      gbar(i,j) = gbar(i,j) + W_pq(i,k)
      ENDIF
      ENDDO
      ENDDO
      ENDDO

      ENDIF
c------------------------------------------------------------------


      return
      end

c################################################################################
c################################################################################


       subroutine do_gbar_4ind(W_pq,gbar,
     *                   a1,a2,b1,b2,c1,c2,d1,d2,
     *                   i1,i2,j1,j2,k1,k2,l1,l2,RL)
      implicit none
      include 'interpreter.h'
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      integer a1,a2,b1,b2,c1,c2,d1,d2
      integer i1,i2,j1,j2,k1,k2,l1,l2
      integer nocca,nvirta 
      double precision W_pq(a1:a2,b1:b2,c1:c2,d1:d2)
      double precision gbar(i1:i2,j1:j2,k1:k2,l1:l2)

      integer i,j,k,l,RL,kk,ll,a,b,c,d
      double precision val, div

      nocca  = Nalpha_occupied 
      nvirta = Nalpha_virtual 

c      if(me.eq.0) write(6,*) a1,a2,'--->>>',b1,b2
c      if(me.eq.0) write(6,*) i1,i2,'<<<>>>',j1,j2

c       call igetrec(1, 'JOBARC', 'OCCUPYA0', 1, nocca)

c---------------G_ik_ld-----------------------------------------|

      IF(RL.eq.5) then
      DO l=l1,l2
      DO k=k1,k2
      DO i=i1,i2
      DO j=j1,j2
      DO a=a1,a2
      IF(j.eq.a) then
      gbar(i,j,k,l) = gbar(i,j,k,l) + W_pq(j,i,k,l)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      ENDIF 


      IF(RL.eq.6) then
      DO i=i1,i2
      DO j=j1,j2
      DO l=l1,l2
      DO k=k1,k2
      DO a=a1,a2
      IF(k.eq.a) then
      gbar(i,j,k,l) = gbar(i,j,k,l) + W_pq(k,l,i,j)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO  

      ENDIF


      IF(RL.eq.7) then
      DO i=i1,i2
      DO j=j1,j2
      DO k=k1,k2
      DO l=l1,l2
      DO a=a1,a2
      kk = a + nocca
      IF(l.eq.kk) then
      gbar(i,j,k,l) = gbar(i,j,k,l) + W_pq(a,k,i,j)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      ENDIF
c---------------G_ac_ld-----------------------------------------|

      IF(RL.eq.8) then
      DO i=i1,i2
      DO k=k1,k2
      DO l=l1,l2
      DO j=j1,j2
      DO a=a1,a2
      kk = a + nocca
      IF(j.eq.kk) then
      gbar(i,j,k,l) = gbar(i,j,k,l) + W_pq(a,i,k,l)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO


      ENDIF


      IF(RL.eq.9) then
      DO i=i1,i2
      DO j=j1,j2
      DO k=k1,k2
      DO l=l1,l2
      DO a=a1,a2
      kk = a + nocca
      IF(l.eq.kk) then
      gbar(i,j,k,l) = gbar(i,j,k,l) + W_pq(a,k,i,j)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO


      ENDIF


      IF(RL.eq.10) then
      DO i=i1,i2
      DO j=j1,j2
      DO l=l1,l2
      DO k=k1,k2
      DO a=a1,a2
      IF(k.eq.a) then
      gbar(i,j,k,l) = gbar(i,j,k,l) + W_pq(k,l,i,j)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO


      ENDIF

c---------------G_ik_jl-----------------------------------------|

      IF(RL.eq.11) then
      DO j=j1,j2
      DO l=l1,l2
      DO k=k1,k2
      DO i=i1,i2
      DO a=a1,a2
      IF(i.eq.a) then
      gbar(i,j,k,l) = gbar(i,j,k,l) + W_pq(i,j,l,k)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO


      ENDIF

      IF(RL.eq.12) then
      DO i=i1,i2
      DO j=j1,j2
      DO l=l1,l2
      DO k=k1,k2
      DO a=a1,a2
      IF(k.eq.a) then
      gbar(i,j,k,l) = gbar(i,j,k,l) + W_pq(k,l,j,i)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      ENDIF

      IF(RL.eq.13) then
      DO j=j1,j2
      DO l=l1,l2
      DO i=i1,i2
      DO k=k1,k2
      DO a=a1,a2
      DO c=c1,c2
      IF(i.eq.a) then
      IF(k.eq.c) then
      gbar(i,j,k,l) = gbar(i,j,k,l) + W_pq(i,j,k,l)
      ENDIF
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      ENDIF
c---------------G_ac_bd-----------------------------------------|      

      IF(RL.eq.14) then
      DO i=i1,i2
      DO k=k1,k2
      DO l=l1,l2
      DO j=j1,j2
      DO a=a1,a2
      kk = a + nocca
      IF(j.eq.kk) then
      gbar(i,j,k,l) = gbar(i,j,k,l) + W_pq(a,i,k,l)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      

      ENDIF


      IF(RL.eq.15) then
      DO i=i1,i2
      DO j=j1,j2
      DO k=k1,k2
      DO l=l1,l2
      DO a=a1,a2
      kk = a + nocca
      IF(l.eq.kk) then
      gbar(i,j,k,l) = gbar(i,j,k,l) + W_pq(a,k,i,j)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO


      ENDIF


      IF(RL.eq.16) then
      DO j=j1,j2
      DO l=l1,l2
      DO i=i1,i2
      DO k=k1,k2
      DO a=a1,a2
      DO c=c1,c2
      kk = nocca + a
      ll = nocca + c
      IF(i.eq.kk) then
      IF(k.eq.ll) then
      gbar(j,i,l,k) = gbar(j,i,l,k) + W_pq(a,j,c,l)
      ENDIF
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO


      ENDIF
c---------------G_ab_ij-----------------------------------------| 


      IF(RL.eq.17) then
      DO i=i1,i2
      DO k=k1,k2
      DO l=l1,l2
      DO j=j1,j2
      DO d=d1,d2
      kk = d + nocca
      IF(j.eq.kk) then
      gbar(i,j,k,l) = gbar(i,j,k,l) + W_pq(k,l,i,d)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO


      ENDIF


      IF(RL.eq.18) then
      DO i=i1,i2
      DO j=j1,j2
      DO k=k1,k2
      DO l=l1,l2
      DO a=a1,a2
      IF(l.eq.a) then
      gbar(i,j,k,l) = gbar(i,j,k,l) + W_pq(l,k,i,j)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      ENDIF


      IF(RL.eq.19) then
      DO i=i1,i2
      DO k=k1,k2
      DO j=j1,j2
      DO l=l1,l2
      DO a=a1,a2
      DO d=d1,d2
      ll = nocca + d
      IF(l.eq.a) then
      IF(j.eq.ll) then
      gbar(i,j,k,l) = gbar(i,j,k,l) + W_pq(l,k,i,d)
      ENDIF
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO


      ENDIF      
c---------------G_ab_ij-----------------------------------------|

      IF(RL.eq.20) then
      DO k=k1,k2
      DO l=l1,l2
      DO j=j1,j2
      DO a=a1,a2
      kk = a + nocca
      DO i=i1,i2
      IF(i.eq.kk) then
      gbar(i,j,k,l) = gbar(i,j,k,l) + W_pq(a,j,l,k)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO


      ENDIF

      IF(RL.eq.21) then
      DO k=k1,k2
      DO l=l1,l2
      DO i=i1,i2
      DO a=a1,a2
      DO j=j1,j2
      IF(j.eq.a) then
      gbar(i,j,k,l) = gbar(i,j,k,l) + W_pq(j,i,l,k)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO


      ENDIF

      IF(RL.eq.22) then
      DO k=k1,k2
      DO l=l1,l2
      DO i=i1,i2
      DO j=j1,j2
      DO a=a1,a2
      DO b=b1,b2
      ll = nocca + b
      IF(j.eq.a) then
      IF(i.eq.ll) then
C      if(me.eq.0) write(6,*) i,j,k,l,'||||====',kk
      gbar(i,j,k,l) = gbar(i,j,k,l) + W_pq(j,b,l,k)
      ENDIF
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO


      ENDIF



      return
      end

