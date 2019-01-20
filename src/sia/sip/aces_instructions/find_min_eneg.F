
      subroutine Find_min_eneg(
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
      integer array_in,array_out, index, nindex_in,nindex_out, ierr
      integer block, blkndx, seg
      integer find_current_block
      integer array_in_type, array_out_type
      integer*8 indblk1, indblk2, get_block_index
      integer stack,ind_in(mx_array_index),ind_out(mx_array_index)
      integer nsend

      integer comm

      integer fop1(mx_array_index), fop2(mx_array_index)
      integer sop1(mx_array_index), sop2(mx_array_index)
      integer sindex(mx_array_index), findex(mx_array_index)
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
      array_in = op(c_result_array)
      array_out = op(c_op1_array)

      array_in_type   = array_table(c_array_type, array_in)
      array_out_type  = array_table(c_array_type, array_out)
C
      if (array_in_type .ne. static_array .or.array_out_type 
     &                  .ne. static_array) then
         Write(6,"(a)") "@-FIND_MIN_ENEG: Both arrays must be static."
         call abort_job()
      endif

      nindex_in  = array_table(c_nindex, array_in)
      nindex_out = array_table(c_nindex,array_out)

      do i = 1, nindex_in
         ind_in(i) = array_table(c_index_array1+i-1,array_in)
      enddo

      do j = 1, nindex_out
         ind_out(j) = array_table(c_index_array1+j-1,array_out)
      enddo

       n1 = index_table(c_index_size, ind_in(1))  ! pick up length of index

       m1 = index_table(c_index_size, ind_out(1))  ! pick up length of index
       m2 = index_table(c_index_size, ind_out(2))  ! pick up length of index

c---------------------------------------------------------------------------
c   Look up each array's address.
c---------------------------------------------------------------------------

      iarray_in  = get_index_from_base(address_table(array_in), x, 2)
      iarray_out = get_index_from_base(address_table(array_out), x, 2)

      Call What_2exclude(x(iarray_in),1,n1,x(iarray_out),1,m1,1,m2)

      return
      end
C
      Subroutine What_2exclude(Eigs_curnt,i1,i2,Roots,a1,a2,b1,b2)

      implicit none 
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'eom_truncate.h'

      Integer i,j,l,m,i1,i2,j1,j2,a1,a2,b1,b2
      Integer Ndim1,Ndim2,Ndim3,Iloc,Ilocmin
      Integer Istart,Ierror 
      Double precision Eigs_curnt(i1:i2),Roots(a1:a2,b1:b2)
      Double Precision X,X2,Vlocmin,Tol

      Parameter(Tol = 1.0D-03)

      Double Precision Tmp(Maxdim)

      Ndim1 = (i2 -i1) + 1
      Ndim2 = (a2 -a1) + 1
      Ndim3 = (b2 -b1) + 1

      If (Ndim1 .Gt. Maxdim) Then
         Write(6,"(a,a)") "@-FIND_MIN_ENEG: Exceeded the max",
     +                    " dimension allowed in FIND_MIN_ENEG"
         Call abort_job()
      Endif 
      
C#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      l=0
      m=0
      do i=1,ndim1
         if (Eigs_curnt(i) .ne. 0.0d0) l=l+1
      enddo
      do i=1,Ndim2
         if (roots(i,1) .ne. 0.0d0) m=m+1
      enddo
      Write(6,"(2a,2i4)") "@-Entry to FIND_MIN_ENEG actual order",
     +                   " of the arrays", l,m
      Write(6,"(a)") "The current eigenvalues"
      Write(6,"(6(1x,F10.5))") (Eigs_curnt(i),i=1,L)
      Write(6,"(a)") "The roots already found"
      Write(6,"(6(1x,F10.5))") ((Roots(i,j),i=1,M),j=1,Ndim3)
      Endif 
C#endif 

      Call Dcopy(Ndim1,Eigs_curnt,1,Tmp,1)

      Do I = 1, Ndim2
         Do J = 1, Ndim1
            Call Fndclose(Ndim1,Tmp,Roots(J,I),X,Iloc)
            X2 =Dabs(Roots(J,I)-Tmp(Iloc))
            If (X2 .LT. Tol) Tmp(Iloc) = 1.0D+30
         Enddo
      Enddo

      Vlocmin =  1.0D+30
      Do I = 1, Ndim1
         X = Tmp(I)
         IF (X .LT. Vlocmin) Then
            Ilocmin = I
            Vlocmin = Tmp(I)
         Endif
      Enddo

      Iselect = Ilocmin 

C#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a,a,I2,1X,F10.5)") "The location and the eigevalue",
     +                            " selected: ", Ilocmin, Vlocmin
      endif 
C#endif 
     
      Return
      End
