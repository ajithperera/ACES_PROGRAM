
      subroutine Xgeev(
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
         Write(6,"(a)") "@-XGEEV: Both arrays must be static."
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
       n2 = index_table(c_index_size, ind_in(2))  ! pick up length of index

       m1 = index_table(c_index_size, ind_out(1))  ! pick up length of index
       m2 = index_table(c_index_size, ind_out(2))  ! pick up length of index

c---------------------------------------------------------------------------
c   Look up each array's address.
c---------------------------------------------------------------------------

      iarray_in  = get_index_from_base(address_table(array_in), x, 2)
      iarray_out = get_index_from_base(address_table(array_out), x, 2)

      Call diagonalize(x(iarray_in),1,n1,1,n2,x(iarray_out),1,m1,1,m2)

      return
      end
C
      Subroutine diagonalize(Hbar,i1,i2,j1,j2,Guess_vec,a1,a2,b1,b2)

      implicit none 
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'aces3_vars.h'

      Integer i,j,k,l,i1,i2,j1,j2,a1,a2,b1,b2,Ndim,Nroots 
      Integer Istart,Maxdim, Ierror 

      Parameter(Maxdim = 100)

      Double precision Hbar(i1:i2,j1:j2),Guess_vec(a1:a2,b1:b2)
      Double precision EvecR(Maxdim*Maxdim),EvecL(Maxdim*Maxdim)
      Double Precision EigR(Maxdim),EigI(Maxdim)
      Double Precision Work(4*Maxdim)

C Hbar is non-symmetric matrix, and we obtain both left/right eigenvectors.

      Ndim = (i2 -i1) + 1
C#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a,i4)") "@-Entry to XGEEV order of the matrix", 
     +                Ndim 
      k = 0
      do i = 1, ndim
         if (Hbar(i,i) .ne. 0.0D0) k=k+1
      enddo 
      call output(Hbar,1,k,1,k,Ndim,Ndim,1)
      Endif 
C#endif 

      If (Ndim .Gt. Maxdim) Then
         Write(6,"(a,a)") "@-XGEEV: Exceeded the max",
     +                    " dimension allowed in XGEEV"
         Call abort_job()
      Endif 
         
      Call Dgeev("V","V",Ndim, Hbar, Ndim, EigR, EigI,EvecL,NDIM,
     +            EvecR,Ndim,Work,4*Ndim,Ierror)

#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a)") "The eigenvectors of the incoming matrix"
      Write(6,"(a)") "The right-hand side"
      call output(EvecR,1,k,1,k,Ndim,Ndim,1)
C      Write(6,"(a)") "The left-hand side"
C      call output(EvecL,1,k,1,k,Ndim,Ndim,1)
      write(6,"(a,2i3)") "Lhs,rhs: ",lhs,rhs 
      endif
#endif 
C Copy the eigenvectors to returned array (Guess_Vec) and replace
C the digonal of incoming matrix (Hbar) with eigenvalues. 
C (by defalut copy the right eigenvector). If RHS or LHS is set
C then we choose. 

      If (RHS .EQ. 1 .AND. LHS .EQ. 0) Then
         Call Dcopy(Ndim*Ndim, EvecR, 1, Guess_vec,1)
      Elseif (RHS .EQ. 0 .AND. LHS .EQ. 1) Then
         Call Dcopy(Ndim*Ndim, EvecL, 1, Guess_vec,1)
      ELse if (RHS .EQ. 0 .AND. LHS .EQ. 0) Then
         Call Dcopy(Ndim*Ndim, EvecR, 1, Guess_vec,1)
      Endif 

      Call Zero(Hbar,Ndim*Ndim)
      Do i = 1, Ndim
         Hbar(i,i) = Eigr(i)
      Enddo 

C#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a)") "Eigenvalues/vectors of the incoming matrix"
      Write(6,"(a)") "The right-hand side"
      Write(6,"(6(1x,F10.5))") (EigR(i),i=1,k)
      Write(6,"(a)") "The imaginary eigenvalues"
      Write(6,"(6(1x,F10.5))") (EigI(i),i=1,k)
      call output(Guess_vec,1,k,1,k,Ndim,Ndim,1)
C      Write(6,"(a)") "The left-hand side"
C      call output(EvecL,1,k,1,k,Ndim,Ndim,1)
      endif 
C#endif 

      Return
      End

