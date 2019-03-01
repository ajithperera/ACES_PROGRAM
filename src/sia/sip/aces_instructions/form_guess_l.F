
      subroutine Form_guess_l(
     *                      array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)

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
      integer nsend,Ndim,Mdim

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
         Write(6,"(a)") "@-FORM_GUESS_L: Both arrays must be static."
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

      Ndim = n1
      Mdim = n2
      If (Ndim .Ne. Mdim) then
         Write(6,"(a)") "@-FORM_GUESS: Must be a square matrix"
         call abort_job()
      endif

      Call Get_guess_l(x(iarray_in),1,n1,1,n2,x(iarray_out),1,m1,1,m2,
     &                 Ndim)

      return
      end
C
      Subroutine Get_guess_l(Hbar,i1,i2,j1,j2,Guess_vec,a1,a2,b1,b2,
     &                       Ndim)

      implicit none 
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      Integer i,j,k,l,i1,i2,j1,j2,a1,a2,b1,b2,Ndim,Nroots 
      Integer Istart,Maxdim, Ierror 

      Double precision Hbar(i1:i2,j1:j2),Guess_vec(a1:a2,b1:b2)
      Double precision EvecR(Ndim*Ndim),EvecL(Ndim*Ndim)
      Double Precision EigR(Ndim),EigI(Ndim)
      Double Precision Work(4*Ndim)

C Hbar is non-symmetric matrix, and we obtain both left/right eigenvectors.

#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a,i4)") "@-Entry to Get_guess order of the matrix", 
     +                Ndim 
#endif 

      Call Dgeev("V","V",Ndim, Hbar, Ndim, EigR, EigI,EvecL,NDIM,
     +            EvecR,Ndim,Work,4*Ndim,Ierror)

#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a)") "H(occ,occ) eigenvalues/vectors"
      Write(6,"(a)") "The right-hand side"
      Write(6,"(6(1x,F10.5))") (EigR(i),i=1,Ndim)
      Write(6,"(6(1x,F10.5))") (EigI(i),i=1,Ndim)
      call output(EvecR,1,Ndim,1,Ndim,Ndim,Ndim,1)
      Write(6,"(a)") "The left-hand side"
      call output(EvecL,1,Ndim,1,Ndim,Ndim,Ndim,1)
      endif 
#endif 
     
      Call Process_vecs_l(Eigr,EvecR,EvecL,Guess_vec,Ndim,
     +                  a1,a2,b1,b2)

#ifdef _DEBUG_LVL0
      if (me.eq.0) then
      Nroots = eom_roots
      Write(6,"(a)") "Initial Right guess vectors for EOM-IP/EA"
      call output(Guess_vec,1,Ndim,1,Nroots,Ndim,Ndim,1)
      Write(6,"(a)") "Initial Left, guess vectors for EOM-IP/EA"
      call output(Guess_vec(1,51),1,Ndim,1,Nroots,Ndim,Ndim,1)
      endif 
#endif 
C Return the eigenvalues in the diagonal. 

      Call zero(Hbar,Ndim*Ndim)
      do i =1, Ndim
         Hbar(i,i) = Eigr(i)
      enddo 

      Return
      End

      Subroutine Process_vecs_l(Eigr,EvecR,EvecL,Guess_vec,Ndim,
     +                        a1,a2,b1,b2)

      Implicit None
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      Integer a1,a2,b1,b2
      Integer i,j,Ndim,k,l,Istart,Nroots 

      Double Precision Eigr(Ndim), EvecR(Ndim,Ndim)
      Double Precision EvecL(Ndim,Ndim)
      Double Precision Guess_vec(a1:a2,b1:b2)
      Double Precision Norm(Ndim)
      Double Precision Dtmp,DtmpR,DtmpL,Ddot

C Reorder the eigenvalues and vectors descending order

      If (Ndim .Eq. 1) Return

      Do i = 1, Ndim-1
         Do j = i+1, Ndim 
         If (EigR(i) .Gt. EigR(j)) Then
              Dtmp    = EigR(i)
              EigR(i) = EigR(j)
              EigR(j) = Dtmp
              Do k = 1, Ndim
                 DtmpR      = EvecR(k,i) 
                 EvecR(k,i) = EvecR(k,j)
                 EvecR(k,j) = DtmpR
       
                 DtmpL      = EvecL(k,i) 
                 EvecL(k,i) = EvecL(k,j)
                 EvecL(k,j) = DtmpL
              Enddo 
          Endif 
          Enddo       
      Enddo       

C#ifdef _DEBUG_LVL0
      if (me.eq.0) then
      Write(6,"(a)") "The reordered eigenvalues and vectors (left)"
      Write(6,"(6(1x,F10.5))") (EigR(i),i=1,Ndim)
      Write(6,"(a)") "The left-hand side"
      call output(EvecL,1,Ndim,1,Ndim,Ndim,Ndim,1)
C      Write(6,"(a)") "The left-hand side"
C      call output(EvecL,1,NDIM,1,NDim,Ndim,Ndim,1)
      Endif 
C#endif 
C 
C Normalized the Eigenvectors so that they are bi-orthonormal.
C
      Do i = 1, Ndim
         Norm(i) = Ddot(Ndim, EvecR(1,i), 1, EvecL(1,i),1)
         If (Dsqrt(Norm(i)) .lt. 1.0D-9) Then 
             Write(6,"(a)") "@-Process_vecs: Division by a small",
     +                      " number. Discontinuing."
             Call abort_job()
         Else
             Norm(i) = 1.0D0/Dsqrt(Norm(i))
         Endif 
      Enddo 
C
#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a)") "<L|R>"
      Write(6,"(6(1x,F10.5))") (Norm(i),i=1,Ndim)
      Endif 
#endif 
C
      Do i = 1, Ndim
         Call Dscal(Ndim, Norm(i), EvecR(1,i), 1)
         Call Dscal(Ndim, Norm(i), EvecL(1,i), 1)
      Enddo
C
C Form the initial Guess in appropriate array. We only use the first
C nroot (number of roots from the int_gen_param common block). The first
C 25 positions are reseved for R vectors and next 25 position is reserved
C for L vectors. 
C
      Nroots = eom_roots 
#ifdef _DEBUG_LVL0
      if (me.eq.0) then
      Write(6,"(a,i4)") "The number of roots:", nroots 
      endif 
#endif 
      If (Nroots .EQ. 0) Then
         Write(6,"(a,a)") "@-PROCESS_VECS_L: The requested number of"
     +                    " roots can not be zero"
         Call abort_job()
      Endif 

      If (Nroots .Gt. 50 .or. Nroots .Gt. Ndim) Then
         Write(6,"(a,a,a)") "@-PROCESS_VECS_L: The requested number of",
     +                  "roots are greater than allowed maximum",
     +                  " or allowed."
         Call abort_job()
      Endif 

      Do j = 1, Nroots 
         Call Dcopy(Ndim, EvecL(1,j), 1, Guess_vec(1,j), 1)
      Enddo 

      Return
      End
