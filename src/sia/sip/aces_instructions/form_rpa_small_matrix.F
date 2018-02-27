
      subroutine Form_rpa_small_matrix(
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
      include 'aces3_vars.h'
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
      integer*8 iarray_in1,iarray_in2,ievec,daddr,iarray_ca

      integer i, j, k,n1,n2,m1,m2,l1,l2,m
      integer array_in,array_out,index,nindex_in,nindex_out, ierr
      integer nindex_ca,array_ca
      integer block, blkndx, seg
      integer find_current_block
      integer array_in_type, array_out_type,array_ca_type
      integer*8 indblk1, indblk2, get_block_index
      integer stack,ind_in(mx_array_index),ind_out(mx_array_index)
      integer ind_ca(mx_array_index)
      integer nsend, Ndim

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
         Write(6,"(a,a)") "FORM_RPA_SMALL_MATRIX: Both arrays must ",
     &                     "be static."
         call abort_job()
      endif

      nindex_in  = array_table(c_nindex,array_in)
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

      iarray_in1 = get_index_from_base(address_table(array_in), x, 2)
      iarray_in2 = get_index_from_base(address_table(array_out), x, 2)

      Ndim = N1

      Call Form_a_matrix(x(iarray_in1),1,n1,1,n2,x(iarray_in2),1,
     &                   m1,1,m2,Ndim)

#ifdef _DEBUG_LVL0
      If (me .eq. 0) then
      m=0
      Do i=1,Ndim
         m1=i+(i-1)*ndim
         if (Dabs(x(m1+iarray_in2-1)) .GT. 0.0) m=m+1
      Enddo 
      Write(6,*) "RHS,LHS", rhs, lhs
      If (rhs .eq. 1) then 
      Write(6,"(a)") "Returning right eigenvectors"
      Call output(X(iarray_in2),1,m,1,m,Ndim,Ndim,1)
      Write(6,"(a)") "Returning eigenvalues"
      Call output(X(iarray_in1),1,m,1,m,Ndim,Ndim,1)
      Write(6,*) "End printing from the SIP"
      Elseif (lhs.eq.1) then
      Write(6,"(a)") "Returning left eigenvectors"
      Call output(X(iarray_in2),1,m,1,m,Ndim,Ndim,1)
      Write(6,"(a)") "Returning eigenvalues"
      Call output(X(iarray_in1),1,m,1,m,Ndim,Ndim,1)
      Write(6,*) "End printing from the SIP"
      Endif
      Endif
#endif

      return
      end
C
      Subroutine Form_A_Matrix(AM,i1,i2,j1,j2,AP,a1,a2,b1,b2,Ndim)

      implicit none 
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'aces3_vars.h'

      Integer i,j,k,l,i1,i2,j1,j2,a1,a2,b1,b2,Ndim
      Integer k1,k2,m1,m2,n,Ierror,Neg

      DOuble Precision Thres,Eig_I
      Double precision AM(i1:i2,j1:j2),AP(a1:a2,b1:b2)
      Double precision Z(Ndim,Ndim),U(Ndim,Ndim)
      Double Precision EigR(Ndim),EigI(Ndim),Work(8*Ndim)
      Double Precision EvecR(Ndim,Ndim),EvecL(Ndim,Ndim)

      Parameter (Thres = 1.0D-09)
C N is only use for printing (and to check wheather the actual 
C dimension is one (the lead dimension is always what it was 
C defined to be in the SIAL code. 

      N = 0
      Do i=1,Ndim
         if (Dabs(ap(i,i)) .GT. 0.0) N=N+1
      Enddo 
#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a,i4)")"@-Entry: Form_A_Matrix,Ndim: ",Ndim
      Write(6,"(a)") "The AM=(A-B) small matrix"
      call output(AM,1,N,1,N,Ndim,Ndim,1)
      Write(6,"(a)") "The AP=(A+B) small matrix"
      call output(AP,1,N,1,N,Ndim,Ndim,1)
      EndiF
#endif
      If (N .Eq. 1) Then
          AP(1,1) = 1.0D0
          Return
      ENdif 

      Call Dgeev("V","V",Ndim,AM,Ndim,EigR,EigI,EvecL,Ndim,
     &            EvecR,Ndim,Work,4*Ndim,Ierror)
      
      Neg = 0
      Do I = 1, Ndim
         If (EIgR(I) .LT. 0.0D0) Neg = Neg + 1
      Enddo
     
      If (Neg .Gt. 0) Then
         If (me .eq. 0) Then
         Write(6,"(a,a)") "@-form_rpa_small_matrix: The RPA (A-B):", 
     &                    " matrix is not positive definite."
         Endif 
         Call abort_job()
      Endif 
      
#ifdef _DEBUG_LVL0
      If (me .eq.0) then
          Write(6,"(a)") "Real and imaginary eigenvalues of AM"
          Write(6,"(6(1x,F12.6))") (EigR(i),i=1,N)
          Write(6,"(6(1x,F12.6))") (EigI(i),i=1,N)  
          Write(6,"(a)") "The left eigenvectors"
          call output(EvecL,1,N,1,N,Ndim,Ndim,1)
          Write(6,"(a)") "The right eigenvectors"
          call output(EvecR,1,N,1,N,Ndim,Ndim,1)
       Endif 
#endif
C Since Ndim is at most 50, lets work with the entire matrix. If
C we want to work with the actuall matrix then SIAL code need to
C pass in the iteration number (which is not easy with the number
C of arguments limited to 2). 

C Form (A-B)^{1/2}

      Call Zero (Z,Ndim*Ndim)
      Do I=1,Ndim
         If (EigR(I) .LT. Thres) Then
             Z(i,i) = 0.0D0
         Else 
             Z(i,i) = Dsqrt(EigR(I))
         Endif 
      Enddo 

      Call Zero(AM,Ndim*Ndim)
      Call Zero(U,Ndim*Ndim)
      Call Dgemm("N","N",Ndim,Ndim,Ndim,1.0D0,EvecL,Ndim,Z,
     &              Ndim,0.0D0,U,Ndim)
      Call Dgemm("N","T",Ndim,Ndim,Ndim,1.0D0,U,Ndim,EvecR,
     &              Ndim,0.0D0,AM,Ndim)

#ifdef _DEBUG_LVL0
      if (me .eq. 0) then
      Write(6,"(a)") "Check AM^(1/2) by AM=AM^(1/2)AM(1/2)matrix "
      Call Dgemm("N","N",Ndim,Ndim,Ndim,1.0D0,AM,Ndim,AM,
     &              Ndim,0.0D0,Z,Ndim)
      Call output(Z,1,N,1,N,Ndim,Ndim,1)
      endif 
#endif 
C If left hand-side save the (A+B). See below.

      If (LHS .EQ. 1) Then
         Call Dcopy(Ndim*Ndim,AP,1,U,1)
      Endif 

C Form AM^{1/2}}(AP^{1/2}){t}AM^{1/2}

      Call Zero(Z,Ndim*Ndim)
      Call Dgemm("N","N",Ndim,Ndim,Ndim,1.0D0,AM,Ndim,AP,
     &            Ndim,0.0D0,Z,Ndim)
      Call Zero(AP,Ndim*Ndim)
      Call Dgemm("N","N",Ndim,Ndim,Ndim,1.0D0,Z,Ndim,AM,
     &            Ndim,0.0D0,AP,Ndim)

#ifdef _DEBUG_LVL0
      If (me .eq. 0) then
      Write(6,"(a)") "AM^(1/2)AP^(1/2) matrix "
      Call output(AP,1,N,1,N,Ndim,Ndim,1)
      Endif
#endif
      Call Dgeev("V","V",Ndim,AP,Ndim,EigR,EigI,EvecL,Ndim,
     &            EvecR,Ndim,Work,4*Ndim,Ierror)

#ifdef _DEBUG_LVL0
        If (me .eq. 0) then
        Write(6,"(a,a)") "Real and imaginary eigenvalues of",
     &                   "AM^{1/2}}(AP^{1/2}){t}AM^{1/2}"
        Write(6,"(6(1x,F12.6))") (EigR(i),i=1,N)
        Write(6,"(6(1x,F12.6))") (EigI(i),i=1,N)
        Endif 
#endif

C Form right eigenvector: R=(A-B)^(-1/2)(X+Y). Multiply this 
C by (A-B)1/2 to obtain X+Y

      Call Zero(Z,Ndim*Ndim)
      Call Dgemm("N","N",Ndim,Ndim,Ndim,1.0D0,AM,Ndim,EvecR,
     &            Ndim,0.0D0,Z,Ndim)

C Form the left eiegnvector.L=(A+B)R, The matrix U contains the
C (A+B), see above.

      If (LHS .EQ. 1) Then

         Call Zero(EvecL,Ndim*Ndim)
         Call Dgemm("N","N",Ndim,Ndim,Ndim,1.0D0,U,Ndim,Z,
     &            Ndim,0.0D0,EvecL,Ndim)
         Do I = 1, Ndim
            Eig_I = 1.0D0/Dsqrt(Dabs(EigR(I)))
            Call Dscal(Ndim,Eig_I,EvecL(1,I),1)
         Enddo
      Endif 

C Return the eigenvalues in the first argument and eigenvectors in the second.

      Call Zero(AM,Ndim*Ndim)
      Do i = 1, Ndim
         AM(i,i) = Eigr(i)
      Enddo

      If (RHS .EQ. 1 .AND. LHS .EQ. 0) Then
         Call Dcopy(Ndim*Ndim, Z, 1, Ap,1)
      Elseif (RHS .EQ. 0 .AND. LHS .EQ. 1) Then
         Call Dcopy(Ndim*Ndim, EvecL, 1, Ap,1)
      Endif

      Return
      End

