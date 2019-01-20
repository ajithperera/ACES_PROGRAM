
      subroutine Truncate_eom_subspace(
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
         Write(6,"(a,a)") "@-runcate_eom_subspace: Both arrays must",
     &                    "be static."
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

      Call Find_subspace(x(iarray_in),1,n1,x(iarray_out),1,m1,1,m2)

      return
      end
C
      Subroutine Find_subspace(Eigs,i1,i2,Eigvecs,a1,a2,b1,b2)

      implicit none 
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'aces3_vars.h'
      include 'eom_truncate.h'

      Integer i,j,k,l,i1,i2,j1,j2,a1,a2,b1,b2
      Integer Ndim1,Ndim2,Ndim3,Ndim_test,Ione 
      Integer Istart,Ierror,Iclose,Icount 

      Double Precision Eigs(i1:i2),Eigvecs(a1:a2,b1:b2)
      Double Precision Residual,Xclose 
      Double Precision Scr(Maxdim)

      Data Ione /1/

      Ndim1 = (i2 -i1) + 1
      Ndim2 = (a2 -a1) + 1
      Ndim3 = (b2 -b1) + 1

      If (Ndim1 .Ne. Ndim2 .and. Ndim1 .Ne. Ndim3 .and. Ndim2 
     +    .Ne. Ndim3) Then
      If (Ndim1 .Ne. Ndim_old) 
     +   Write(6,"(a)") " The incoming dimension mismatch"
         Call Abort_job()
      Endif 

      If (Ndim_new .Gt. Maxdim) Then
         Write(6,"(a,I3)") " The truncated subspace dimension must",
     +                      " be less than 20. Currently it is: ",
     +                        Ndim_new
         Call Abort_job()
      Endif 

#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a,I3,1x,I3,1X,F15.7)") "eom_truncate_vars: ", 
     +          ndim_new,Ndim_old,omega 
      Write(6,"(a)") "Eigenvalues/vectors of incoming matrix"
      Write(6,"(6(1x,F10.5))") (Eigs(i),i=1,Ndim1)
      call output(Eigvecs,1,Ndim2,1,Ndim3,Ndim2,Ndim3,1)
      Write(6,"(a,F10.5)") " The current eigenvalue:",omega 
      endif
#endif
c
C Frist find the current eigenvalue and the corresponding 
C vector 

      Icount = 0
      Call Fndclose(Ndim_old,Eigs,Omega,Xclose,Iclose)
      If (Iclose .NE. 0) Then
          Icount = Icount + 1
          Call Dcopy(Ndim_old,Eigvecs(1,Iclose),1,Evec_selct(1,Icount),
     +               1)
          Eval_selct(Icount) = Xclose
          Eigs(Iclose) = 1.0D30
      Else
          If (me .eq. 0) then
          Write(6,"(a,a)") " The current eigenvalue is not in"
     +                     " the list of eigenvalues"
          Endif 
          Call Abort_job()
      Endif 

C Retain the eigenvalues and vectors that are close to omega. Note
C that we only keep the number of vectors of the truncarte subspace
C dimensions,

      Do I = 1, (Ndim_new -1)
         Call Fndclose(Ndim_old,Eigs,Omega,Xclose,Iclose)
         Icount = Icount + 1
         Call Dcopy(Ndim_old,Eigvecs(1,Iclose),1,Evec_selct(1,Icount),
     +              1) 
         Eval_selct(Icount) = Xclose
         Eigs(Iclose) = 1.0D+30
      Enddo 

C Makesure the all the selecteD eigenvalues make sense.

      Ndim_test = Ndim_new
      Do I = Ndim_test, 1, -1
         If (Eval_selct(I) .Gt. 1000.0D0) Then
             Ndim_new = Ndim_new - 1
        Endif
      Enddo 
      If (Ndim_new .Lt. Ione) Then 
        If (Me .eq. 0) Then
        Write(6,"(a)") "The reduced dimension can not be negative.",
     +     "The selected eigenvalues are (all >1000)"
        Write(6,"(6(1X,F15.7))") (Eval_selct(i),i=1,Ndim_old)
        Endif 
        Call Abort_job()
      Endif 

#ifdef _DEBUG_LVL0
      If (Me.eq.0) Then 
      Write(6,"(a)") " The subspace eigenvalues and vectors",
     +               " before orthog."
      Write(6,"(6(1X,F15.7))") (Eval_selct(i),i=1,Ndim_new)
      call output(Evec_selct,1,Ndim_old,1,Ndim_new,Ndim_old,Ndim_new,1)
      endif 
#endif 

C Schmidt Orthogonalize the subspace.

      Do I = 2, Ndim_new
         Call GSchmidt_v2(Evec_selct(1,I),Evec_selct,Ndim_old,
     +                    I-1,Scr,Residual)
      Enddo 

      Call Dcopy(Ndim_new,Eval_selct,1,Eigs,1)
      Do I = 1, Ndim_new
         Call Dcopy(Ndim_old,Evec_selct(1,I),1,Eigvecs(1,I),1)
      Enddo 

#ifdef _DEBUG_LVL0
      if (Me.eq.0) Then 
      Write(6,*) " The subspace wigenvalues and vectors after orthog."
      Write(6,"(6(1X,F15.7))") (Eigs(i),i=1,Ndim_New)
      call output(Eigvecs,1,Ndim_old,1,Ndim_new,Ndim2,Ndim3,1)
      Call Dzero(Scr,Ndim_new*Ndim_new)
      call xgemm("T","N",Ndim_new,Ndim_new,Ndim_old,1.0D0,Evec_selct,
     +            Ndim_old,Evec_selct,Ndim_old,0.0D0,Scr,Ndim_new)
      call output(scr,1,Ndim_new,1,Ndim_new,Ndim_new,Ndim_new,1)
      endif 
#endif `

      Return
      End

      Subroutine Fndclose(Length,V,Test,Vfound,Ifound)

      Implicit None

      Integer Ifound,I,Length
      Double Precision V(Length),Vfound,Test,Xclose
      Double Precision Diff 
      
      Xclose = 1.0D0+30
      Ifound = 0
      Do I = 1, Length
         Diff = Dabs(V(i)-Test)
         If (Diff .LT. Xclose) Then
            Ifound = I
            Vfound = V(I)
            Xclose = Diff 
         Endif
      Enddo
      
      Return
      End

      Subroutine GSchmidt_v2(Vec,Vorth,Nsize,Ndim,Tmp,Residual)

      Implicit None
      Double Precision Vec(Nsize),Vorth(Nsize,Ndim)
      Double Precision Tmp(Nsize)
      Double Precision Residual,One,Fact,Empty,X
      Double Precision Ddot,Dsqrt
      Integer Nsize,Ndim,I

      include 'parallel_info.h'

      Data One,Empty/1.0D0,0.0D0/

      If (Ndim .EQ. 0) Then 
         Residual = Dsqrt(Ddot(Nsize,Vec,1,Vec,1))
         If (Residual .NE. Empty) Then
             X = One/Residual
             Call Dscal(Nsize,X,Vec,1)
          Else
             If (Me.eq.0) Then
             Write(6,"(a)") "Divsion by zero error occurs. The",
     +                      "eigenvectors are not proper"
             Endif 
             call abort_job
          Endif 
          Return
      Endif 

      Call Dcopy(Nsize,Vec,1,Tmp,1)
      Do I = Ndim, 1, -1
         Fact = Ddot(Nsize,Vorth(1,I),1,Vec,1)
         Call Daxpy(Nsize,-Fact,Vorth(1,I),1,Tmp,1)
      Enddo
      Call Dcopy(Nsize,Tmp,1,Vec,1)
      
      Residual = Dsqrt(Ddot(Nsize,Vec,1,Vec,1))
      X = One/Residual
      Call Dscal(Nsize,X,Vec,1)

      Return
      End
