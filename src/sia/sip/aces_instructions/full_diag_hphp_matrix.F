
      subroutine Full_diag_HPHP_matrix(
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
      include 'parallel_info.h'
      include 'int_gen_parms.h'
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

      integer i, j, k,n1,n2,n3,n4,m1,m2,m3
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
      integer Idim,Jdim,Kdim,Ldim,Mdim,Ndim,Nroots 
      Integer nocc_nvrt_a,nocc_nvrt_b
      Integer nocc_alpha,nvrt_alpha,nocc_beta,nvrt_beta
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
         Write(6,"(a,a)") "@-FULL_DIAG_HPHP_MATRIX: Both arrays",
     &                  " must be static."
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

       n1 = index_table(c_index_size, ind_in(1))  
       n2 = index_table(c_index_size, ind_in(2)) 

       m1 = index_table(c_index_size, ind_out(1))  
       m2 = index_table(c_index_size, ind_out(2)) 
       m3 = index_table(c_index_size, ind_out(3)) 

c---------------------------------------------------------------------------
c   Look up each array's address.
c---------------------------------------------------------------------------

      iarray_in  = get_index_from_base(address_table(array_in ), x, 2)
      iarray_out = get_index_from_base(address_table(array_out), x, 2)

      Idim   = n1
      Jdim   = n2

      Kdim   = m1
      Ldim   = m2
      Mdim   = m3
      Nroots = eom_roots 

      nocc_alpha = nalpha_occupied
      nvrt_alpha = nalpha_virtual
      nocc_beta  = nbeta_occupied
      nvrt_beta  = nbeta_virtual

      nocc_nvrt_a = nocc_alpha * nvrt_alpha
      nocc_nvrt_b = nocc_beta  * nvrt_beta
      Ndim        = nocc_nvrt_a + nocc_nvrt_b
     
#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,*)
      Write(6,"(a,5(1x,i4))") "@-Entry to full_diagonalize variables:",
     +                   Idim,jdim,KDim,Ldim,Mdim
      Endif 
#endif 
      If (nocc_nvrt_a .NE. nocc_nvrt_b) Then
         Write(6,"(a,a)") "Open shell CIS method is not implemented:", 
     +                    "code modifications are required to process"
         Write(6,"(a)")   " open shell CIS calculations."
         Call abort_job()
      Endif

      If (Idim .NE. Jdim) Then
         Write(6,"(a,a)") "Array dimension mismatch: The source", 
     +                    " matrix must be square." 
         Call abort_job()
      Endif 

      If (Mdim .Lt. Nroots) Then 
         Write(6,"(2a)") "Column dimension of the target array must", 
     +                  " be grater or equal to the number of roots."
         Call abort_job()
      Endif 

      Call Full_Diag_phph(x(iarray_in),Idim,Jdim,x(iarray_out),
     +                    Kdim,Ldim,Mdim,Ndim,Nroots)

      Return
      End

      Subroutine Full_diag_phph(AMAT,IJdim,KLdim,Evecs,Kdim,Ldim,
     +                          Mdim,Ndim,Nroots)

      implicit none 
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      Integer Nroots,IJdim,KLdim,Kdim,Ldim,Ndim,Mdim,MN
      Integer I,J,K,M,N,Idim,Jdim,Ierror 
      Integer nocc_alpha,nvrt_alpha,nocc_beta,nvrt_beta
      Integer T_count,S_count 

      Double precision AMAT(IJdim,KLdim)
      Double Precision EVECS(Kdim,Ldim,MDim)
      Double precision EvecR(Ndim,Ndim),EvecL(Ndim,Ndim)
      Double Precision EigR(Ndim),EigI(Ndim)
      Double Precision Work(4*Ndim)
      Double Precision Esave,Csave,Dtmpr,Dtmpl,dtmp,sum
      Double Precision Threshold

C Lets use generalized eigenvector solver. 

      nocc_alpha = nalpha_occupied
      nvrt_alpha = nalpha_virtual
      nocc_beta  = nbeta_occupied
      nvrt_beta  = nbeta_virtual
      Threshold = eom_tol 

      If (Ndim .GT. Ijdim) Then
         Write(6, "(2a,I8,2a,I8)") " The actual dimension of the CIS",
     +                    " matrix",Ndim," is greater than the lead",
     +                    " dimension",Ijdim
         Call abort_job()
      Endif 

C#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a,a,1x,i5,1x,i5)") "@-Entry to Full_diag_phph",
     +                       " order of the matrix", IJDim, Ndim
C      call output(Amat,1,Ndim,1,Ndim,IJdim,IJdim,1)
      Endif 
C#endif 

      Call Dgeev("V","V",Ndim, Amat, IJdim, EigR, EigI,
     +            EvecL,Ndim,EvecR,Ndim,Work,4*Ndim,Ierror)

C Copy very small eigenvalues to the end. Do the same thing for
c eigenvectors. This is really to handle nearly zero eigenvalues.

#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,*) 
      Write(6,"(a)") "AMAT(IJ,KL) eigenvalues/vectors"
      Write(6,"(a)") "The right-hand side"
      Write(6,"(6(1x,F10.5))") (EigR(i),i=1,Ndim)
C      Write(6,"(6(1x,F10.5))") (EigI(i),i=1,Ndim)
C      call output(EvecR,1,Ndim,1,Ndim,Ndim,Ndim,1)
C      Write(6,"(a)") "The left-hand side"
C      call output(EvecL,1,Ndim,1,Ndim,Ndim,Ndim,1)
      endif 
#endif 

      If (Ndim .Ne. 1) Then

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
    
      Endif 

#ifdef _DEBUG_LVL0
      if (me.eq.0) then
      Write(6,*)
      Write(6,"(a,a)") "@-Diagonalize: The reordered eigenvalues",
     +               " and vectors"
      Write(6,"(6(1x,F10.5))") (EigR(i),i=1,Ndim)
      Write(6,"(a)") "The right-hand side"
      call output(EvecR,1,Ndim,1,Ndim,IJdim,Kldim,1)
      Endif
#endif
     

C Locate the nroot number of singlet and triplet lowest energy 
C root. Store the singlets in the first Nroot diagonals of the
c source matrix  and copy the corresponding alpha eigenvectors 
C to the first nroot vectors of target matrix. The triplets
C would be stored in next nroot elements. 

      S_count = 0
      T_count = 0
      Do I =1, Ndim 
          
         Sum = 0.D0
         Do J =1, Ndim 
            Sum = Sum +  EvecR(J,I)
         Enddo 

         If (Dabs(Sum) .Lt. Threshold) Then

C Triplet states 

            T_count = T_count + 1
            If (T_count .LE. Nroots) Then
               Amat(T_count+Nroots, T_count+Nroots) = EigR(I)
               MN = 0
               DO M=1, nocc_alpha
               Do N=1, nvrt_alpha 
                  MN = MN + 1
#ifdef _DEBUG_LVL0
         If (me.eq.0) write(6,"(3i3,F15.10)")m,n,mn,EvecR(MN,I)
#endif 
                  Evecs(N,M,T_count+Nroots) = EvecR(MN,I)
               Enddo
               Enddo
            Endif 

         Else

C Singlet states

            S_count = S_count + 1
            If (S_count .LE. Nroots) Then
               Amat(S_count,S_count) = EigR(I)
               MN = 0
               DO M=1, nocc_alpha
               Do N=1, nvrt_alpha 
                  MN = MN + 1
                  Evecs(N,M,S_count) = EvecR(MN,I)
#ifdef _DEBUG_LVL0
         If (me.eq.0) write(6,"(3i3,F15.10)")m,n,mn,EvecR(MN,I)
#endif 
               Enddo
               Enddo
            Endif

         Endif 
         If (S_count .EQ. Nroots .and. T_count .EQ. Nroots) 
     +       Go to 10
      Enddo

  10  Continue 

#ifdef _DEBUG_LVL0
      if (me.eq.0) then
          Write(6,"(a)") "The selected CIS eigenvalues"
          Write(6,"(6(1x,F15.6))") (Amat(i,i),i=1,2*Nroots)
CSSS          Write(6,"(a)") "The selected CIS eigenvectors"
CSSS          Do i =1, 2*nroots 
CSSS          Call output(Evecs(1,1,i),1,Nvrt_alpha,1,Nocc_alpha,Kdim,
CSSS     +                Ldim,1)
CSSS          Enddo
      endif 
#endif 

      Return
      End

