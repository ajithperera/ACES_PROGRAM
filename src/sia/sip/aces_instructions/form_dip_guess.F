
      subroutine Form_DIP_guess(
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
      integer Idim,Jdim,Kdim,Ldim,IJdim,KLdim
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
         Write(6,"(a)") "@-FORM_DIP_GUESS: Both arrays must be static."
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
       n3 = index_table(c_index_size, ind_in(3)) 
       n4 = index_table(c_index_size, ind_in(4)) 

       m1 = index_table(c_index_size, ind_out(1))  
       m2 = index_table(c_index_size, ind_out(2)) 
       m3 = index_table(c_index_size, ind_out(3)) 

c---------------------------------------------------------------------------
c   Look up each array's address.
c---------------------------------------------------------------------------

      iarray_in  = get_index_from_base(address_table(array_in), x, 2)
      iarray_out = get_index_from_base(address_table(array_out), x, 2)

      Idim = n1
      Jdim = n2
      Kdim = n3
      Ldim = n4
     
      IJdim = n1*n2
      KLdim = n3*n4
#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a,9i4)") "@-Entry to pack4d_to_2D variables",
     +                IJDim,KLDim,Idim,jdim,kdim,ldim,m1,m2,m3
      Endif 
#endif 

      Call pack4d_to_2D(x(iarray_in),Idim,Jdim,Kdim,Ldim,IJdim,KLdim,
     &                  x(iarray_out),m1,m2,m3)

      Return
      End

      Subroutine pack4d_to_2D(Hbar,Idim,Jdim,Kdim,Ldim,IJdim,KLdim,
     &                        Guess_vecs,m1,m2,m3)

      implicit none 
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      Integer I,J,K,L,P,Idim,Jdim,Kdim,Ldim,IJdim,KLdim
      Integer IJ,KL,m1,m2,m3,Nroots 

      Double Precision Hbar(Idim,Jdim,kdim,Ldim) 
      Double Precision Hbar_packed(IJdim,KLdim)
      Double Precision Guess_Vecs(m1,m2,m3)

      If (IJdim .NE. KLdim) Then
         If (me. eq.0) then
         Write(6,"(a,a)") "@-FORM_DIP_GUESS: The matrix is ",
     +                    " not a square matrix"
         Endif 
         Call abort_job()
      Endif 
#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a,9i4)") "@-pack4d_to_2D variables;before packing",
     +                IJDim,KLDim,Idim,jdim,kdim,ldim,m1,m2,m3
      Endif 
#endif 
      
      IJ = 0
      KL = 0
      DO J =1, Jdim
      DO I =1, Idim
         IJ = IJ + 1
         KL = 0
         DO L =1, Ldim
         DO K =1, Kdim
            KL = KL + 1
            Hbar_packed(KL,IJ) = Hbar(K,L,I,J)
         ENDDO
         ENDDO
      ENDDO
      ENDDO

      Nroots = eom_roots

#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a,7i4)") "@-Entry to Get_dip_guess order of the matrix",
     +                IJDim,KLDim,Nroots,Idim,jdim,kdim,ldim
      Endif
#endif 

      Call Get_dip_guess(Hbar_packed,IJdim,KLdim,Nroots,Hbar,
     &                   Guess_Vecs,m1,m2,m3,Idim,Jdim)
#ifdef _DEBUG_LVL0
      if (me.eq.0) then
      do p=1,Nroots
      write(6,*)
      Write(6,"(a)") "Initial Right guess vectors for EOM-DIP/EA"
      call output(Guess_Vecs(1,1,p),1,Idim,1,Jdim,Idim,Jdim,1)
      enddo
      endif
#endif
          
      Return
      End
C
      Subroutine Get_dip_guess(Hbar_packed,IJdim,KLdim,Nroots,
     &                         Guess_intmd,Guess_vecs,m1,m2,m3,
     &                         Idim,Jdim)

      implicit none 
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      Integer Nroots,IJdim,KLdim,Ndim,Nill,Nbgn,Nend
      Integer I,J,K,P,IJ,m1,m2,m3,Idim,Jdim,Ierror 

      Double precision Hbar_packed(IJdim,KLdim)
      Double Precision Guess_vecs(m1,m2,m3)
      Double Precision Guess_intmd(IJdim,KLdim)
      Double precision EvecR(IJdim,KLdim),EvecL(IJdim,KLdim)
      Double Precision EigR(IJdim),EigI(IJdim)
      Double Precision Work(4*IJdim)
      Double Precision Thres,Esave,Csave,Dtmpr,Dtmpl,dtmp

C Hbar is non-symmetric matrix, and we obtain both left/right eigenvectors.
C It is square (always). 

      Ndim=IJdim
      Thres=1.0D-06

#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a,1x,i4)") "@-Entry to Get_guess order of the matrix", 
     +                IJDim
      Endif 
#endif 

      Call Dgeev("V","V",Ndim, Hbar_packed, Ndim, EigR, EigI,
     +            EvecL,Ndim,EvecR,Ndim,Work,4*Ndim,Ierror)

C Copy very small eigenvalues to the end. Do the same thing for
c eigenvectors. This is really to handle nearly zero eigenvalues.

#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a)") "H(IJ,KL) eigenvalues/vectors"
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
CSSS     
CSSS      Call Process_vecs(Eigr,EvecR,EvecL,Guess_intmd,Ndim,
CSSS     +                  1,Ndim,1,Ndim)

#ifdef _DEBUG_LVL0
      if (me.eq.0) then
      Write(6,*)
      Write(6,"(a,a)") "@-FORM_DIP_GUESS: The reordered eigenvalues",
     +               " and vectors"
      Write(6,"(6(1x,F10.5))") (EigR(i),i=1,Ndim)
      Write(6,"(a)") "The right-hand side"
      call output(EvecR,1,Ndim,1,Ndim,Ndim,Ndim,1)
      Endif
#endif
C
      Nill = 0
      Do i = 1, Ndim
         If (Dabs(EigR(i)) .Lt. Thres) Then
            Nill  = Nill + 1
         Endif
      Enddo

      Nbgn = Nill+1 
      Nend = Nbgn+Nroots 

CSSS      if (me.eq.0) write(6,*) Nbgn,Nend,Nill,Ndim-Nill,nroots 

      If (Nend .Gt. Nroots) Then

      k=0
      Do P=Nbgn, Nend 

#ifdef _DEBUG_LVL0
         if (me.eq.0) write(6,*)
#endif 
         IJ = 0
         K  = K + 1
         DO I=1, Jdim 
         DO J=1, Idim 
            IJ = IJ + 1
            Guess_vecs(I,J,k) = EvecR(IJ,P)

#ifdef _DEBUG_LVL0
       if (me.eq.0) Write(6,"(1x,2F10.5,3i4)") Guess_intmd(IJ,P),
     +              Guess_vecs(I,J,P),ij,i,j
#endif 
         ENDDO
         ENDDO
      ENDDO 

      Else
         Write(6,"(a)") "No non-zero guess eigenvalues" 
      Endif 

#ifdef _DEBUG_LVL0
      if (me.eq.0) then
      do p=1,Nroots
      Write(6,"(a)") "Initial Right guess vectors for EOM-DIP/EA"
         DO I=1, Jdim 
         DO J=1, Idim 
            if (me.eq.0) Write(6,"(1x,F10.5)") Guess_vecs(I,J,p)
         ENDDO
         ENDDO
      enddo 
      endif 
#endif 

      Return
      End

