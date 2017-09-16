
      subroutine Genrl_eig_solver(
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
      integer*8 iarray_in1,iarray_in2,ievec,daddr,iarray_ca

      integer i, j, k,n1,n2,m1,m2,l1,l2
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
         Write(6,"(a)") "GNRL_EIG_SOLVER: Both arrays must be static."
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
c----------------------------------------------------------------------------
c   Get the Ca coeffcients. 
c----------------------------------------------------------------------------

      array_ca = 1
      array_ca_type = array_table(c_array_type, array_ca)

      if (array_ca_type .ne. static_array) then
         Write(6,"(a)") "GNRL_EIG_SOLVER: CA array must be static."
         call abort_job()
      endif

      nindex_ca = array_table(c_nindex, array_ca)

      if (nindex_ca .ne. 2) then
         Write(6,*) "GNRL_EIG_SOLVER: Ca array must have two indices"
         call abort_job()
      endif

      do i = 1, nindex_ca
         ind_ca(i) = array_table(c_index_array1+i-1,array_ca)
      enddo
      l1 = index_table(c_index_size, ind_ca(1)) 
      l2 = index_table(c_index_size, ind_ca(2))  

      iarray_ca = get_index_from_base(address_table(array_ca), x, 2)

      Call Solve_gnrl_eigvalue(x(iarray_in1),1,n1,1,n2,x(iarray_in2),1,
     &                         m1,1,m2,x(iarray_ca),1,l1,1,l2,Ndim)

      return
      end
C
      Subroutine Solve_gnrl_eigvalue(H,i1,i2,j1,j2,D,a1,a2,b1,b2,
     &                               Ca,k1,k2,m1,m2,Ndim)
      implicit none 
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      Integer i,j,k,l,i1,i2,j1,j2,a1,a2,b1,b2,Ndim,Nroots 
      Integer k1,k2,m1,m2,izero,Nfin
      Integer Istart, Ndim_modf, Ierror,Ilo,Ilh,Neg,Nill
      Integer Vec_max(Ndim),IfailR(Ndim),failL(Ndim)

      Logical Select(Ndim)

      DOuble Precision Thres,Esave,Csave,Trace,MaxR,Conv
      Double precision H(i1:i2,j1:j2),D(a1:a2,b1:b2),Ca(k1:k2,m1:m2)
      Double precision Z(Ndim,Ndim),U(Ndim,Ndim)
      Double precision Q(Ndim,Ndim)
      Double Precision Tau(Ndim)
      Double Precision Lscale(Ndim),Rscale(Ndim),Work(8*Ndim)
      Double Precision Work1((Ndim+2)*Ndim),EigR(Ndim),EigI(Ndim)
      Double Precision EvecR(Ndim,Ndim),EvecL(Ndim,Ndim)
      Double Precision Beta(Ndim),Exkt_ip(Ndim)
      Logical Useqz
 
      Parameter (Thres = 1.0D-07,Conv=27.21138505)

C Hbar is non-symmetric matrix, and we obtain both left/right eigenvectors.

#ifdef _DEBUG_LVL0
      If(me.eq.0) Then
      Write(6,"(a,i4)")"@-Entry to Gnrl_eig_solver order of the matrix", 
     +                  Ndim 
      Write(6,"(a)") " The F_ext(p,q)"
      call output(H,1,Ndim,1,Ndim,Ndim,Ndim,1)
      Write(6,*) 
      Write(6,"(a)") " The Density matrix D(p,q)"
      call output(D,1,Ndim,1,Ndim,Ndim,Ndim,1)
      Endif 
#endif
C
C Symmetrize the EXKT matrix (in place)
      
      Do I = 2, Ndim
         Do J = 1, (I-1)
            H(I,J) = 0.50D0*(H(I,J) + H(J,I))
         Enddo
      Enddo
      Do I = 1, Ndim
         Do J= 1, (I-1)
            H(J,I) = H(I,J)
         Enddo
      Enddo 

#ifdef _DEBUG_LVL0
      Trace = 0.0d0
      do i=1, Ndim
      do j=1,Ndim
         if (i .eq. j) Trace = Trace + H(i,j)
      enddo 
      enddo 
      If (me.eq.0) Write(6,"(a,F15.8)") "The trace of F(p,q)", Trace
#endif
C
C Filter out the Noise in Density and Fock matrix.

      Call Filter_out(D,Ndim*Ndim,1.0D-8)
      Call Filter_out(H,Ndim*Ndim,1.0D-8)

      Call Dcopy(Ndim*Ndim,D,1,Q,1)

      Call Dgeev("V","V",Ndim,D,Ndim,EigR,EigI,EvecL,Ndim,
     &            U,Ndim,Work,4*Ndim,Ierror)

      if (me .eq.0) then
      Write(6,*)
      Write(6,"(a,a,a)")" The eigenvalues of the D matrix. The",
     &               " negative or zero means"
      Write(6,"(a)") " density is singular."
      Write(6,"(6(1x,F12.6))") (EigR(i),i=1,Ndim)
      Write(6,*)
      Endif 

      Useqz = .False.
      Neg = 0
      Do I = 1, NDim
         If (EigR(I) .lt. 0.0D0) Neg = Neg + 1
      Enddo 
      If (Neg .Ne. 0) Then
         Write(6,"(a)") " The Density matrix is not positive definite"
         Useqz = .True.
      Endif 

      Useqz = .False.
C
      if (.NOT. Useqz) Then
         If (Me .eq. 0) Write(6,"(a,a)") "The Density matrix is ",  
     &      "not singular: Follow the traditional path."

C Form D^(-1/2)

         Call Zero (Z,Ndim*Ndim)
         Izero = 0
         Do I=1,Ndim
               If (EigR(I) .LT. Thres) Then
                   Z(i,i) = 0.0D0
                   Izero  = Izero + 1
               Else 
                   Z(i,i) = 1.0D0/Dsqrt(EigR(I))
               Endif 
         Enddo 
         Ndim_modf = Ndim - Izero

         If (me .eq.0) Write(6,"(a,a,i4)") "The number of", 
     &                             " singular eig values: ",
     &                      Izero 
#ifdef _DEBUG_LVL0
         if (me .eq. 0) then
         Write(6,"(a)") "D^(-1/2) matrix "
         Call output(Z,1,Ndim,1,Ndim,Ndim,Ndim,1)
         endif 
#endif 

C Form X=UD^(-1/2)

         Call Dgemm("N","N",Ndim,Ndim,Ndim,1.0D0,U,Ndim,Z,
     &              Ndim,0.0D0,EvecL,Ndim)

#ifdef _DEBUG_LVL0
         if (me .eq. 0) then
         Write(6,"(a)") "X matrix "
         Call output(EvecL,1,Ndim,1,Ndim,Ndim,Ndim,1)
         endif
#endif 

C Form F'=X^tFX

         Call Dgemm("T","N",Ndim,Ndim,Ndim,1.0D0,EvecL,Ndim,H,
     &              Ndim,0.0D0,EvecR,Ndim)
         Call Dgemm("N","N",Ndim,Ndim,Ndim,1.0D0,EvecR,Ndim,EvecL,
     &              Ndim,0.0D0,U,Ndim)
         Call Zero (H,Ndim*Ndim)
         Call Dcopy(Ndim*Ndim,U,1,H,1)

#ifdef _DEBUG_LVL0
         if (me .eq. 0) then
         Write(6,"(a,i4)") "The transformed F-matrix "
         Call output(h,1,Ndim,1,Ndim,Ndim,Ndim,1)
         endif 
#endif
         Call Dgeev("V","V",Ndim,H,Ndim,EigR,EigI,U,Ndim,
     &            EvecR,Ndim,Work,4*Ndim,Ierror)

C Transform the eigenvectors to the original basis.

        Call Dgemm("N","N",Ndim,Ndim,Ndim,1.0D0,EvecL,Ndim,EvecR,Ndim,
     &            0.0D0,U,Ndim)
C Form DC 

        Call Dgemm("N","N",Ndim,Ndim,Ndim,1.0D0,Q,Ndim,U,Ndim,
     &            0.0D0,Z,Ndim)
C Form (DC)^t(DC)

        Call Dgemm("T","N",Ndim,Ndim,Ndim,1.0D0,Z,Ndim,Z,Ndim,
     &            0.0D0,EvecR,Ndim)

        Do i=1, Ndim
           EigI(i) = EvecR(i,i)
        Enddo 

#ifdef _DEBUG_LVL0
         if (me .eq. 0) then
         Write(6,*)
         Write(6,"(a)") " Eigenvalues and vectors before reordering"
         Write(6,"(6(1x,F12.6))") (EigR(i),i=1,Ndim)
C         call output(EvecR,1,Ndim,1,Ndim,Ndim,Ndim,1)
         Write(6,*)
         Write(6,"(a)") " (DC)^tDC: Pole strength"
         Write(6,"(6(1x,F12.6))") (EigI(i),i=1,Ndim)
         endif 
#endif

      Else

         if (me .eq. 0) then
         Write(6, "(a,a)")" The generalized eigenvalue problem",
     &                   " for IP is solved by the QZ algorithm"
         endif 
         Call Dcopy(Ndim*Ndim,Q,1,D,1)
         Call Dggev("V","V",Ndim,H,Ndim,D,Ndim,EigR,EigI,Beta,
     &               EvecL,Ndim,EvecR,Ndim,Work,8*Ndim,Ierror)

         Do i = 1, Ndim
            If (Dabs(Beta(I)) .Gt. Thres) Then
                EigR(i) = EigR(i)/Beta(i)
            Else
                EigR(i) = 0.0D0
            Endif
         Enddo

C Form DC 
        Call Dgemm("N","N",Ndim,Ndim,Ndim,1.0D0,Q,Ndim,U,Ndim,
     &            0.0D0,Z,Ndim)
C Form (DC)^t(DC)

        Call Dgemm("T","N",Ndim,Ndim,Ndim,1.0D0,Z,Ndim,Z,Ndim,
     &            0.0D0,EvecR,Ndim)

        Do i=1, Ndim
           EigI(i) = EvecR(i,i)
        Enddo 

#ifdef _DEBUG_LVL0
         if (me .eq. 0) then
         Write(6,*)
         Write(6,"(a)") " Eigenvalues and vectors before reordering"
         Write(6,"(6(1x,F12.6))") (EigR(i),i=1,Ndim)
         call output(EvecR,1,Ndim,1,Ndim,Ndim,Ndim,1)
         Write(6,*)
         Write(6,"(a)") " (DC)^tDC: Pole strength"
         Write(6,"(6(1x,F12.6))") (EigI(i),i=1,Ndim)
         endif 
#endif

      Endif 

C Reorder the eigenvectors based on the ascending order of (DC)^t(DC) 

      Do i = 1, Ndim
      Do j = i+1, Ndim
         If (EigI(i) .Gt. EigI(j)) Then
            Esave   = EigI(j)
            EigI(j) = EigI(i)
            EigI(i) = Esave
            Esave   = EigR(j)
            EigR(j) = EigR(i)
            EigR(i) = Esave
         Endif
      Enddo
      Enddo

C      If (Useqz) then
C     
C          Do i = 1, Ndim
C          Do j = i+1, Ndim
C             If (EigI(i) .Gt. EigI(j)) Then
C                Esave   = EigR(j)
C                EigR(j) = EigR(i)
C                EigR(i) = Esave
C             Endif
C          Enddo
C          Enddo
C      Endif 

#ifdef _DEBUG_LVL0
         if (me .eq. 0) then
         Write(6,*)
         Write(6,"(a)") " Eigenvalues and PS after reordering"
         Write(6,"(6(1x,F12.6))") (EigR(i),i=1,Ndim)
         Write(6,*)
         Write(6,"(6(1x,F12.6))") (EigI(i),i=1,Ndim)
         endif 
#endif
C Zero out all the orbital energies except for the occupied orbitals 

      Do i =1, Noccupied 
         EigR(i) = EigR(Ndim-i+1)
         EigI(i) = EigI(Ndim-i+1)
      Enddo
C 
      Do i = 1, Noccupied 
      Do j = i+1, Noccupied 
         If (EigR(i) .Gt. EigR(j)) Then
            Esave   = EigR(j)
            EigR(j) = EigR(i)
            EigR(i) = Esave
            Esave   = EigI(j)
            EigI(j) = EigI(i)
            EigI(i) = Esave
         Endif
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
         if (me .eq. 0) then
         Write(6,*)
         Write(6,"(a)") " Eigenvalues and PS after reordering"
         Write(6,"(6(1x,F12.6))") (EigR(i),i=1,Noccupied)
         Write(6,*)
         Write(6,"(6(1x,F12.6))") (EigI(i),i=1,Noccupied)
         endif 
#endif

#ifdef _TRIAL
C Solve the generalized eigenvalue problem

       Call Dgeqrf(Ndim,Ndim,D,Ndim,Tau,Work,Ndim,Ierror)
       Call Dggbal("N",Ndim,H,Ndim,D,Ndim,Ilo,Ilh,Lscale,Rscale,
     &              Work,Ierror) 
       Call Dgghrd("I","I",Ndim,Ilo,Ilh,H,Ndim,D,Ndim,Q,Ndim,
     &              Z, Ndim, Ierror)
       Call Dhseqr("E","N",Ndim,Ilo,Ilh,H,Ndim,EigR,EigI,Z,Ndim,
     &             Work, Ndim,Ierror)
       Call Dhsein("B","Q","N",Select,Ndim,H,EigR,EigI,EvecL,
     &              Ndim,EvecR,Ndim,Ndim,Ndim,Work1,IfailL,IfailR,
     &              Ierror)
#endif 
      If(me.eq.0) Then

         Write(6,*)
         Write(6,"(19x,A)")" The EXKT first ionization Potentials"
         Write(6,"(19x,a)")" ------------------------------------"
         Write(6,5)
         Write(6,10)
         Write(6,5)

         Do i=1, Noccupied
            Write(6,20) -EigR(i),-EigR(i)*Conv,EigI(i)
         Enddo
      Endif 
      Write(6,5)

   5  format(12x,63("-"))
  10  format(12x,'Exkt-IP(a.u) ',7x,'Exkt-IP (eV)',11x,
     +                               "Pole Strength (a.u.)")
  20  format(3x,F20.8,x,F20.8,x,F20.8)
      Write(6,*)

      Return
      End

      SUBROUTINE FILTER_OUT(A,LENGTH,TOL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(LENGTH)
      if (length.lt.1) return
      DO I = 1, LENGTH
         X = ABS(A(I))
         IF (X.LT.TOL) A(I) = 0.0
      END DO
      RETURN
      END
