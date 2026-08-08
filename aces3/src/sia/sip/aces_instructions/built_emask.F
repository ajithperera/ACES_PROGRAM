
      subroutine built_emask(
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
      integer array_in,index, nindex_in,nindex_out, ierr
      integer block, blkndx, seg
      integer find_current_block
      integer array_in_type, array,array_type
      integer*8 indblk1, indblk2, get_block_index
      integer stack,ind_in(mx_array_index),ind_out(mx_array_index)
      integer nsend,sind,nindex 

      integer comm

      integer fop1(mx_array_index), fop2(mx_array_index)
      integer sop1(mx_array_index), sop2(mx_array_index)
      integer sindex(mx_array_index), findex(mx_array_index)
      integer type
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
C The incomming array and the scalar is not used. They are only to
C make the sial compiler happy (i.e. to avoid having to call SIPs
C without arguments)

      array_in = op(c_result_array)
      array_in_type   = array_table(c_array_type, array_in)
C
      if (array_in_type .ne. static_array) then
         Write(6,"(a)") "@-BULID_EMASK: Array must be static."
         call abort_job()
      endif

      nindex_in  = array_table(c_nindex, array_in)

      do i = 1, nindex_in
         ind_in(i) = array_table(c_index_array1+i-1,array_in)
      enddo

       n1 = index_table(c_index_size, ind_in(1))  ! pick up length of index
       n2 = index_table(c_index_size, ind_in(2))  ! pick up length of index

c   Look up each array's address.

      iarray_in  = get_index_from_base(address_table(array_in), x, 2)

      array       = op(c_op1_array)
      array_type  = array_table(c_array_type, array)
      nindex      = array_table(c_nindex, array)

      if (array_type .ne. scalar_value) then
         Write(6,"(a,a)") "@-BUILD_EMASK Second argument must be",
     +                  " an scalar."
         call abort_job()
      endif

      sind =  array_table(c_scalar_index, array)
      type = INT (scalar_table(sind))

      If (Type .EQ. -1) Then
         Call Setup_emask_1(x(iarray_in),1,n1,1,n2,Type)
         scalar_table(sind) = Dble(Type)
      Elseif (Type .GT. 0) then
         Call Setup_emask_2(x(iarray_in),1,n1,1,n2,Type)
      Endif 

      return
      end
C
      Subroutine Setup_emask_1(Data,i1,i2,j1,j2,Ndim)

      implicit none 
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'ex_mask_vars.h'
      include 'proto_events.h'
      include 'epsilon.h'
      include 'mpif.h'

      Integer i,i1,i2,j1,j2,Ndim
      Integer Ipick,Nocc,Iend,Ispin,Ncalc_a,Ncalc_b
      Integer Ncalc
      Double Precision Thres,Emin_prev,Emin_curr

      Double Precision Data(i1:i2,j1:j2)
     
      Nocc = Max(Nalpha_occupied,Nbeta_occupied)
      If (Nocc .Gt. Max_nocc) Then
         Write(6,"(a,a,i4,a)") " The current max no. of occupied",
     +           "orbitals", Nocc, " is greater than the maximum"
         Write(6,"(a,i4,a)") "allowed", Max_nocc, "occupied functions"
         call abort_job()
      Endif

#ifdef _DEBUG_LVL0
      If (me .eq.0) Then
      Write(6,"(a,2I4)") " Variables, Nbasis, Nalpha_occupied: ",
     +                    Nbasis,Nalpha_Occupied
      Endif
#endif 
C Determine the total number of core sites in the molecule. The number
C of core sites are determined by identifying how many orbital 
C energies are below 4 eV. Start from the lowest. 

      Thres     = 2.0D0/27.2113957D0
      Do Ispin = 1, 1+Iuhf
         If (Ispin .EQ. 1) Ncalc_a = 0
         If (Ispin .EQ. 2) Ncalc_b = 0
         Emin_prev = -1.0D08
 100     Emin_curr =  1.0D02
         If (Ispin .Eq. 1) Iend = Nalpha_occupied 
         If (Ispin .Eq. 2) Iend = Nbeta_occupied
         Do I = 1, Iend 
            If (Ispin .Eq. 1) Then
               If (Epsilon(I) .Gt. Emin_prev .And. Epsilon(I) 
     +         .lt. Emin_curr) Emin_curr = Epsilon(I)
            Else 
               If (Epsilonb(I) .Gt. Emin_prev .And. Epsilonb(I) 
     +        .lt. Emin_curr) Emin_curr = Epsilonb(I)
            Endif 
         Enddo

         If (Emin_curr .Lt. -4.0D0) Then
            If (Ispin .EQ. 1) Ncalc_a = Ncalc_a + 1
            If (Ispin .EQ. 2) Ncalc_b = Ncalc_b + 1
            Emin_prev = Emin_curr + Thres
         Endif
         If (Emin_curr .Lt. -4.0D0) Go to 100
      Enddo 

      If (Iuhf .EQ. 0) Ncalc_b = Ncalc_a

#ifdef _DEBUG_LVL0
      If (me .eq. 0) Then
      Write(6,"(a,2I3)") " The number of alpha and Beta core-sites: ",
     +                    Ncalc_a, Ncalc_b
      Endif
#endif 
      
      If (Iuhf .EQ. 0 .AND. (Ncalc_a .NE. Ncalc_b)) Then
         If (Me .EQ. 0) Then
            Write(6,*) 
            Write(6,"(3a)") " Unequal alpha and beta core sites.", 
     +                      " Maximum is used as the number of",
     +                      " core sites."
         Endif 
      Endif 

      If (Iuhf .EQ. 0) Then
         Nsites = Ncalc_a
      Else
         Nsites = Max(Ncalc_a,Ncalc_b)
      Endif 

      If (Me .Eq. 0) Then
      Write(6,"(a,I3)")' Number of different core-sites: ', Nsites
      Write(6,"(a,a)") " Number of core-sites are determined by",
     +                 " counting the number of consecutive"
      Write(6,"(a,a)") " eigenvalues that are within 2 eV. Note",
     +                 " that counting start from the lowest"
      Write(6,"(a)")   " eigenvalue."
      Write(6,*)
      Endif

C Return the number of sites to the sial code 

      Ndim = Nsites 

      Return
      End 

      Subroutine Setup_emask_2(Data,i1,i2,j1,j2,Entry)


      implicit none
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'ex_mask_vars.h'
      include 'proto_events.h'
      include 'epsilon.h'
      include 'mpif.h'

      Integer i,i1,i2,j1,j2,Entry 
      Integer Ipick,Nocc,Iend,Ispin,Ncalc_a,Ncalc_b
      Integer Ncalc
      Double Precision Thres,Emin_prev,Emin_curr

      Double Precision Data(i1:i2,j1:j2)

      Nocc = Max(Nalpha_occupied,Nbeta_occupied)
      If (Nocc .Gt. Max_nocc) Then
         Write(6,"(a,a,i4,a)") " The current max no. of occupied",
     +           "orbitals", Nocc, " is greater than the maximum"
         Write(6,"(a,i4,a)") "allowed", Max_nocc, "occupied functions"
         call abort_job()
      Endif

#ifdef _DEBUG_LVL0
      If (me .eq.0) Then
      Write(6,"(a,2I4)") " Variables, Nbasis, Nalpha_occupied: ",
     +                    Nbasis,Nalpha_Occupied
      Endif
#endif
      Do I = 1, Max_nocc
         Ibelong_a(I) = 0
         Ibelong_b(I) = 0
      Enddo

C Pick the core-site that the user wants to be processed. This must
C eventually need to be read from the input file. There are several
C ways one can think of building the excitation mask. At this point
C this simply use the principal hole to particle excitation.
    
      Ipick  = Entry
     
C Now determine the the core orbitals belongs to this site. The orbitals
C belong to the same site if they are withing 2eV.

      Thres     = 2.0D0/27.2113957D0
      Do Ispin = 1, 1+Iuhf
         Ncalc = 0
         Emin_prev = -1.0D08
 200     Emin_curr =  1.0D02
         If (Ispin .Eq. 1) Iend = Nalpha_occupied 
         If (Ispin .Eq. 2) Iend = Nbeta_occupied 
         Do I = 1, Iend 
            If (Ispin .Eq. 1) Then
               If (Epsilon(I) .Gt. Emin_prev .And. Epsilon(I)
     +         .lt. Emin_curr) Emin_curr = Epsilon(I)
            Else
               If (Epsilonb(I) .Gt. Emin_prev .And. Epsilonb(I)
     +        .lt. Emin_curr) Emin_curr = Epsilonb(I)
            Endif
         Enddo
         If (Emin_curr .Lt. -4.0D0) Then
             Ncalc = Ncalc + 1
             If (Ncalc .Lt. Ipick) Then
                Emin_prev = Emin_curr + Thres
             Endif 
          Endif
         If (me.eq.0) write(6,"(I2,F15.10)") Ncalc, Emin_curr
         If (Ncalc .Lt. Ipick) Go to 200

         Do I = Ipick, Iend 
            If (Ispin .Eq. 1) Then
               If ((Epsilon(I) - Emin_curr) .Lt. Thres) Then
                  Ibelong_a(I) = 1
                  If (Iuhf .Eq. 0) Ibelong_b(I) = 1
               Endif 
            Elseif (Ispin .Eq. 2) Then
               If ((Epsilonb(I) - Emin_curr) .Lt. Thres) Then
                  Ibelong_b(I) = 1
               Endif 
            Endif 
         Enddo
      Enddo 

C#ifdef _DEBUG_LVL0
      If (me .eq.0) then
      Write(6,*) 
      Write(6,"(a)") " The markers for the excitation mask"
      Write(6,"(10I4)") (Ibelong_a(i),i=1, Nalpha_occupied)
      Write(6,"(10I4)") (Ibelong_b(i),i=1, Nbeta_occupied)
      write(6,*)
      endif 
C#endif 

      Return
      End

