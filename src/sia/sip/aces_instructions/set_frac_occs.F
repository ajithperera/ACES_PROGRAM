
      subroutine Set_frac_occs(
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
      array_in = op(c_result_array)
      array_in_type   = array_table(c_array_type, array_in)
C
      if (array_in_type .ne. static_array) then
         Write(6,"(a)") "@-SET-FRACS_OCCS: Array must be static."
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
         Write(6,"(a)") "@-SET_FRACS_OCCS Must be an scalar."
         call abort_job()
      endif

      sind =  array_table(c_scalar_index, array)
      type = INT (scalar_table(sind))

      Call Load_frac_occs(x(iarray_in),1,n1,1,n2,Type)

      return
      end
C
      Subroutine Load_frac_occs(Data,i1,i2,j1,j2,Ndim)

      implicit none 
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'frac_occ_vars.h'
      include 'proto_events.h'
      include 'mpif.h'

      Integer i,j,k,l,i1,i2,j1,j2
      Integer Lda,Ldb,Ndim,Nocc,Ierr 
      integer company_comm, pst_get_company_comm,master,pst_get_master
      Logical Found 

      Double Precision Data(i1:i2,j1:j2)

      Lda = i2 -i1 + 1
      Ldb = j2 -j1 + 1

      If (Nbasis .Gt. Maxbfns) Then
         Write(6,"(a,i4,a)") " The current no. of basis functions",
     +                       Nbasis, " is greater than the maximum"
         Write(6,"(a,i4,a)") "allowed", Maxbfns, "functions"
         call abort_job()
      Endif

#ifdef _DEBUG_LVL0
      If (me .eq.0) Then
      Write(6,"(a,2I4)") "Variables, Nbasis, Nalpha_occupied:",
     +                    Nbasis,Nalpha_Occupied
      Endif
#endif 

#ifdef _DEBUG_LVL0
C As a debugging example test, lets consider H2O Dzp basis. The
C system has 
      Nocc = 0
      do i=1, Nalpha_occupied 
         Occ_numsa(i) = 1.0D0
         Nocc = Nocc + 1
      Enddo 

      Do j=Nocc+1, Nbasis 
         Occ_numsa(j) = 0.0D0
      Enddo 

      Nocc = 0
      do i=1, Nbeta_occupied 
         Occ_numsb(i) = 1.0D0
         Nocc = Nocc + 1
      Enddo 

      Do j=Nocc+1, Nbasis 
         Occ_numsb(j) = 0.0D0
      Enddo 
#endif 
C
C Use this call to load the fractional occupation numbers to the 
C Occ_nums varible. At the moment this is in a common block (and
C need to put on to heap.
     
      Master = Pst_get_master()

      If (Me .Eq. Master) Then
      Call Get_external_occnums(Occ_numsa,Occ_numsb,Nbasis,Found,
     +                          Me)
   
      If (.Not. Found) Then
CSSS         If (me .eq. 0) Then
         Write(6,"(2a)")  "               ************Warning*********",
     +                                   "*******"
         Write(6,"(3a)") " The fractional occupation number file",
     +                   " FRAC_OCCNUMS is not found. The canonical"
         Write(6,"(a)")  " integer occupations are ussed!"
         Write(6,"(2a)")  "               ************Warning*********",
     +                                   "*******"
CSS         Endif 

         Nocc = 0
         do i=1, Nalpha_occupied
            Occ_numsa(i) = 1.0D0
            Nocc = Nocc + 1
         Enddo

         Do j=Nocc+1, Nbasis
            Occ_numsa(j) = 0.0D0
         Enddo

         Nocc = 0
         do i=1, Nbeta_occupied
            Occ_numsb(i) = 1.0D0
            Nocc = Nocc + 1
         Enddo

         Do j=Nocc+1, Nbasis
            Occ_numsb(j) = 0.0D0
         Enddo
      Endif 
      Endif 
 
      Company_comm = Pst_get_company_comm(me)
      Call Mpi_bcast(Occ_numsa,Nalpha_occupied,Mpi_double_precision,
     +               Master,Company_comm,Ierr)
      Call Mpi_bcast(Occ_numsb,Nbeta_occupied,Mpi_double_precision,
     +               Master,Company_comm,Ierr)

C#ifdef _DEBUG_LVL0
      If (me .eq.0) then
      Write(6,"(a)") "The Occ_number arrays"
      
      Write(6,"(6F7.4)") (Occ_numsa(i),i=1, Nbasis)
      Write(6,"(6F7.4)") (Occ_numsb(i),i=1, Nbasis)
      endif 
C#endif 

      Return
      End

