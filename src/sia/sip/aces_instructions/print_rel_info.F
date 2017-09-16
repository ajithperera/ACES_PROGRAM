C  Copyright (c) 2003-2010 University of Florida

C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
      subroutine print_rel_info
     +                     (array_table, 
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)

c--------------------------------------------------------------------------
c Scalar value    Code    Printed  
c          -10    SCF     exit for no convergence 
c           -2    CCSD    exit for no convergence 
c
c            0    SCF     ca printed to file ca.data 
c           -1    SCF     cb printed to file cb.data 
c            1    SCF     iteration, energy, error 
c            2    SCF     final energy, error 
c            3    SCF     alpha orbital energies  
c            4    SCF     beta orbital energies  
c
c            5    RHF CCSD iteration, energy  
c            6    RHF CCSD/iguess iteration, energy  
c            7    RHF CCSD iteration, energy  
c
c            8    RHF lambda iteration pseudo energy 
c            9    RHF lambda/iguesd iteration pseudo energy 
c            10   RHF lambda iteration pseudo energy Summary  
c
c            11   SCF gradient 
c            12   MP2 gradient 
c            13   CCSD gradient after 1-p  
c            14   CCSD gradient final  
c
c            15      EOM-CCSD(T) RHF ENERGY
c            16      EOM-CCSD(T) UHF ENERGY; Singlets 
c            16,17   EOM-CCSD(T) UHF ENERGY; Triplets 
c
c            21   RHF T2 amplitude header 
c            22   RHF T1 amplitude header 
c            23   UHF T2ab amplitude header 
c            24   UHF T2aa amplitude header 
c            25   UHF T2bb amplitude header 
c            26   RHF T1a amplitude header 
c            27   RHF T1b amplitude header 
c
c            28   RHF CIS energies 
c            29   UHF CIS energies 
C
c            30   RHF CIS Singles amplitude header 
c            31   UHF CIS Singles amplitude header 
c
c            28   RHF CIS energies 
c            32   RHF LCCD iteration, energy  
c            33   RHF LCCD iteration, energy  
C
c            70   RHF EOM-CCSD energies (singlets)
c            71   UHF EOM-CCSD energies (triplets)
c           911   EOM warnings 
c   
c--------------------------------------------------------------------------
C Changes to make GNU complier happy, C. Taylor, 08/2013.
C
      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer narray_table, nindex_table, nsegment_table, 
     *        nblock_map_table 
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer type(mx_array_index) 
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table)

      integer*8 iarray, ievec, get_index_from_base, daddr
      integer i, j, n1, n2, n3, n4  
      integer array, evec_array
      integer array_type, evec_type
      integer nindex, nindex_evec
      integer sind, ind(mx_array_index)
      integer seg, index(mx_array_index), val1(mx_array_index), 
     *        val2(mx_array_index), junk

      integer block, blkndx, st_flag
      integer find_current_block
      integer*8 get_block_index
      integer stack

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

      if (dryrun) return
#ifdef ALTIX
      dptr = dshptr
#endif

       
      open(66,file='summary.out') 

c----------------------------------------------------------------------------
c   Locate the scalar. --> what do we want to print??   
c---------------------------------------------------------------------------

      evec_array  = op(c_result_array)
      evec_type   = array_table(c_array_type, evec_array)
      nindex_evec = array_table(c_nindex, evec_array)

      if (evec_type .ne. scalar_value) return

      if (evec_array .lt. 1 .or. evec_array .gt. narray_table) then
         print *,'Error: Invalid scalar in print_rel_info, line ',
     *     current_line
         print *,'Array index is ',array,' Allowable values are ',
     *      ' 1 through ',narray_table
         call abort_job()
      endif 

      sind =  array_table(c_scalar_index, evec_array)
      if (sind .lt. 1 .or. sind .gt. nscalar_table) then
         print *,'Scalar table index out of range in print_rel_info, ',
     *           'line ',current_line
         print *,'Index for array ',evec_array,' is ',ind,' should be ',
     *           'between 1 and ',nscalar_table
         call abort_job()
      endif

c----------------------------------------------------------------------------
c   Exiting for no SCF convergence  
c---------------------------------------------------------------------------

      if (dabs(scalar_table(sind)+10.0) .le. 0.001) then 
         if (me .eq. 0) then 
         write(66,10) 
         call flush(66) 
         call abort_job() 
         endif 
      endif 

c----------------------------------------------------------------------------
c   SCF Print Flag = 0 or -1 Print Ca/Cb to file   
c---------------------------------------------------------------------------

      if ((dabs(scalar_table(sind)) .le. 0.00001) .or.
     *    (dabs(scalar_table(sind)+1.0)) .le. 0.00001) then 

c----------------------------------------------------------------------------
c   Locate the data for the input matrix which is to be analyses.
c---------------------------------------------------------------------------

      array      = op(c_op1_array)

      array_type = array_table(c_array_type, array)
      nindex     = array_table(c_nindex, array)

      daddr = address_table(array)
      iarray  = get_index_from_base(daddr, x, 2)

      if (array .eq. 0) then
         print *,'Error: The first argument in print_rel_info 
     *            must be an array.'  
         print *,(op(i),i=1,loptable_entry)
         call abort_job()
      endif

      do i = 1, nindex
         ind(i) = array_table(c_index_array1+i-1,array)
      enddo

      n1 = index_table(c_index_size, ind(1))  ! pick up length of index
      n2 = index_table(c_index_size, ind(2))  ! pick up length of index

          call print_rhf_orbitals(x(iarray),scalar_table(sind),
     *             1,n1,1,n2)
      return 
      endif 


c----------------------------------------------------------------------------
c   SCF Print Flag = 1, 2, 3, 4     
c---------------------------------------------------------------------------

      if (scalar_table(sind) .le. 4.0) then 

c----------------------------------------------------------------------------
c   Locate the data for the input matrix which is to be analyses.
c---------------------------------------------------------------------------

      array      = op(c_op1_array)

      array_type = array_table(c_array_type, array)
      nindex     = array_table(c_nindex, array)

      daddr = address_table(array)
      iarray  = get_index_from_base(daddr, x, 2)

      if (array .eq. 0) then
         print *,'Error: The first argument in print_rel_info
     *            must be an array.'  
         print *,(op(i),i=1,loptable_entry)
         call abort_job()
      endif

      do i = 1, nindex
         ind(i) = array_table(c_index_array1+i-1,array)
      enddo

      n1 = index_table(c_index_size, ind(1))  ! pick up length of index
      n2 = index_table(c_index_size, ind(2))  ! pick up length of index

          call print_scf_cinfo(x(iarray),scalar_table(sind),
     *             1,n1,1,n2)
      endif 

c----------------------------------------------------------------------------
c   CCSD/lambda Print Flag = 5, 6, 7, 8, 9, 10   
c---------------------------------------------------------------------------

      if ((scalar_table(sind) .ge. 5.0) .and. 
     *    (scalar_table(sind) .le. 100.0)) then 

c----------------------------------------------------------------------------
c   Locate the data for the input matrix which is to be analyses.
c---------------------------------------------------------------------------

      array      = op(c_op1_array)

      array_type = array_table(c_array_type, array)
      nindex     = array_table(c_nindex, array)

      daddr = address_table(array)
      iarray  = get_index_from_base(daddr, x, 2)

      if (array .eq. 0) then
         print *,'Error: The first argument in return_sval  
     *            must be an array.'  
         print *,(op(i),i=1,loptable_entry)
         call abort_job()
      endif

      do i = 1, nindex
         ind(i) = array_table(c_index_array1+i-1,array)
      enddo

c---------------------------------------------------------------------------
c   Look up its address.
c---------------------------------------------------------------------------

c     block  = find_current_block(array, array_table(1,array),
c    *                            index_table, nindex_table,
c    *                            segment_table, nsegment_table,
c    *                            block_map_table, blkndx)
c     stack = array_table(c_array_stack,array) 
c     iarray = get_block_index(array, block, stack, blkndx, x, .true.)  

c---------------------------------------------------------------------------
c   Find the indices of the array block.
c---------------------------------------------------------------------------

c     do i = 1, 1 !nindex
c        index(i) = array_table(c_index_array1+i-1,array)
c        type(i)  = index_table(c_index_type, index(i))
c        seg      = index_table(c_current_seg,index(i))

c-------------------------------------------------------------------------
c   Get segment ranges.
c-------------------------------------------------------------------------

c        call get_index_segment(index(i), seg, segment_table,
c    *                             nsegment_table, index_table,
c    *                             nindex_table, val1(i), val2(i))
c     enddo

      n1 = index_table(c_index_size, ind(1))  ! pick up length of index
      n2 = index_table(c_index_size, ind(2))  ! pick up length of index

      if (scalar_table(sind) .ge. 5 .and. scalar_table(sind) .le. 7) 
     *    then  
          call print_ccsd_cinfo(x(iarray),scalar_table(sind),
     *             1,n1,1,n2)
      endif 

      if (scalar_table(sind) .ge. 8 .and. scalar_table(sind) .le. 10) 
     *    then  
          call print_lambda_cinfo(x(iarray),scalar_table(sind),
     *             1,n1,1,n2)
      endif 

      if (scalar_table(sind) .eq. 15.0) then  
          call print_ccsdpt_rhf(x(iarray),scalar_table(sind),
     *             1,n1,1,n2)
      endif 

      if (scalar_table(sind) .eq. 16.0 .OR. 
     *    scalar_table(sind) .eq. 17.0) then  
          call print_ccsdpt_uhf(x(iarray),scalar_table(sind),
     *             1,n1,1,n2)
      endif 

      endif 

c................Prakash Start...........................................
      if ((scalar_table(sind) .ge. 700.0) .and.
     *    (scalar_table(sind) .le. 720.0)) then

c----------------------------------------------------------------------------
c   Locate the data for the input matrix which is to be analyses.
c---------------------------------------------------------------------------

      array      = op(c_op1_array)

      array_type = array_table(c_array_type, array)
      nindex     = array_table(c_nindex, array)

      daddr = address_table(array)
      iarray  = get_index_from_base(daddr, x, 2)

      if (array .eq. 0) then
         print *,'Error: The first argument in return_sval  
     *            must be an array.'
         print *,(op(i),i=1,loptable_entry)
         call abort_job()
      endif

      do i = 1, nindex
         ind(i) = array_table(c_index_array1+i-1,array)
      enddo

      n1 = index_table(c_index_size, ind(1))  ! pick up length of index
      n2 = index_table(c_index_size, ind(2))  ! pick up length of index

      if (scalar_table(sind) .ge. 700 .and. scalar_table(sind) .le. 710)
     *   then
          call print_ccsdTx_cinfo(x(iarray),scalar_table(sind),
     *             1,n1,1,n2)
       endif
c Note that above is not used in printing anything?

       if (scalar_table(sind) .eq. 711) then
          call print_pol_info(x(iarray),scalar_table(sind),
     *             1,n1,1,n2)
       endif

C Polarizablity

       if (scalar_table(sind) .eq. 712) then
          call print_pol_info(x(iarray),scalar_table(sind),
     *             1,n1,1,n2)
       endif

C g-tensor: mass correction

       if (scalar_table(sind) .eq. 713) then
          call print_gten_mc_info(x(iarray),scalar_table(sind),
     *             1,n1,1,n2)
       endif
C g-tensor: 1-el gauge correction

       if (scalar_table(sind) .eq. 714) then
          call print_gten_1egc_info(x(iarray),scalar_table(sind),
     *             1,n1,1,n2)
        endif

C g-tensor: 1-el spin-orbit

       if (scalar_table(sind) .eq. 715) then
          call print_gten_1bl_info(x(iarray),scalar_table(sind),
     *             1,n1,1,n2)
        endif

C g-tensor: 2-el spin-orbit

       if (scalar_table(sind) .eq. 716) then
          call print_gten_2bl_info(x(iarray),scalar_table(sind),
     *             1,n1,1,n2)
        endif

C g-tensor: 2-el gauge correction

       if (scalar_table(sind) .eq. 717) then
          call print_gten_2egc_info(x(iarray),scalar_table(sind),
     *             1,n1,1,n2)
        endif

       if (scalar_table(sind) .eq. 718) then
          call print_dipole_info(x(iarray),scalar_table(sind),
     *             1,n1,1,n2)
        endif

      endif
c................... END_prakash......................................       

c----------------------------------------------------------------------------
c   SCF gradient Print Flag = 11     
c---------------------------------------------------------------------------

c     if (scalar_table(sind) .eq. 11.0) then 
c         call print_scfgrad_cinfo(scalar_table(sind))
c     endif 

c----------------------------------------------------------------------------
c   CCSD one gradient Print Flag = 13     
c---------------------------------------------------------------------------

c     if (scalar_table(sind) .eq. 13.0) then 
c         call print_ccsdonegrad_cinfo(scalar_table(sind))
c     endif 

c----------------------------------------------------------------------------
c   CCSD final gradient Print Flag = 14     
c---------------------------------------------------------------------------

c     if (scalar_table(sind) .eq. 14.0) then 
c         call print_ccsdgrad_cinfo(scalar_table(sind))
c     endif 

c----------------------------------------------------------------------------
c   CCSD Maximum T2 amplitudes header Print Flag = 21     
c---------------------------------------------------------------------------

      if (scalar_table(sind) .eq. 21.0) then 
          call print_rhfamp_cinfo(scalar_table(sind))
      endif 

c----------------------------------------------------------------------------
c   CCSD Maximum T1 amplitudes header Print Flag = 22     
c---------------------------------------------------------------------------

      if (scalar_table(sind) .eq. 22.0) then 
          call print_rhfamp_cinfo(scalar_table(sind))
      endif 

c----------------------------------------------------------------------------
c   CCSD Maximum UHF  amplitudes header Print Flag = 23-27      
c---------------------------------------------------------------------------

      if (scalar_table(sind) .ge. 23.0 .and. 
     *    scalar_table(sind) .le. 27.0) then 
          call print_rhfamp_cinfo(scalar_table(sind))
      endif 

c----------------------------------------------------------------------------
c   Exiting for no CCSD convergence  
c---------------------------------------------------------------------------

      if (dabs(scalar_table(sind)+2.0) .le. 0.001) then 
         if (me .eq. 0) then 
         write(66,20) 
         call flush(66) 
         call abort_job() 
         endif 
      endif 
c----------------------------------------------------------------------------
c   Exiting for no Brueckner CCSD convergence
c---------------------------------------------------------------------------

      if (dabs(scalar_table(sind)+20.0) .le. 0.001) then
         if (me .eq. 0) then
         write(66,30)
         call flush(66)
         call abort_job()
         endif
      endif

10    FORMAT(5x,'Maximum iteration reached and The SCF has not 
     *           converged') 
20    FORMAT(5x,'Maximum iteration reached and The CCSD has not 
     *           converged') 
30    FORMAT(5x,'Maximum iteration reached and The Brueckner CCSD  
     *           has not converged') 

c----------------------------------------------------------------------------
c   CIS excitation energy(RHF,UHF) Print Flag = 28,29
c---------------------------------------------------------------------------
      st_flag = scalar_table(sind)
      if (st_flag .eq. 28.0 .or. st_flag .eq. 29.0) then 
          call print_cis_cinfo(x(iarray),scalar_table(sind),
     *             1,n1,1,n2)
      endif 
c----------------------------------------------------------------------------
c   EOM-CCSD excitation energy(RHF,UHF) Print Flag = 70,71 (Right)
c                                                    80,81 (Left)
c---------------------------------------------------------------------------
      st_flag = scalar_table(sind)
      if (st_flag .eq. 70.0 .or. st_flag .eq. 71.0 .or. 
     *    st_flag .eq. 80.0 .or. st_flag .eq. 81.0) then
          call print_eomccsd_cinfo(x(iarray),scalar_table(sind),
     *             1,n1,1,n2)
      endif

c----------------------------------------------------------------------------
c   HF stability analysis flag=55;singlet 56;triplet 
c----------------------------------------------------------------------------
      st_flag = scalar_table(sind)
      If (st_flag .eq. 55 .or. st_flag .eq. 56) then 
          call print_hfstability_matrix(x(iarray),st_flag,
     *             1,n1,1,n2,)
      endif 


c----------------------------------------------------------------------------
c   CIS singles amplitude header *RHF) Print Flag = 30   
c---------------------------------------------------------------------------

      if (scalar_table(sind) .ge. 100.0 .and. 
     *    scalar_table(sind) .le. 600.0) then 
          call print_cis_rhf_sheader(scalar_table(sind))
      endif 

c----------------------------------------------------------------------------
c   LCCD RHF Energy Flag = 32
c---------------------------------------------------------------------------

      if (scalar_table(sind) .ge. 32.0 .and. scalar_table(sind) .le.
     *    33.0) then  
          call print_lccd_cinfo(x(iarray),scalar_table(sind),
     *             1,n1,1,n2)
      endif 

c----------------------------------------------------------------------------
c   DEA total energies (singlet and triplet) Flag = 40, 41
c---------------------------------------------------------------------------

     
      if (scalar_table(sind) .eq. 40.0 .or. 
     *    scalar_table(sind) .eq. 41.0) then 
          call print_dea_einfo(x(iarray),scalar_table(sind),
     *             1,n1,1,n2)
      endif 

c----------------------------------------------------------------------------
c   DIP total energies (singlet and triplet) Flag = 50, 51
c---------------------------------------------------------------------------

     
      if (scalar_table(sind) .eq. 50.0 .or. 
     *    scalar_table(sind) .eq. 51.0) then 
          call print_dip_einfo(x(iarray),scalar_table(sind),
     *             1,n1,1,n2)
      endif 
c----------------------------------------------------------------------------
c   IP/EA total energies (doublet) Flag =60, 61
c---------------------------------------------------------------------------

     
      if (scalar_table(sind) .eq. 60.0 .or. 
     *    scalar_table(sind) .eq. 61.0) then 
          call print_ip_ea_einfo(x(iarray),scalar_table(sind),
     *             1,n1,1,n2)
      endif 
c----------------------------------------------------------------------------
c   Warning flags (911:EOM)
c---------------------------------------------------------------------------
      if (scalar_table(sind) .eq. 911.0) then 
          call print_warning(x(iarray),scalar_table(sind),
     *             1,n1,1,n2)
      endif 
c----------------------------------------------------------------------------

      call flush(66) 

      return 
      end 

      subroutine print_ccsd_cinfo(array1,sval,a1,a2,b1,b2) 

c---------------------------------------------------------------------------
c
c 
c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
c
      integer i, ival, amax, a, a1, a2, b, b1, b2   
      double precision array1(a1:a2,b1:b2), sval, stemp   
      double precision efinal, conv   

      open(66,file='summary.out') 

       ival = sval 
       conv = 27.21138505 
c      if (me .eq. my_company_rank)then

       if (me .eq. 0)then

       if (ival .eq. 6) then  
          do a = a1, a2   
            if (a .eq. 1.0) then 
               write(66,*) 
               write(66,*) 
               write(66,108) 
               write(66,101) 
               write(66,110) 
               write(66,101) 
            endif  
c           if (dabs(array1(a,b1)) .gt. 0.0001) then 
               write(66,99) a-1, (array1(a,i), i=b1, b2) 
               write(66,101) 
               write(66,*) 
               write(66,*) 
               return 
c           endif 
          enddo 
       endif 

       if (ival .eq. 5) then  
          do a = a2, a1, -1  
            if (a .eq. 1.0) then 
               write(66,*) 
               write(66,*) 
               write(66,98) 
               write(66,101) 
               write(66,100) 
               write(66,101) 
            endif  
            if (dabs(array1(a,b1)) .gt. 0.0001) then 
               write(66,99) a, (array1(a,i), i=b1, b2) 
               return 
            endif 
          enddo 
       endif 

       if (ival .eq. 7) then  
          do a = a1, a2 
            if (a .eq. 1.0) then 
               write(66,*) 
               write(66,*) 
               write(66,97) 
               write(66,101) 
               write(66,100) 
               write(66,101) 
            endif  
            if (dabs(array1(a,b1)) .lt. 0.0001) then 
c              write(66,99) a-1, (array1(a-1,i), i=b1, b2) 
               amax = a - 1  
               efinal = array1(a-1,2) 
            go to 10 
            endif 
            write(66,99) a, (array1(a,i), i=b1, b2) 
          enddo 
10        continue 
          write(66,101) 
          write(66,102) efinal, amax 
          write(66,101) 
          write(66,*) 
          write(66,*) 

       endif 

      endif

97    format(2x,"SUMMARY OF CCSD CALCULATION") 
98    format(2x,"START OF CCSD CALCULATION") 
108   format(2x,"SUMMARY OF MBPT(2) CALCULATION") 
99    format(2x,I8,5x,F20.10,2x,F20.10) 
100   format(2x,"Iteration#",9x,"CCSD corr energy",7x,"Total energy")
110   format(2x,"Iteration#",9x,"MBPT(2) corr energy",7x,"Total energy")
101   format(2x,60('-'))
102   format(2x,"Final CCSD energy = ",F20.10,2x,"After",2x,I8,2x,
     *       "iterations")

      return
      end


      subroutine print_scf_cinfo(array1,sval,a1,a2,b1,b2) 

c---------------------------------------------------------------------------
c
c 
c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
c
      integer i, ival, amax, a, a1, a2, b, b1, b2   
      double precision array1(a1:a2,b1:b2), sval, stemp   
      double precision efinal, conv   

      open(66,file='summary.out') 

       ival = sval 
       conv = 27.21138505 
c      if (me .eq. my_company_rank)then

       if (me .eq. 0)then

       if (ival .eq. 1) then  
          do a = a2, a1, -1  
            if (a .eq. 1.0) then 
               write(66,*) 
               write(66,*) 
               write(66,98) 
               write(66,101) 
               write(66,100) 
               write(66,101) 
            endif  
            if (dabs(array1(a,b1)) .gt. 0.0001) then 
               write(66,99) a, (array1(a,i), i=b1, b2) 
               return 
            endif 
          enddo 
       endif 

       if (ival .eq. 2) then  
          do a = a1, a2 
            if (a .eq. 1.0) then 
               write(66,*) 
               write(66,*) 
               write(66,97) 
               write(66,101) 
               write(66,100) 
               write(66,101) 
            endif  
            if (dabs(array1(a,b1)) .lt. 0.0001) then 
c              write(66,99) a-1, (array1(a-1,i), i=b1, b2) 
               amax = a - 1  
               efinal = array1(a-1,1) 
            go to 10 
            endif 
            write(66,99) a, (array1(a,i), i=b1, b2) 
          enddo 
10        continue 
          write(66,101) 
          write(66,102) efinal, amax 
          write(66,101) 
          write(66,*) 
          write(66,*) 

       endif 

       if (ival .eq. 3) then  
          do a = a1, nalpha_occupied   
            if (a .eq. 1.0) then 
               write(66,*) 
               write(66,*) 
               write(66,96) 
               write(66,101) 
               write(66,103) 
               write(66,101) 
            endif  
            write(66,105) a, array1(a,a), array1(a,a)*conv   
          enddo  
               write(66,101) 
          do a = nalpha_occupied +1, a2  
               write(66,105) a, array1(a,a), array1(a,a)*conv   
          enddo 
               write(66,101) 
       endif 

       if (ival .eq. 4) then  
          do a = a1, nbeta_occupied  
            if (a .eq. 1.0) then 
               write(66,*) 
               write(66,*) 
               write(66,96) 
               write(66,101) 
               write(66,104) 
               write(66,101) 
            endif  
            write(66,105) a, array1(a,a), array1(a,a)*conv   
          enddo  
               write(66,101) 
          do a = nbeta_occupied +1, a2  
               write(66,105) a, array1(a,a), array1(a,a)*conv   
          enddo 
               write(66,101) 
       endif 

      endif

96    format(2x,"Orbital energies") 
97    format(2x,"SUMMARY OF SCF CALCULATION") 
98    format(2x,"START OF SCF CALCULATION") 
99    format(2x,I8,5x,F20.10,2x,F12.8) 
100   format(2x,"Iteration#",9x,"SCF energy",7x,"Density difference")
101   format(2x,60('-'))
102   format(2x,"Final SCF energy = ",F20.10,2x,"After",2x,I8,2x,
     *       "iterations")
103   format(2x,"Alpha orbital energies",4x,"Hartree",9x,"eV")
104   format(2x,"Beta orbital energies",4x,"Hartree",9x,"eV")
105   format(2x,I8,5x,F20.10,2x,F20.10) 

      return
      end


      subroutine print_lambda_cinfo(array1,sval,a1,a2,b1,b2) 

c---------------------------------------------------------------------------
c
c 
c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
c
      integer i, ival, amax, a, a1, a2, b, b1, b2   
      double precision array1(a1:a2,b1:b2), sval, stemp   
      double precision efinal, conv   

      open(66,file='summary.out') 

       ival = sval 
       conv = 27.21138505 
c      if (me .eq. my_company_rank)then

       if (me .eq. 0)then

       if (ival .eq. 9) then  
          do a = a1, a2   
            if (a .eq. 1.0) then 
               write(66,*) 
               write(66,*) 
               write(66,108) 
               write(66,95) 
               write(66,101) 
               write(66,110) 
               write(66,101) 
            endif  
c           if (dabs(array1(a,b1)) .gt. 0.0001) then 
               write(66,99) a-1, (array1(a,i), i=b1, b2) 
               write(66,101) 
               write(66,*) 
               write(66,*) 
               return 
c           endif 
          enddo 
       endif 

       if (ival .eq. 8) then  
          do a = a2, a1, -1  
            if (a .eq. 1.0) then 
               write(66,*) 
               write(66,*) 
               write(66,98) 
               write(66,95) 
               write(66,101) 
               write(66,100) 
               write(66,101) 
            endif  
            if (dabs(array1(a,b1)) .gt. 0.0001) then 
               write(66,99) a, (array1(a,i), i=b1, b2) 
               return 
            endif 
          enddo 
       endif 

       if (ival .eq. 10) then  
          do a = a1, a2 
            if (a .eq. 1.0) then 
               write(66,*) 
               write(66,*) 
               write(66,97) 
               write(66,95) 
               write(66,101) 
               write(66,100) 
               write(66,101) 
            endif  
            if (dabs(array1(a,b1)) .lt. 0.0001) then 
c              write(66,99) a-1, (array1(a-1,i), i=b1, b2) 
               amax = a - 1  
               efinal = array1(a-1,2) 
            go to 10 
            endif 
            write(66,99) a, (array1(a,i), i=b1, b2) 
          enddo 
10        continue 
          write(66,101) 
          write(66,102) efinal, amax 
          write(66,101) 
          write(66,*) 
          write(66,*) 

       endif 

      endif

95    format(2x,"Note that Energies are Pseudo energies")  
97    format(2x,"SUMMARY OF LAMBDA CALCULATION") 
98    format(2x,"START OF LAMBDA CALCULATION") 
108   format(2x,"SUMMARY OF LAMBDA(2) CALCULATION") 
99    format(2x,I8,5x,F20.10,2x,F20.10) 
100   format(2x,"Iteration#",9x,"lambda corr energy",7x,"Total energy")
110   format(2x,"Iteration#",9x,"lambda(2) corr energy",7x,
     *       "Total energy")
101   format(2x,60('-'))
102   format(2x,"Final lambda energy = ",F20.10,2x,"After",2x,I8,2x,
     *       "iterations")

      return
      end

      subroutine print_scfgrad_cinfo(sval)

c---------------------------------------------------------------------------
c
c 
c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
      include 'saved_data.h' 
c
      integer iatom, ival, component, igrad1, igrad2, igrad3  
      double precision sval

      open(66,file='summary.out') 

      ival = sval 

      if (me .eq. 0) then 

      write(66,*) ' ' 
      write(66,10) 
      write(66,11) 
      write(66,10) 
      write(66,12) 
      write(66,10) 
      do iatom = 1, ncenters

c-----------------------------------------------------------------------------
c   Print the gradient component  
c-----------------------------------------------------------------------------

            igrad1 = (iatom-1)*3 + 1
            igrad2 = (iatom-1)*3 + 2
            igrad3 = (iatom-1)*3 + 3
            write(66,100) iatom, gradient_data(igrad1),  
     *                   gradient_data(igrad2), gradient_data(igrad3) 

      enddo
      write(66,10) 
      write(66,*) ' ' 

10    format(2x,60('-'))
11    format(5x,"FINAL GRADIENT INFORMATION IN CARTESIAN COORDINATES")
12    format(10x,"Center",10x,"X",12x,"Y",14x,"Z")
100   Format(3x,I10,5x,3F14.10) 

      endif 

      return 
      end 

      subroutine print_ccsdonegrad_cinfo(sval)

c---------------------------------------------------------------------------
c
c 
c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
      include 'saved_data.h' 
c
      integer iatom, ival, component, igrad1, igrad2, igrad3  
      double precision sval

      open(66,file='summary.out') 

      ival = sval 

      if (me .eq. 0) then 

      write(66,*) ' ' 
      write(66,10) 
      write(66,11) 
      write(66,10) 
      write(66,12) 
      write(66,10) 
      do iatom = 1, ncenters

c-----------------------------------------------------------------------------
c   Print the gradient component  
c-----------------------------------------------------------------------------

            igrad1 = (iatom-1)*3 + 1
            igrad2 = (iatom-1)*3 + 2
            igrad3 = (iatom-1)*3 + 3
            write(66,100) iatom, gradient_data(igrad1),  
     *                   gradient_data(igrad2), gradient_data(igrad3) 

      enddo
      write(66,10) 
      write(66,*) ' ' 

10    format(2x,60('-'))
11    format(5x,"CCSD ONE GRADIENT IN CARTESIAN COORDINATES")
12    format(10x,"Center",10x,"X",12x,"Y",14x,"Z")
100   Format(3x,I10,5x,3F14.10) 

      endif 

      return 
      end 

      subroutine print_ccsdgrad_cinfo(sval)

c---------------------------------------------------------------------------
c
c 
c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
      include 'saved_data.h' 
c
      integer iatom, ival, component, igrad1, igrad2, igrad3  
      double precision sval

      open(66,file='summary.out') 

      ival = sval 

      if (me .eq. 0) then 

      write(66,*) ' ' 
      write(66,10) 
      write(66,11) 
      write(66,10) 
      write(66,12) 
      write(66,10) 
      do iatom = 1, ncenters

c-----------------------------------------------------------------------------
c   Print the gradient component  
c-----------------------------------------------------------------------------

            igrad1 = (iatom-1)*3 + 1
            igrad2 = (iatom-1)*3 + 2
            igrad3 = (iatom-1)*3 + 3
            write(66,100) iatom, gradient_data(igrad1),  
     *                   gradient_data(igrad2), gradient_data(igrad3) 

      enddo
      write(66,10) 
      write(66,*) ' ' 

10    format(2x,60('-'))
11    format(5x,"FINAL CCSD GRADIENT IN CARTESIAN COORDINATES")
12    format(10x,"Center",10x,"X",12x,"Y",14x,"Z")
100   Format(3x,I10,5x,3F14.10) 

      endif 

      return 
      end 
c---------------------------------------------------------------------------

      subroutine print_ccsdpt_rhf(array1,sval,a1,a2,b1,b2) 

c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
c
      integer i, ival, amax, a, a1, a2, b, b1, b2, ieom_roots
      double precision array1(a1:a2,b1:b2), sval, stemp   
      double precision Eccsd,Eccsdt,Etotal,Conv
      conv = 27.21138505

      open(66,file='summary.out') 

      ieom_roots = eom_roots 

      If (me .eq. 0) then 
         write(66,10)
         write(66,100)
         write(66,30)
         write(66,100)

          do i = 1, ieom_roots

CSSS#ifdef _DEBUG_LVL0
        Write(6,"(a,2i)")     "Root   =", i
        Write(6,"(a,F15.10)") "E-CCSD =", EE_CCSD(i)
        Write(6,"(a,F15.10)") "E4_aaa =", E4_aaa(i)
        Write(6,"(a,F15.10)") "E5_aaa =", E5_aaa(i)
        Write(6,"(a,F15.10)") "E4_aab =", E4_aab(i)
        Write(6,"(a,F15.10)") "E5_aab =", E5_aab(i)
        Write(6,"(a,F15.10)") "Total  =", E5_aab(i)+E4_aab(i)+
     &                         E5_aaa(i)+E4_aaa(i)
CSSS#endif 
             Eccsd  = EE_CCSD(i) 
             Eccsdt = E4_aaa(i)+E5_aaa(i)+E4_aab(i)+E5_aab(i)
             Etotal = Eccsd + Eccsdt
             write(66,20) i, Eccsd*Conv, Etotal*Conv

          enddo
          write(66,100)

      endif

10    format(2x,"EOM-CCSD(T) singlet excitation energies")  
100   format(2x,60('-'))
20    format(2x,I4,x,F20.8,x,F20.8)
30    format(2x,'ROOT',12x,'EOM-CCSD ',12x,'EOM-CCSD(T) (eV)')

      return
      end

c---------------------------------------------------------------------------

      subroutine print_ccsdpt_uhf(array1,sval,a1,a2,b1,b2) 

c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
c
      integer i, ival, amax, a, a1, a2, b, b1, b2,ieom_roots 
      double precision array1(a1:a2,b1:b2), sval, stemp   
      double precision efinal,Two
      double precision Eccsd,Eccsdt,Etotal,Conv

      conv = 27.21138505
      Two  = 2.0D0

      ieom_roots = eom_roots

      If (me .eq. 0) then
         If (Sval .Eq. 16) Write(66,10)
         If (Sval .Eq. 17) Write(66,40)

         write(66,100)
         write(66,30)
         write(66,100)

          do i = 1, ieom_roots

CSSS#ifdef _DEBUG_LVL0
        If (sval .eq. 16) Then
           Write(6,"(a,2i)")     "Root   =", i
           Write(6,"(a,F15.10)") "E-CCSD =", EE_CCSD(i)
           Write(6,"(a,F15.10)") "E4_aaa =", E4_aaa(i)
           Write(6,"(a,F15.10)") "E5_aaa =", E5_aaa(i)
           Write(6,"(a,F15.10)") "E4_aab =", E4_aab(i)
           Write(6,"(a,F15.10)") "E5_aab =", E5_aab(i)
           Write(6,"(a,F15.10)") "E5_bbb =", E4_bbb(i)
           Write(6,"(a,F15.10)") "E4_bbb =", E5_bbb(i)
           Write(6,"(a,F15.10)") "E5_bba =", E4_bba(i)
           Write(6,"(a,F15.10)") "E4_bba =", E5_bba(i)

           Write(6,"(a,F15.10)")"Total  =", (E5_aab(i)+E4_aab(i)+
     &                            E5_aaa(i)+E4_aaa(i))*Two
        else
           Write(6,"(a,2i)")     "Root   =", i
           Write(6,"(a,F15.10)") "E-CCSD =", EE_CCSD(i)
           Write(6,"(a,F15.10)") "E4_aaa =", E4_aaa(30+i)
           Write(6,"(a,F15.10)") "E5_aaa =", E5_aaa(30+i)
           Write(6,"(a,F15.10)") "E4_aab =", E4_aab(30+i)
           Write(6,"(a,F15.10)") "E5_aab =", E5_aab(30+i)
           Write(6,"(a,F15.10)") "E5_bbb =", E4_bbb(30+i)
           Write(6,"(a,F15.10)") "E4_bbb =", E5_bbb(30+i)
           Write(6,"(a,F15.10)") "E5_bba =", E4_bba(30+i)
           Write(6,"(a,F15.10)") "E4_bba =", E5_bba(30+i)

           Write(6,"(a,F15.10)")"Total  =", E5_aab(30+i)+
     &                                      E4_aab(30+i)+
     &                                      E5_aaa(30+i)+
     &                                      E4_aaa(30+i)+
     &                                      E4_bbb(30+i)+
     &                                      E5_bbb(30+i)+
     &                                      E4_bba(30+i)+
     &                                      E5_bba(30+i)
        endif 
CSSS#endif
            If (Sval .EQ. 16) Then
               Eccsd  = EE_CCSD(i)
               Eccsdt = E4_aaa(i)+E5_aaa(i)+E4_aab(i)+E5_aab(i)
               Eccsdt = Eccsdt * Two
               Etotal = Eccsd + Eccsdt
               write(66,20) i, Eccsd*Conv, Etotal*Conv
            Endif 

            If (Sval .EQ. 17) Then
               Eccsd  = EE_CCSD(30+i)
               Eccsdt = E4_aaa(30+i)+E5_aaa(30+i)+
     &                  E4_aab(30+i)+E5_aab(30+i)+
     &                  E4_bbb(30+i)+E5_bbb(30+i)+
     &                  E4_bba(30+i)+E5_bba(30+i)
               Etotal = Eccsd + Eccsdt
               write(66,20) i, Eccsd*Conv, Etotal*Conv
            Endif 

            enddo
            write(66,100)

      endif

10    format(2x,"EOM-CCSD(T) singlet excitation energies")  
40    format(2x,"EOM-CCSD(T) triplet excitation energies")  
100   format(2x,60('-'))
20    format(2x,I4,x,F20.8,x,F20.8)
30    format(2x,'ROOT',12x,'EOM-CCSD ',12x,'EOM-CCSD(T) (eV)')

      return
      end

c---------------------------------------------------------------------------
      subroutine print_rhfamp_cinfo(sval)

c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
      include 'saved_data.h' 
c
      integer iatom, ival, component, igrad1, igrad2, igrad3  
      double precision sval

      open(66,file='summary.out') 

      ival = sval 

      if (me .eq. 0)then

c-----------------------------------------------------------------------------
c   Print the header  
c-----------------------------------------------------------------------------


      write(66,*) ' ' 
      write(66,10) 
      write(66,10) 
c RHF T2ab 
      if (ival .eq. 21) then 
         write(66,11) 
         write(66,10) 
         write(66,12) 
      endif 
c RHF T1a 
      if (ival .eq. 22 .or. ival .eq. 26) then 
         write(66,13) 
         write(66,10) 
         write(66,14) 
      endif 
c UHF T2ab 
      if (ival .eq. 25) then 
         write(66,11) 
         write(66,10) 
         write(66,12) 
      endif 
c UHF T2aa 
      if (ival .eq. 23) then 
         write(66,11) 
         write(66,10) 
         write(66,15) 
      endif 
c UHF T2bb 
      if (ival .eq. 24) then 
         write(66,11) 
         write(66,10) 
         write(66,16) 
      endif 
c UHF T1b 
      if (ival .eq. 27) then 
         write(66,13) 
         write(66,10) 
         write(66,17) 
      endif 
      write(66,10) 

      endif 

10    format(2x,60('-'))
11    format(5x,"LARGEST      T2      AMPLITUDES")
12    format(5x,'A',8x,'I',8x,'B',8x,'J',8x,'T2(A,I,B,J)')
13    format(5x,"LARGEST      T1      AMPLITUDES")
14    format(5x,'A',8x,'I',8x,'T1(A,I)')
15    format(5x,'A',8x,'I',8x,'A',8x,'I',8x,'T2(A,I,A,I)')
16    format(5x,'B',8x,'J',8x,'B',8x,'J',8x,'T2(B,J,B,J)')
17    format(5x,'B',8x,'B',8x,'T1(B,J)')
100   Format(3x,I10,5x,3F14.10) 

      return 
      end 

      subroutine print_cis_cinfo(array1,sval,a1,a2,b1,b2) 

c---------------------------------------------------------------------------
c
c 
c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
c
      integer i, ival, amax, a, a1, a2, b, b1, b2   
      integer ieom_roots 
      double precision array1(a1:a2,b1:b2), sval, stemp   
      double precision efinal, conv   
      logical singlet, triplet 

      open(66,file='summary.out') 

      ival = sval 
      conv = 27.21138505 
      ieom_roots = eom_roots 
      amax = max(ieom_roots,6) 
      singlet = (Sval .eq. 28)
      triplet = (Sval .eq. 29)

      if (me .eq. 0)then

         If (singlet) then
            write(66,10) 
            write(66,100) 
            write(66,30) 
            write(66,100) 
            do i = 1, ieom_roots
               write(66,20) i, array1(i,1), array1(i,1)*conv
            enddo 
            write(66,100) 
            Call Dputrec(20,"JOBARC","CIS_EXES",ieom_roots,Array1)

         Elseif (Triplet) then

            write(66,*)
            write(66,40) 
            write(66,100) 
            write(66,30) 
            write(66,100) 
            do i = 1, ieom_roots
               write(66,20) i, array1(i,1), array1(i,1)*conv
            enddo 
            Call Dputrec(20,"JOBARC","CIS_EXES",ieom_roots,Array1)
            write(66,100) 

         Endif 

      endif 


100   format(2x,60('-'))
10    format(2x,"CIS singlet excitation energies")  
40    format(2x,"CIS Triplet excitation energies")  
20    format(2x,I4,5x,F20.8,x,F20.8)  
30    format(2x,'ROOT',15x,'Hartree',15x,'eV')  
C 

      return
      end

c---------------------------------------------------------------------------
      subroutine print_eomccsd_cinfo(array1,sval,a1,a2,b1,b2) 

c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
c
      integer i, ival, amax, a, a1, a2, b, b1, b2   
      integer ieom_roots 
      double precision array1(a1:a2,b1:b2), sval, stemp   
      double precision efinal, conv   
      logical singlet, triplet, left_hand, right_hand 

      open(66,file='summary.out') 

      ival = sval 
      conv = 27.21138505 
      ieom_roots = eom_roots 
      amax = max(ieom_roots,6) 

      singlet = .False.
      triplet = .False.
      Right_hand = .False.
      Left_hand  = .False.
      singlet = (Sval .eq. 70 .or. Sval .eq. 80)
      triplet = (Sval .eq. 71 .or. Sval .eq. 81)

      Right_hand = (Sval .eq. 70 .or. Sval .eq. 71)
      Left_hand  = (Sval .eq. 80 .or. Sval .eq. 81)

      if (me .eq. 0)then

         If (singlet) then
            If (Right_hand) then
                write(66,10) 
            Elseif (Left_hand) then
                write(66,11) 
            Endif 
            write(66,100) 
            write(66,30) 
            write(66,100) 
            do i = 1, ieom_roots
               write(66,20) i, array1(i,1), array1(i,1)*conv
            enddo 
            write(66,100) 
            Call Dputrec(20,"JOBARC","EOM_EXES",ieom_roots,Array1(1,1))

         Elseif (Triplet) then

            write(66,*)
            If (Right_hand) then
               write(66,40) 
            Elseif (Left_hand) then
               write(66,41) 
            Endif
            write(66,100) 
            write(66,30) 
            write(66,100) 
            do i = 1, ieom_roots
               write(66,20) i, array1(i,1), array1(i,1)*conv
            enddo 
            Call Dputrec(20,"JOBARC","EOM_EXES",ieom_roots,Array1(1,1))
            write(66,100) 

         Endif 

      endif 

100   format(2x,60('-'))
10    format(2x,"Right EOM-CCSD singlet excitation energies")  
11    format(2x,"Left  EOM-CCSD singlet excitation energies")  
40    format(2x,"Right EOM-CCSD Triplet excitation energies")  
41    format(2x,"Left  EOM-CCSD Triplet excitation energies")  
20    format(2x,I4,5x,F20.8,x,F20.8)  
30    format(2x,'ROOT',15x,'Hartree',15x,'eV')  
C 
      return
      end

c---------------------------------------------------------------------------
      subroutine print_warning (array1,sval,a1,a2,b1,b2)
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
c
      integer a1, a2, b, b1, b2
      integer ieom_roots
      double precision array1(a1:a2,b1:b2), sval

      If (me   .eq. 0) then
      If (Sval .eq.911) then

      Write(66,*)
      Write(66,"(2x,a)") "Warnings!!"
      Write(66,"(2x,a)") "----------"
      Write(66,"(2x,a)") "Zero entry: Did not converged!"
      Write(66,"(2x,a)") "Left and right energies should be identical."
      Write(66,"(2x,a,a)") "Right converged Left did not ",
     &                     "or vice versa: rare but possible."
      Write(66,"(2x,a,a,a)") "Possible solutions:"
      Write(66,"(2x,a)")"1. Increase the number of iterations (<50)."
      Write(66,"(2x,a)")"2. Change the denominator from Fock to Hbar."

      Endif
      Endif

      return
      end

c---------------------------------------------------------------------------
      subroutine print_hfstability_matrix(array1,sval,a1,a2,b1,b2)
c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
c
      integer i, ival, amax, a, a1, a2, b, b1, b2
      integer ieom_roots,sval
      double precision array1(a1:a2,b1:b2), stemp
      double precision efinal, conv
      logical Singlet, Triplet 

      open(66,file='summary.out')

      conv = 27.21138505

      singlet = (Sval .eq. 55)
      triplet = (Sval .eq. 56)
      ieom_roots = eom_roots 
    
      if (me .eq. 0)then

        If (Singlet)  Then

         write(66,10)
         write(66,100)
         write(66,30)
         write(66,100)
         do i = 1, ieom_roots
            write(66,20) i, array1(i,1), array1(i,1)*conv
         enddo
         write(66,100)
 
        Elseif (Triplet) then

         Write(66,*)
         write(66,40)
         write(66,100)
         write(66,30)
         write(66,100)
         do i = 1, ieom_roots
            write(66,20) i, array1(i,1), array1(i,1)*conv
         enddo
         write(66,100)

        Endif

      endif

100   format(2x,60('-'))
10    format(2x,"HF-stability analysis: Singlet ---> Singlet")
40    format(2x,"HF-stability analysis: Singlet ---> Triplet")
20    format(2x,I4,5x,F20.8,x,F20.8)
30    format(2x,'ROOT',15x,'Hartree',15x,'eV')
C
      return
      end
c---------------------------------------------------------------------------


      subroutine print_dea_einfo(array1,sval,a1,a2,b1,b2)

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
c
      integer i, ival, amax, a, a1, a2, b, b1, b2
      integer ieom_roots
      double precision array1(a1:a2,b1:b2), sval, stemp
      double precision efinal, conv

      open(66,file='summary.out')

      ival = sval
      conv = 27.21138505
      ieom_roots = eom_roots

      if (me .eq. 0)then

         if (sval .eq. 40) write(66,10)
         if (sval .eq. 41) write(66,20)
         write(66,100)
         write(66,40)
         write(66,100)
         do i = 1, ieom_roots
            write(66,30) i, array1(i,1), array1(i,1)*conv 
         enddo
         write(66,100)
      endif


100   format(2x,60('-'))
10    format(2x,"DEA Singlet energies")
20    format(2x,"DEA Triplet energies")
30    format(2x,I4,5x,F20.8,x,F20.8)
40    format(2x,'ROOT',15x,'Hartree',15X,"eV")

      return
      end
c---------------------------------------------------------------------------

      subroutine print_dip_einfo(array1,sval,a1,a2,b1,b2)

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
c
      integer i, ival, amax, a, a1, a2, b, b1, b2
      integer ieom_roots
      double precision array1(a1:a2,b1:b2), sval, stemp
      double precision efinal, conv

      open(66,file='summary.out')

      ival = sval
      conv = 27.21138505
      ieom_roots = eom_roots

      if (me .eq. 0)then

         if (sval .eq. 50) write(66,10)
         if (sval .eq. 51) write(66,20)
         write(66,100)
         write(66,40)
         write(66,100)
         do i = 1, ieom_roots 
            write(66,30) i, array1(i,1), array1(i,1)*Conv 
         enddo
         write(66,100)
      endif


100   format(2x,60('-'))
10    format(2x,"DIP Singlet energies")
20    format(2x,"DIP Triplet energies")
30    format(2x,I4,5x,F20.8,x,F20.8)
40    format(2x,'ROOT',15x,'Hartree',15X,"eV")

      return
      end

c---------------------------------------------------------------------------
      subroutine print_ip_ea_einfo(array1,sval,a1,a2)
c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
c
      integer i, ival, amax, a, a1, a2, b, b1, b2
      integer ieom_roots
      double precision array1(a1:a2), sval, stemp
      double precision efinal, conv

      open(66,file='summary.out')

      ival = sval
      conv = 27.21138505
      ieom_roots = eom_roots

      if (me .eq. 0)then
         If (Sval .eq. 60) write(66,10)
         If (Sval .eq. 61) write(66,11)
         write(66,100)
         write(66,30)
         write(66,100)
         do i = 1, ieom_roots
            write(66,20) i, array1(i),  array1(i)*conv
         enddo
         write(66,100)
      endif

100   format(2x,60('-'))
10    format(2x,"IP-EOM-CCSD energies")
11    format(2x,"EA-EOM-CCSD energies")
20    format(2x,I4,5x,F20.8,x,F20.8)
30    format(2x,'ROOT',15x,'Hartree',15x,'eV')

      return
      end

c---------------------------------------------------------------------------
      subroutine print_eom_rhf_cinfo(array1,sval,a1,a2,b1,b2)
c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
c
      integer i, ival, amax, a, a1, a2, b, b1, b2
      integer ieom_roots
      double precision array1(a1:a2,b1:b2), sval, stemp
      double precision efinal, conv

      open(66,file='summary.out')

      ival = sval
      conv = 27.21138505
      ieom_roots = eom_roots
      amax = max(ieom_roots,6)

      if (me .eq. 0)then

         write(66,10)
         write(66,100)
         write(66,30)
         write(66,100)
         do i = 1, 12
            if (i .le. amax) then
            write(66,20) i, array1(i,1)/conv, array1(i,1)
            endif
         enddo
         Write(6,*) "Writing to JOBARC"
         Call Dputrec(20,"JOBARC","EOM_EXES",ieom_roots,Array1)

      endif

100   format(2x,60('-'))
10    format(2x,"EOM-CCSD excitation energies")
20    format(2x,I4,5x,F20.8,x,F20.8)
30    format(2x,'ROOT',15x,'Hartree',15x,'eV')

      return
      end

      subroutine print_cis_uhf_cinfo(array1,sval,a1,a2,b1,b2) 

c---------------------------------------------------------------------------
c
c 
c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
c
      integer i, ival, amax, a, a1, a2, b, b1, b2   
      double precision array1(a1:a2,b1:b2), sval, stemp   
      double precision efinal, conv   

      open(66,file='summary.out') 

      ival = sval 
      conv = 27.21138505 

      if (me .eq. 0)then

         write(66,10) 
         do i = 1, 12 
            write(66,20) i, array1(i,1) 
         enddo 

      endif 


10    format(2x,"CIS excitation energies")  
20    format(2x,I4,F10.8)  

      return
      end


      subroutine print_cis_rhf_sheader(sval) 

c---------------------------------------------------------------------------
c
c 
c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
c
      double precision sval, stemp   
      integer i, ival, jval  

      open(66,file='summary.out') 

      ival = sval 
      do i = 1, 5 
         ival = ival - 100 
         if (ival .lt. 0) ival = ival + 100 
      enddo 
             

10    FORMAT(2x,'Largest EOM amplitudes for root',x,I8) 
20    FORMAT(2x,60('-'))
30    FORMAT(2x,I8,x,I8,x,F16.8)
40    FORMAT(8x,'A',10x,'B',8x,'C1A(A,B)')
50    FORMAT(8x,'A',10x,'B',8x,'C1B(A,B)')
60    FORMAT(7x,'A',11x,'B',11x,'C',11x,'D',8x,'C2AA(A,B,C,D)')
70    FORMAT(7x,'A',11x,'B',11x,'C',11x,'D',8x,'C2BB(A,B,C,D)')
80    FORMAT(7x,'A',11x,'B',11x,'C',11x,'D',8x,'C2AB(A,B,C,D)')

       if (me .eq. 0)then

         write(66,*) ' '
         write(66,20)
         write(66,10) ival 
         write(66,20)
         if (sval .gt. 500.0) then 
            write(66,80)
            write(66,20)
            return 
         endif 
         if (sval .gt. 400.0) then 
            write(66,70)
            write(66,20)
            return 
         endif 
         if (sval .gt. 300.0) then 
            write(66,60)
            write(66,20)
            return 
         endif 
         if (sval .gt. 200.0) then 
            write(66,50)
            write(66,20)
            return 
         endif 
         if (sval .lt. 200.0) write(66,40)
         write(66,20)
c
       endif

      return
      end


      subroutine print_rhf_orbitals(array1,sval,a1,a2,b1,b2) 

c---------------------------------------------------------------------------
c
c 
c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
c
      integer i, ival, amax, a, a1, a2, b, b1, b2   
      double precision array1(a1:a2,b1:b2), sval, stemp   
      double precision efinal, conv   
      logical breturn 

      common /orbread/breturn 

       ival = sval 
       conv = 27.21138505 
c      if (me .eq. my_company_rank)then

       if (me .eq. 0)then

       if (ival .eq. 0) then  
          write(66,100) 
          open(111,file='ca.data') 
          write(111,*) a2-a1+1  
          do a = a1, a2   
          do b = b1, b2   
             write(111,*) array1(a,b) 
          enddo 
          enddo 
          close (111) 
       endif 

       if (ival .eq. -1) then  
          write(66,200) 
          open(222,file='cb.data') 
          write(222,*) a2-a1+1  
          do a = a1, a2   
          do b = b1, b2   
             write(222,*) array1(a,b) 
          enddo 
          enddo 
          close (222) 
       endif 

       endif 

       breturn = .false.

100    FORMAT(5x,'PRINTING RHF COEFFICIENTS TO FILE') 
200    FORMAT(5x,'PRINTING UHF COEFFICIENTS TO FILE') 

      return
      end
c.......................................................
c......................Prakash.............................................

      subroutine print_ccsdTx_cinfo(array1,sval,a1,a2,b1,b2)
      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
      integer i, ival, amax, a, a1, a2, b, b1, b2
      double precision array1(a1:a2,b1:b2), sval, stemp
      double precision efinal, conv
      open(66,file='summary.out')
       ival = sval
       conv = 27.21138505
c       if (me .eq. my_company_rank)then
       if (me .eq. 0)then
       if (ival .eq. 6) then
          do a = a1, a2
            if (a .eq. 1.0) then
               write(66,*)
               write(66,*)
               write(66,108)
               write(66,101)
               write(66,110)
               write(66,101)
            endif
c           if (dabs(array1(a,b1)) .gt. 0.0001) then 
               write(66,99) a-1, (array1(a,i), i=b1, b2)
               write(66,101)
               write(66,*)
               write(66,*)
               return
c           endif 
          enddo
       endif
       if (ival .eq. 700) then
          do a = a2, a1, -1
            if (a .eq. 1.0) then
               write(66,*)
               write(66,*)
               write(66,98)
               write(66,101)
               write(66,100)
               write(66,101)
            endif
            if (dabs(array1(a,b1)) .gt. 0.0001) then
               write(66,99) a, (array1(a,i), i=b1, b2)
               return
            endif
          enddo
       endif
       if (ival .eq. 701) then
          do a = a1, a2
            if (a .eq. 1.0) then
               write(66,*)
               write(66,*)
               write(66,97)
               write(66,101)
               write(66,100)
               write(66,101)
            endif
            if (dabs(array1(a,b1)) .lt. 0.0001) then
c              write(66,99) a-1, (array1(a-1,i), i=b1, b2) 
               amax = a - 1
               efinal = array1(a-1,2)
            go to 10
            endif
            write(66,99) a, (array1(a,i), i=b1, b2)
          enddo
10        continue
          write(66,101)
          write(66,102) efinal, amax
          write(66,101)
          write(66,*)
          write(66,*)
       endif
      endif
97    format(2x,"SUMMARY OF CCSD T X-derivative CALCULATION")
98    format(2x,"START OF CCSD T X-derivative CALCULATION")
108   format(2x,"SUMMARY OF MBPT(2) CALCULATION")
99    format(2x,I8,5x,F20.10,2x,F20.10)
100   format(2x,"Iteration#",9x,"CCSD Tx corr energy",7x,"Total energy")
110   format(2x,"Iteration#",9x,"MBPT(2) corr energy",7x,"Total energy")
101   format(2x,60('-'))
102   format(2x,"Final CCSD energy = ",F20.10,2x,"After",2x,I8,2x,
     *       "iterations")
      return
      end
c...............................................................
      subroutine print_pol_info(array1,sval,a1,a2,b1,b2)
      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
      integer i, ival, amax, a, a1, a2, b, b1, b2
      double precision array1(a1:a2,b1:b2), sval, stemp
      double precision efinal, conv
      open(66,file='summary.out')
       ival = sval
       conv = 27.21138505
c       if (me .eq. my_company_rank)then
       if (me .eq. 0)then
       if (ival .eq. 711) then
          do a = a2, a1, -1
            if (a .eq. 1.0) then
               write(66,*)
               write(66,*)
               write(66,101)
               write(66,*)
               write(66,98)
               write(66,*)
               write(66,101)
            endif
c            if (dabs(array1(a,b1)) .gt. 0.0001) then
c               write(66,99) a, (array1(a,i), i=b1, b2)
c               return
c            endif
          enddo
       endif
       if (ival .eq. 712) then
          write(66,97)
c          do a = a1, a2
               write(66,102)  (array1(1,i), i=b1, b2)
               write(66,103)  (array1(2,i), i=b1, b2)
               write(66,104)  (array1(3,i), i=b1, b2)
               write(66,*)
               Write(66,105)  (array1(1,1)+array1(2,2)+array1(3,3))/3
c          enddo 
       endif
      endif
97    format(17x," X                        Y                      Z")
98    format(26x,"CCSD Polarizability Results (a.u.)")
99    format(2x,I8,5x,F20.10,2x,F20.10)
101   format(2x,75('-'))
102   format(2x,"X",F20.7,4x,F20.7,4x,F20.7)
103   format(2x,"Y",F20.7,4x,F20.7,4x,F20.7)
104   format(2x,"Z",F20.7,4x,F20.7,4x,F20.7)
105   format(2x,"Isotropic =",F20.7)
      return
      end
c
c...............................................................
c
      subroutine print_gten_mc_info(array1,sval,a1,a2,b1,b2)
      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
      integer i, ival, amax, a, a1, a2, b, b1, b2
      double precision array1(a1:a2,b1:b2), sval, stemp, FNS, PPM
      double precision efinal, conv, Fac, Two, Gfac, Refspin,
     &                  Nunpaired 
      open(66,file='summary.out')
      ival = sval
      conv = 27.21138505

c Fine structure constant and parts per million

      FNS = 7.29735308D-03
      PPM = 1.0D+06
      GFAC= 2.0023193044D0
      Two = 2.0D0
c
c External multiplication constant for Relativistic Mass correction 
C  
      Nunpaired = Nalpha_occupied - Nbeta_occupied

      If (Nunpaired .EQ. 0) Then
         Write(6, "(a)") "Invalid input: A closed shell EPR calculation"
         Call abort_job()
      Endif

      Refspin   = Dble(Nunpaired)*0.50D0
      Fac = -1.0D0*FNS*FNS*PPM*(1.0D0/(Refspin))

       if (me .eq. 0)then
       if (ival .eq. 713) then
          do a = a2, a1, -1
            if (a .eq. 1.0) then
               write(66,*)
               write(66,*)
               write(66,101)
               write(66,*)
               write(66,98)
               write(66,*)
               write(66,101)
            endif
          enddo
       endif

       if (ival .eq. 713) then
          write(66,97)
c          do a = a1, a2
               write(66,102)  (Fac*array1(1,i), i=b1, b2)
               write(66,103)  (Fac*array1(2,i), i=b1, b2)
               write(66,104)  (Fac*array1(3,i), i=b1, b2)
               write(66,*)
c          enddo
       endif
      endif
97    format(17x," X                        Y                      Z")
98    format(26x,"CCSD g-tensor: Relativistic Mass Correction (ppm)")
99    format(2x,I8,5x,F20.10,2x,F20.10)
101   format(2x,75('-'))
102   format(2x,"X",F20.7,4x,F20.7,4x,F20.7)
103   format(2x,"Y",F20.7,4x,F20.7,4x,F20.7)
104   format(2x,"Z",F20.7,4x,F20.7,4x,F20.7)
      Return
      End
c
c-------------------------------------------------------------------------
c
      subroutine print_gten_1egc_info(array1,sval,a1,a2,b1,b2)
      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
      integer i, ival, amax, a, a1, a2, b, b1, b2
      double precision array1(a1:a2,b1:b2), sval, stemp, FNS, PPM
      double precision efinal, conv, Fac, Gfac,  Refspin,
     &                  Nunpaired
      open(66,file='summary.out')
      ival = sval
      conv = 27.21138505

c Fine structure constant and parts per million

      FNS = 7.29735308D-03
      GFAC= 2.0023193044D0
      PPM = 1.0D+06

c External multiplication constant for 1-el gauge correction 

      Nunpaired = Nalpha_occupied - Nbeta_occupied

      If (Nunpaired .EQ. 0) Then
         Write(6, "(a)") "Invalid input: A closed shell EPR calculation"
         Call abort_job()
      Endif 

      Refspin   = Dble(Nunpaired)*0.50D0
      Fac = FNS*FNS*PPM*(1.0D0/(4.0D0*Refspin))

       if (me .eq. 0)then
       if (ival .eq. 714) then
          do a = a2, a1, -1
            if (a .eq. 1.0) then
               write(66,*)
               write(66,*)
               write(66,101)
               write(66,*)
               write(66,98)
               write(66,*)
               write(66,101)
            endif
          enddo
       endif

       if (ival .eq. 714) then
          write(66,97)
c          do a = a1, a2
               write(66,102)  (Fac*array1(1,i), i=b1, b2)
               write(66,103)  (Fac*array1(2,i), i=b1, b2)
               write(66,104)  (Fac*array1(3,i), i=b1, b2)
               write(66,*)
c          enddo
       endif
      endif
97    format(17x," X                        Y                      Z")
98    format(26x,"CCSD g-tensor: 1el-Diamagnetic Spin-orbit ",
     &       58x, "Gauge Correction (ppm)")
99    format(2x,I8,5x,F20.10,2x,F20.10)
101   format(2x,75('-'))
102   format(2x,"X",F20.7,4x,F20.7,4x,F20.7)
103   format(2x,"Y",F20.7,4x,F20.7,4x,F20.7)
104   format(2x,"Z",F20.7,4x,F20.7,4x,F20.7)
      Return
      End
c
c---------------------------------------------------------------------
c
      subroutine print_gten_1bl_info(array1,sval,a1,a2,b1,b2)
      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
      integer i, ival, amax, a, a1, a2, b, b1, b2
      integer Nunpaired 
      double precision array1(a1:a2,b1:b2), sval, stemp, FNS, PPM
      double precision efinal, conv, Fac, Gfac, Refspin
      open(66,file='summary.out')
      ival = sval
      conv = 27.21138505

c Fine structure constant and parts per million

       FNS = 7.29735308D-03
       GFAC= 2.0023193044D0
       PPM = 1.0D+06
c
c External multiplication constant for paramagnetic spin-orbit 
c contribution.

      Nunpaired = Nalpha_occupied - Nbeta_occupied

      If (Nunpaired .EQ. 0) Then
         Write(6, "(a)") "Invalid input: A closed shell EPR calculation"
         Call abort_job()
      Endif

      Refspin   = Dble(Nunpaired)*0.50D0
  
      Fac = FNS*FNS*PPM*(1.0D0/(4.0D0*Refspin))

       if (me .eq. 0)then
       if (ival .eq. 715) then
          do a = a2, a1, -1
            if (a .eq. 1.0) then
               write(66,*)
               write(66,*)
               write(66,101)
               write(66,*)
               write(66,98)
               write(66,*)
               write(66,101)
            endif
          enddo
       endif

       if (ival .eq. 715) then
          write(66,97)
c          do a = a1, a2
               write(66,102)  (Fac*array1(1,i), i=b1, b2)
               write(66,103)  (Fac*array1(2,i), i=b1, b2)
               write(66,104)  (Fac*array1(3,i), i=b1, b2)
               write(66,*)
c          enddo
       endif
      endif
97    format(17x," X                        Y                      Z")
98    format(26x,"CCSD g-tensor: 1el-Paramagnetc Spin-orbit (ppm)")
99    format(2x,I8,5x,F20.10,2x,F20.10)
101   format(2x,75('-'))
102   format(2x,"X",F20.7,4x,F20.7,4x,F20.7)
103   format(2x,"Y",F20.7,4x,F20.7,4x,F20.7)
104   format(2x,"Z",F20.7,4x,F20.7,4x,F20.7)
      Return
      End
c
c---------------------------------------------------------------------
c
      subroutine print_gten_2bl_info(array1,sval,a1,a2,b1,b2)
      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
      integer i, ival, amax, a, a1, a2, b, b1, b2
      integer Nunpaired
      double precision array1(a1:a2,b1:b2), sval, stemp, FNS, PPM
      double precision efinal, conv, Fac, Gfac, Refspin
      open(66,file='summary.out')
      ival = sval
      conv = 27.21138505

c Fine structure constant and parts per million

       FNS = 7.29735308D-03
       GFAC= 2.0023193044D0
       PPM = 1.0D+06
c
c External multiplication constant for paramagnetic spin-orbit
c contribution.

      Nunpaired = Nalpha_occupied - Nbeta_occupied

      If (Nunpaired .EQ. 0) Then
         Write(6, "(a)") "Invalid input: A closed shell EPR calculation"
         Call abort_job()
      Endif

      Refspin   = Dble(Nunpaired)*0.50D0

      Fac = -1.0D0*FNS*FNS*PPM*(1.0D0/(2.0D0*Refspin))

       if (me .eq. 0)then
       if (ival .eq. 716) then
          do a = a2, a1, -1
            if (a .eq. 1.0) then
               write(66,*)
               write(66,*)
               write(66,101)
               write(66,*)
               write(66,98)
               write(66,*)
               write(66,101)
            endif
          enddo
       endif

       if (ival .eq. 716) then
          write(66,97)
c          do a = a1, a2
               write(66,102)  (Fac*array1(1,i), i=b1, b2)
               write(66,103)  (Fac*array1(2,i), i=b1, b2)
               write(66,104)  (Fac*array1(3,i), i=b1, b2)
               write(66,*)
c          enddo
       endif
      endif
97    format(17x," X                        Y                      Z")
98    format(26x,"CCSD g-tensor: 2el-Paramagnetc Spin-orbit (ppm)")
99    format(2x,I8,5x,F20.10,2x,F20.10)
101   format(2x,75('-'))
102   format(2x,"X",F20.7,4x,F20.7,4x,F20.7)
103   format(2x,"Y",F20.7,4x,F20.7,4x,F20.7)
104   format(2x,"Z",F20.7,4x,F20.7,4x,F20.7)
      Return
      End
c
c---------------------------------------------------------------------
c
      subroutine print_gten_2egc_info(array1,sval,a1,a2,b1,b2)
      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
      integer i, ival, amax, a, a1, a2, b, b1, b2
      integer Nunpaired
      double precision array1(a1:a2,b1:b2), sval, stemp, FNS, PPM
      double precision efinal, conv, Fac, Gfac, Refspin
      open(66,file='summary.out')
      ival = sval
      conv = 27.21138505

c Fine structure constant and parts per million

       FNS = 7.29735308D-03
       GFAC= 2.0023193044D0
       PPM = 1.0D+06
c
c External multiplication constant for diamagentic spin-orbit
c contribution.

      Nunpaired = Nalpha_occupied - Nbeta_occupied

      If (Nunpaired .EQ. 0) Then
         Write(6, "(a)") "Invalid input: A closed shell EPR calculation"
         Call abort_job()
      Endif

      Refspin   = Dble(Nunpaired)*0.50D0

      Fac = -1.0D0*FNS*FNS*PPM*(1.0D0/(2.0D0*Refspin))

       if (me .eq. 0)then
       if (ival .eq. 717) then
          do a = a2, a1, -1
            if (a .eq. 1.0) then
               write(66,*)
               write(66,*)
               write(66,101)
               write(66,*)
               write(66,98)
               write(66,*)
               write(66,101)
            endif
          enddo
       endif

       if (ival .eq. 717) then
          write(66,97)
c          do a = a1, a2
               write(66,102)  (Fac*array1(1,i), i=b1, b2)
               write(66,103)  (Fac*array1(2,i), i=b1, b2)
               write(66,104)  (Fac*array1(3,i), i=b1, b2)
               write(66,*)
c          enddo
       endif
      endif
97    format(17x," X                        Y                      Z")
98    format(26x,"CCSD g-tensor: 2el-diamagnetic shielding ", 
     &       58x, "Gauge Correction (ppm)")

99    format(2x,I8,5x,F20.10,2x,F20.10)
101   format(2x,75('-'))
102   format(2x,"X",F20.7,4x,F20.7,4x,F20.7)
103   format(2x,"Y",F20.7,4x,F20.7,4x,F20.7)
104   format(2x,"Z",F20.7,4x,F20.7,4x,F20.7)
      Return
      End
c
c...............................................................
c`
      subroutine print_dipole_info(array1,sval,a1,a2,b1,b2)
      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
      integer i, ival, amax, a, a1, a2, b, b1, b2
      double precision array1(a1:a2,b1:b2),sval
      open(66,file='summary.out')
      ival = sval

       if (me .eq. 0)then

C Write the DIPOLE file; Ajith Perera, 12/2014.

         if (ival .eq. 718) then
            open (unit=80,file='DIPOL',status='unknown',
     &            form="formatted")
            write(6,"(a)") ' Archiving Dipole moments'
            rewind(80)
            write(80,"(3F20.10)") array1(1,1), array1(2,2),
     &                            array1(3,3)
            close(80)

         Write(66,*) 
         Write(66,*) 
         write(66,101)
         write(66,98)
         write(66,101)
         write(66,97)
         write(66,102) array1(1,1), array1(2,2), array1(3,3)
         Write(66,*) 
         Write(66,*) 
         write(66,101)
         endif

       endif
97    format(17x," X                        Y                      Z")
98    format(17x,"Correlated (CCSD or MBPT(2) Dipole Moments (a.u.)")
101   format(2x,75('-'))
102   format(2x,F20.7,4x,F20.7,4x,F20.7)

       Return
       end
C----------------------------------------------------------

      subroutine print_lccd_cinfo(array1,sval,a1,a2,b1,b2) 

c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'parallel_info.h'
c
      integer i, ival, amax, a, a1, a2, b, b1, b2   
      double precision array1(a1:a2,b1:b2), sval, stemp   
      double precision efinal, conv   

      open(66,file='summary.out') 

       ival = sval 
       conv = 27.21138505 
c      if (me .eq. my_company_rank)then

       if (me .eq. 0)then

       if (ival .eq. 32) then  
          do a = a2, a1, -1  
            if (a .eq. 1.0) then 
               write(66,*) 
               write(66,*) 
               write(66,98) 
               write(66,101) 
               write(66,100) 
               write(66,101) 
            endif  
            if (dabs(array1(a,b1)) .gt. 0.0001) then 
               write(66,99) a, (array1(a,i), i=b1, b2) 
               return 
            endif 
          enddo 
       endif 

       if (ival .eq. 33) then  
          do a = a1, a2 
            if (a .eq. 1.0) then 
               write(66,*) 
               write(66,*) 
               write(66,97) 
               write(66,101) 
               write(66,100) 
               write(66,101) 
            endif  
            if (dabs(array1(a,b1)) .lt. 0.0001) then 
c              write(66,99) a-1, (array1(a-1,i), i=b1, b2) 
               amax = a - 1  
               efinal = array1(a-1,2) 
            go to 10 
            endif 
            write(66,99) a, (array1(a,i), i=b1, b2) 
          enddo 
10        continue 
          write(66,101) 
          write(66,102) efinal, amax 
          write(66,101) 
          write(66,*) 
          write(66,*) 

       endif 

      endif

97    format(2x,"SUMMARY OF LCCD CALCULATION") 
98    format(2x,"START OF LCCD CALCULATION") 
99    format(2x,I8,5x,F20.10,2x,F20.10) 
100   format(2x,"Iteration#",9x,"LCCD corr energy",7x,"Total energy")
101   format(2x,60('-'))
102   format(2x,"Final LCCD energy = ",F20.10,2x,"After",2x,I8,2x,
     *       "iterations")

      return
      end
