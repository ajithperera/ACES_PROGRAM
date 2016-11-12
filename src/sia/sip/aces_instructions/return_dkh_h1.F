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
      subroutine return_dkh_h1(array_table, 
     *                         narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c   Computes the Douglas-Kroll-Hess transformed H1 array.
c
c   execute array1 DKH_order ; DKH_order is order of transformation
c 
c--------------------------------------------------------------------------

      implicit none

      include 'mpif.h'
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'trace.h'
      include 'parallel_info.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer narray_table, nindex_table, nsegment_table,
     *        nblock_map_table, array
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table), get_index_from_base
      integer*8 get_block_data_index
      integer find_current_block, block

      integer array_type, nindex, inda(2), indb(2)

      logical direct_flag

      integer i, j, k
      integer nind, nsend
      integer ierr, handle
      integer blk, blkndx, maxblk
      integer*8 iarray, pvp_array, indblk, get_block_index
      integer*8 ibarray, tmp_array
      integer*8 nai_array, kin_array, ovl_array
      integer*8 integral_scr
      integer ind, seg, val1(mx_array_index), val2(mx_array_index)
      integer msg(len_sip_server_message)
      integer flag
      integer stack
      integer sind

      integer*8  dkh_array
      integer  evec_array, evec_type, nindex_evec, order

      double precision x(1)
      double precision dbuf(1), temp
      integer ibuf(1)
#ifdef ALTIX
      pointer (dptr, dbuf)
      pointer (xptr, x)
      pointer (iptr, ibuf)
#else
      common dbuf, x
      equivalence (ibuf(1), dbuf(1))
#endif

      common /d2int_com/jatom, jx, jcenter
      integer jatom, jx, jcenter

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
      xptr = dshptr
#endif

      direct_flag = .false. 

c----------------------------------------------------------------------------
c   Set up the integral pointer 
c----------------------------------------------------------------------------

      integral_scr = dpkgscr
      maxblk  = zmax

#ifdef ALTIX
      dptr = dshptr
#endif

c----------------------------------------------------------------------------
c   Locate the scalar that defines the order of transformation
c---------------------------------------------------------------------------
      array      = op(c_result_array)
      array_type = array_table(c_array_type,array)
      nindex     = array_table(c_nindex,array)

      if (array_type .ne. scalar_value) return

      if (array .lt. 1 .or. array .gt. narray_table) then
         print *,'Error: Invalid scalar in return_dkh , line ',
     *     current_line
         print *,'Array index is ',array,' Allowable values are ',
     *      ' 1 through ',narray_table
         call abort_job()
      endif 

      sind =  array_table(c_scalar_index, array)
      if (sind .lt. 1 .or. sind .gt. nscalar_table) then
         print *,'Scalar table index out of range in return_dkh , ',
     *           'line ',current_line
         print *,'Index for array ',array,' is ',ind,' should be ',
     *           'between 1 and ',nscalar_table
         call abort_job()
      endif
C
C   Define the order of the calculation
C
      order = int(scalar_table(sind))

c----------------------------------------------------------------------------
c   Get the Nuclear attraction integrals.
c----------------------------------------------------------------------------

      array      = 11 
      array_type = array_table(c_array_type, array)
      nindex     = array_table(c_nindex, array)

      if (array_type .ne. static_array) then
         print *,'Error: Nuclear attraction integrals must be static.'
         call abort_job()
      endif

      if (nindex .ne. 2) then
         print *,'Error: NAI array must have two indices.'
         call abort_job()
      endif

      do i = 1, nindex
         na = array_table(c_index_array1+i-1,array)
         inda(i) = index_table(c_index_size, na)  ! pick up length of index
      enddo

      nai_array = get_index_from_base(address_table(array), x, 2)

c----------------------------------------------------------------------------
c   Get the Kinetic energy integrals.
c----------------------------------------------------------------------------

      array      = 12 
      array_type = array_table(c_array_type, array)
      nindex     = array_table(c_nindex, array)


      if (array_type .ne. static_array) then
         print *,'Error: Kinetic energy integrals must be static.'
         call abort_job()
      endif

      if (nindex .ne. 2) then
         print *,'Error: KIN array must have two indices.'
         call abort_job()
      endif

      do i = 1, nindex
         na = array_table(c_index_array1+i-1,array)
         inda(i) = index_table(c_index_size, na)  ! pick up length of index
      enddo

      kin_array = get_index_from_base(address_table(array), x, 2)

c----------------------------------------------------------------------------
c   Get the Overlap integrals.
c----------------------------------------------------------------------------

      array      = 13 
      array_type = array_table(c_array_type, array)
      nindex     = array_table(c_nindex, array)

      if (array_type .ne. static_array) then
         print *,'Error: overlap integrals must be static.'
         call abort_job()
      endif

      if (nindex .ne. 2) then
         print *,'Error: ovl array must have two indices.'
         call abort_job()
      endif

      do i = 1, nindex
         na = array_table(c_index_array1+i-1,array)
         inda(i) = index_table(c_index_size, na)  ! pick up length of index
      enddo

      ovl_array = get_index_from_base(address_table(array), x, 2)

c----------------------------------------------------------------------------
c   Get the alpha Fock array.  This will be PVP for now!
c----------------------------------------------------------------------------

      array = 17 ! Fockrohf_a for a placeholder   
      array_type = array_table(c_array_type, array)
      nindex      = array_table(c_nindex, array)

      if (array_type .ne. static_array) then
         print *,'Error: Fock arrays in integral must be static.'
         call abort_job()
      endif

      if (nindex .ne. 2) then
         print *,'Error: Fock array must have two indices.'
         call abort_job()
      endif

      do i = 1, nindex
         na = array_table(c_index_array1+i-1,array)
         inda(i) = index_table(c_index_size, na)  ! pick up length of index
      enddo

      pvp_array = get_index_from_base(address_table(array), x, 2)

c----------------------------------------------------------------------------
c   Get the beta Fock array. This is used as temp OVL array, since 
c   those integrals change.
c----------------------------------------------------------------------------

      array = 18 ! Fockrohf_b for a placeholder   
      array_type = array_table(c_array_type, array)
      nindex      = array_table(c_nindex, array)

      if (array_type .ne. static_array) then
         print *,'Error: Fock arrays in integral must be static.'
         call abort_job()
      endif

      if (nindex .ne. 2) then
         print *,'Error: Fock array must have two indices.'
         call abort_job()
      endif

      do i = 1, nindex
         na = array_table(c_index_array1+i-1,array)
         indb(i) = index_table(c_index_size, na)  ! pick up length of index
      enddo

      tmp_array = get_index_from_base(address_table(array), x, 2)

C---------------------------------------------------------------------------
C
C
C     Do the Douglas Kroll Hess transformation!
C
C
      call do_dkh_trans ( x(nai_array),x(kin_array),
     *                    x(pvp_array),x(ovl_array),
     *                    x(tmp_array),
     *                    1,inda(1),1,inda(2),order,
     *                    zmax,dbuf(dpkgscr))

      return
      end
C
C
C-------------------------------------------------------------------------
C
C
      subroutine do_dkh_trans(nai,kin,pvp,ovl,ovl_tmp,
     *                        nc1,nc2,nd1,nd2,order,zmax,zcore)
      implicit none
      integer nc1,nc2,nd1,nd2,nbas, order
      integer zmax,nbas1,nbas2
      double precision  nai(nc1:nc2,nd1:nd2)
      double precision  kin(nc1:nc2,nd1:nd2)
      double precision  pvp(nc1:nc2,nd1:nd2)
      double precision  ovl(nc1:nc2,nd1:nd2)
      double precision  ovl_tmp (nc1:nc2,nd1:nd2)
      double precision  zcore(1:ZMAX)

      integer I, J, iprint, ierr
C
C
C     Check to see that the matrices given are square!
C
C
      NBAS1 = NC1 + NC2 - 1
      NBAS2 = ND1 + ND2 - 1
      IF (NBAS1 .NE. NBAS2) THEN
         WRITE (*,*) "Matrices are not square in return_dkh_h1."
         CALL  abort_job()
      ENDIF
      NBAS = NBAS1
C
C
C     Zero out PVP since we do not have those integrals yet.
C     Also, put OVL in OVL_TMP, since the OVL integrals are modified.
C
C
      DO I = NC1, NC2
      DO J = ND1, ND2
         OVL_TMP (I,J) = OVL (I,J)
      ENDDO
      ENDDO
C
C
C     Perform the DKH transformation
C
C
C      WRITE (*,*) ' PVP integrals have not been debugged yet! ',
C     +            ' Do not use DKH yet! Aborting job! '
C      call abort_job ()
C
C
C
      iprint = 99
      Write(6,*) "At entry ZMAX:" , ZMAX
      CALL  DKH__MAIN_TRANS (ORDER,NBAS,iprint,ZMAX,
     +                       zcore,OVL,PVP,ierr,KIN,NAI )

      if (ierr .eq. 1) call abort_job()
C
C
C     Put the OVL integrals back into proper form!
C
C
      DO I = NC1, NC2
      DO J = ND1, ND2
         OVL (I,J) = OVL_TMP (I,J)
      ENDDO
      ENDDO

      return
      end
