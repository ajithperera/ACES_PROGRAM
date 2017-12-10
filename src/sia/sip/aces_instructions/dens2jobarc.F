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
      subroutine dens2jobarc(array_table,
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)
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

      integer array_type, nindex, inda(2), indb(2), n1, n2

      logical direct_flag

      integer i, j, k
      integer nind, nsend
      integer ierr, handle
      integer blk, blkndx, maxblk
      integer*8 iarray, dens_array, indblk, get_block_index
      integer*8 ibarray, fbarray 
      integer ind, seg, val1(mx_array_index), val2(mx_array_index)
      integer msg(len_sip_server_message)
      integer flag
      integer stack
      integer sind,type  

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

c----------------------------------------------------------------------------
c   Get the density array
c----------------------------------------------------------------------------

      array = op(c_result_array)
      array_type = array_table(c_array_type, array)

      if (array_type .ne. static_array) then
         print *,'Error: density array must be static in ',
     *           'dens2jobarc'
         call abort_job()
      endif

      nindex      = array_table(c_nindex, array)

      if (nindex .ne. 2) then
         print *,'Error: DENS array must have ',
     *           'two indices.'
         call abort_job()
      endif

      do i = 1, nindex
         na = array_table(c_index_array1+i-1,array)
         inda(i) = index_table(c_index_size, na)  ! pick up length of index
      enddo

      dens_array = get_index_from_base(address_table(array), x, 2)

c----------------------------------------------------------------------------
c   Get the type of density
c----------------------------------------------------------------------------
      array       = op(c_op1_array)
      array_type  = array_table(c_array_type, array)
      nindex      = array_table(c_nindex, array)

      if (array_type .ne. scalar_value) then
         print *,'Error: scalar not sent in dens2jobarc.'
         print *,(op(i),i=1,loptable_entry)
         call abort_job()
      endif

      sind =  array_table(c_scalar_index, array)
      if (sind .lt. 1 .or. sind .gt. nscalar_table) then
         print *,'Scalar table index out of range in dens2jobarc, ',
     *           'line ',current_line
         print *,'Index for array ',array,' is ',sind,' should be ',
     *           'between 1 and ',nscalar_table
         call abort_job()
      endif
      type = INT (scalar_table(sind))

      CALL  DO_THE_DENS_WRITING (x(dens_array),inda(1),inda(2),type,
     +                      ibuf,dbuf)

      return
      end
c---------------------------------------------------------------------------

      SUBROUTINE  DO_THE_DENS_WRITING (DENS,NDIM1,NDIM2,TYPE,ibuf,dbuf)

      IMPLICIT    NONE

      INCLUDE 'interpreter.h'
      INCLUDE 'int_gen_parms.h'
      INCLUDE 'parallel_info.h'

      INTEGER    NDIM,NDIM1,NDIM2,MU,NU,IEOM_ROOTS 
      INTEGER    ierd,ierd1,ierd2
      INTEGER    TYPE,FLAG,LENGTH
      INTEGER    ibuf (*)
      INTEGER    NROOTS 

      CHARACTER*8  JALABEL

      character*8 denstyper_aa(24)
      character*8 denstypel_aa(24)
      character*8 denstyper_bb(24)
      character*8 denstypel_bb(24)
      character*8 denstyper(24)
      character*8 denstypel(24)
      character*8 denstype_tl (24)

      data denstyper_aa /'REOMDN1A','REOMDN2A','REOMDN3A','REOMDN4A',
     +                   'REOMDN5A','REOMDN6A','REOMDN7A','REOMDN8A',
     +                   'REOMDN9A','REOMD10A','REOMD11A','REOMD12A',
     +                   'REOMD13A','REOMD14A','REOMD15A','REOMD16A',
     +                   'REOMD17A','REOMD18A','REOMD19A','REOMD20A',
     +                   'REOMD21A','REOMD22A','REOMD23A','REOMD24A'/
      data denstyper_bb /'REOMDN1B','REOMDN2B','REOMDN3B','REOMDN4B',
     +                   'REOMDN5B','REOMDN6B','REOMDN7B','REOMDN8B',
     +                   'REOMDN9B','REOMD10B','REOMD11B','REOMD12B',
     +                   'REOMD13B','REOMD14B','REOMD15B','REOMD16B',
     +                   'REOMD17B','REOMD18B','REOMD19B','REOMD20B',
     +                   'REOMD21B','REOMD22B','REOMD23B','REOMD24B'/
      data denstypel_aa /'LEOMDN1A','LEOMDN2A','LEOMDN3A','LEOMDN4A',
     +                   'LEOMDN5A','LEOMDN6A','LEOMDN7A','LEOMDN8A',
     +                   'LEOMDN9A','LEOMD10A','LEOMD11A','LEOMD12A',
     +                   'LEOMD13A','LEOMD14A','LEOMD15A','LEOMD16A',
     +                   'LEOMD17A','LEOMD18A','LEOMD19A','LEOMD20A',
     +                   'LEOMD21A','LEOMD22A','LEOMD23A','LEOMD24A'/
      data denstypel_bb /'LEOMDN1B','LEOMDN2B','LEOMDN3B','LEOMDN4B',
     +                   'LEOMDN5B','LEOMDN6B','LEOMDN7B','LEOMDN8B',
     +                   'LEOMDN9B','LEOMD10B','LEOMD11B','LEOMD12B',
     +                   'LEOMD13B','LEOMD14B','LEOMD15B','LEOMD16B',
     +                   'LEOMD17B','LEOMD18B','LEOMD19B','LEOMD20B',
     +                   'LEOMD21B','LEOMD22B','LEOMD23B','LEOMD24B'/
      data denstype_tl / 'EXCDEN01','EXCDEN02','EXCDEN03','EXCDEN04',
     +                   'EXCDEN05','EXCDEN06','EXCDEN07','EXCDEN08',
     +                   'EXCDEN09','EXCDEN10','EXCDEN11','EXCDEN12',
     +                   'EXCDEN13','EXCDEN14','EXCDEN15','EXCDEN16',
     +                   'EXCDEN17','EXCDEN18','EXCDEN19','EXCDEN20',
     +                   'EXCDEN21','EXCDEN22','EXCDEN23','EXCDEN24'/

      data denstyper   / 'REOMDN01','REOMDN02','REOMDN03','REOMDN04',
     +                   'REOMDN05','REOMDN06','REOMDN07','REOMDN08',
     +                   'REOMDN09','REOMDN10','REOMDN11','REOMDN12',
     +                   'REOMDN13','REOMDN14','REOMDN15','REOMDN16',
     +                   'REOMDN17','REOMDN18','REOMDN19','REOMDN20',
     +                   'REOMDN21','REOMDN22','REOMDN23','REOMDN24'/ 

      data denstypel   / 'LEOMDN01','LEOMDN02','LEOMDN03','LEOMDN04',
     +                   'LEOMDN05','LEOMDN06','LEOMDN07','LEOMDN08',
     +                   'LEOMDN09','LEOMDN10','LEOMDN11','LEOMDN12',
     +                   'LEOMDN13','LEOMDN14','LEOMDN15','LEOMDN16',
     +                   'LEOMDN17','LEOMDN18','LEOMDN19','LEOMDN20',
     +                   'LEOMDN21','LEOMDN22','LEOMDN23','LEOMDN24'/ 

      DOUBLE PRECISION  dscale

      DOUBLE PRECISION  DENS (1:NDIM1,1:NDIM2)
      DOUBLE PRECISION  SCR  (1:NDIM1,1:NDIM2)
      DOUBLE PRECISION  EOM_EXE_S(12), EOM_EXE_T(12)
      DOUBLE PRECISION  dbuf (*)

#ifdef _DEBUG_LVL0
      if (me .EQ. 0) then
CSSS         call output(DENS, 1, NDIM1, 1, NDIM2, NDIM1, NDIM2, 1)
         write(6,*) "TYPE, iuhf :", TYPE, iuhf 
      Endif 
    
#endif 
C Block for the UHF wave-functios

      If (Iuhf .Eq. 1) Then 

      FLAG = -1 
      IF (me .EQ. 0) THEN
 
         IF (type .eq. -1) CALL DPUTREC (1, 'JOBARC',"CIS     ",
     &                                   1,SCR)
         IF (type .eq. -2) THEN
             NROOTS = EOM_ROOTS 
             CALL DGETREC(1, 'JOBAC', 'EOM_EXES', NROOTS, EOM_EXE_S)
             CALL DGETREC(0, 'JOBAC', 'EOM_EXET', LENGTH, EOM_EXE_T)

C Check whether EOM_EXET exist before try to read it. This can happen
C if the user specified UHF reference and ran the RHF sial codes.There
C is no reason to crash!

             IF (LENGTH .GT. 0) 
     &       CALL DGETREC(1, 'JOBAC', 'EOM_EXET', NROOTS, EOM_EXE_T)
             CALL ACES_JA_TRUNCATE('MARKER  ',1)
             CALL DPUTREC (1, 'JOBARC',"EOM     ", 1, SCR)
             CALL DPUTREC(1, 'JOBAC', 'EOM_EXES', NROOTS, EOM_EXE_S)
             IF (LENGTH .GT. 0) 
     &       CALL DPUTREC(1, 'JOBAC', 'EOM_EXET', NROOTS, EOM_EXE_T)
             RETURN
         ENDIF 

         IF      (type .EQ.  0) THEN
            JALABEL='HFDENSTY'
         ELSE IF (type .EQ.  1) THEN
            JALABEL='TDENSITY'
         ELSE IF (Type .EQ. 2 .OR. Type .EQ. 3 .OR. Type .EQ. 4 .OR.
     +            Type .EQ. 5 .OR. Type .EQ. 6) THEN
            If (Type .EQ. 2) JALABEL  = denstyper_aa(1)
            If (Type .EQ. 3) JALABEL  = denstyper_bb(1)
            If (Type .EQ. 4) JALABEL  = denstypel_aa(1)
            If (Type .EQ. 5) JALABEL  = denstypel_bb(1)
            If (Type .EQ. 6) JALABEL  = denstype_tl (1)
         ELSE IF (Type .EQ. 7 .OR. Type .EQ. 8 .OR. Type .EQ. 9 .OR.
     +            Type .EQ. 10 .OR. Type .EQ. 11) THEN
            If (Type .EQ. 7)  JALABEL  = denstyper_aa(2)
            If (Type .EQ. 8)  JALABEL  = denstyper_bb(2)
            If (Type .EQ. 9)  JALABEL  = denstypel_aa(2)
            If (Type .EQ. 10) JALABEL  = denstypel_bb(2)
            If (Type .EQ. 11) JALABEL  = denstype_tl (2)
         ELSE IF (Type .EQ. 12 .OR. Type .EQ. 13 .OR. Type .EQ. 14 .OR.
     +            Type .EQ. 15 .OR. Type .EQ. 16) THEN
            If (Type .EQ. 12) JALABEL  = denstyper_aa(3)
            If (Type .EQ. 13) JALABEL  = denstyper_bb(3)
            If (Type .EQ. 14) JALABEL  = denstypel_aa(3)
            If (Type .EQ. 15) JALABEL  = denstypel_bb(3)
            If (Type .EQ. 16) JALABEL  = denstype_tl (3)
         ELSE IF (Type .EQ. 17 .OR. Type .EQ. 18 .OR. Type .EQ. 19 .OR.
     +            Type .EQ. 20 .OR. Type .EQ. 21) THEN
            If (Type .EQ. 17) JALABEL  = denstyper_aa(4)
            If (Type .EQ. 18) JALABEL  = denstyper_bb(4)
            If (Type .EQ. 19) JALABEL  = denstypel_aa(4)
            If (Type .EQ. 20) JALABEL  = denstypel_bb(4)
            If (Type .EQ. 21) JALABEL  = denstype_tl (4)
         ELSE IF (Type .EQ. 22 .OR. Type .EQ. 23 .OR. Type .EQ. 24 .OR.
     +            Type .EQ. 25 .OR. Type .EQ. 26) THEN
            If (Type .EQ. 22) JALABEL  = denstyper_aa(5)
            If (Type .EQ. 23) JALABEL  = denstyper_bb(5)
            If (Type .EQ. 24) JALABEL  = denstypel_aa(5)
            If (Type .EQ. 25) JALABEL  = denstypel_bb(5)
            If (Type .EQ. 26) JALABEL  = denstype_tl (5)
         ELSE IF (Type .EQ. 27 .OR. Type .EQ. 28 .OR. Type .EQ. 29 .OR.
     +            Type .EQ. 30 .OR. Type .EQ. 31) THEN
            If (Type .EQ. 27) JALABEL  = denstyper_aa(6)
            If (Type .EQ. 28) JALABEL  = denstyper_bb(6)
            If (Type .EQ. 29) JALABEL  = denstypel_aa(6)
            If (Type .EQ. 30) JALABEL  = denstypel_bb(6)
            If (Type .EQ. 31) JALABEL  = denstype_tl (6)
         ELSE IF (Type .EQ. 32 .OR. Type .EQ. 33 .OR. Type .EQ. 34 .OR.
     +            Type .EQ. 35 .OR. Type .EQ. 36) THEN
            If (Type .EQ. 32) JALABEL  = denstyper_aa(7)
            If (Type .EQ. 33) JALABEL  = denstyper_bb(7)
            If (Type .EQ. 34) JALABEL  = denstypel_aa(7)
            If (Type .EQ. 35) JALABEL  = denstypel_bb(7)
            If (Type .EQ. 36) JALABEL  = denstype_tl (7)

         ELSE IF (Type .EQ. 37 .OR. Type .EQ. 38 .OR. Type .EQ. 39 .OR.
     +            Type .EQ. 40 .OR. Type .EQ. 41) THEN
            If (Type .EQ. 37) JALABEL  = denstyper_aa(8)
            If (Type .EQ. 38) JALABEL  = denstyper_bb(8)
            If (Type .EQ. 39) JALABEL  = denstypel_aa(8)
            If (Type .EQ. 40) JALABEL  = denstypel_bb(8)
            If (Type .EQ. 41) JALABEL  = denstype_tl (8)

         ELSE IF (Type .EQ. 42 .OR. Type .EQ. 43 .OR. Type .EQ. 44 .OR. 
     +            Type .EQ. 45 .OR. Type .EQ. 46) THEN
            If (Type .EQ. 42) JALABEL  = denstyper_aa(9)
            If (Type .EQ. 43) JALABEL  = denstyper_bb(9)
            If (Type .EQ. 44) JALABEL  = denstypel_aa(9)
            If (Type .EQ. 45) JALABEL  = denstypel_bb(9)
            If (Type .EQ. 46) JALABEL  = denstype_tl (9)

         ELSE IF (Type .EQ. 47 .OR. Type .EQ. 48 .OR. Type .EQ. 49 .OR.
     +            Type .EQ. 50 .OR. Type .EQ. 51) THEN
            If (Type .EQ. 47) JALABEL  = denstyper_aa(10)
            If (Type .EQ. 48) JALABEL  = denstyper_bb(10)
            If (Type .EQ. 49) JALABEL  = denstypel_aa(10)
            If (Type .EQ. 50) JALABEL  = denstypel_bb(10)
            If (Type .EQ. 51) JALABEL  = denstype_tl (10)

         ELSE IF (Type .EQ. 52 .OR. Type .EQ. 53 .OR. Type .EQ. 54 .OR.
     +            Type .EQ. 55 .OR. Type .EQ. 56) THEN
            If (Type .EQ. 52) JALABEL  = denstyper_aa(11)
            If (Type .EQ. 53) JALABEL  = denstyper_bb(11)
            If (Type .EQ. 54) JALABEL  = denstypel_aa(11)
            If (Type .EQ. 55) JALABEL  = denstypel_bb(11)
            If (Type .EQ. 56) JALABEL  = denstype_tl (11)

         ELSE IF (Type .EQ. 57 .OR. Type .EQ. 58 .OR. Type .EQ. 59 .OR.
     +            Type .EQ. 60 .OR. Type .EQ. 61) THEN
            If (Type .EQ. 57) JALABEL  = denstyper_aa(12)
            If (Type .EQ. 58) JALABEL  = denstyper_bb(12)
            If (Type .EQ. 59) JALABEL  = denstypel_aa(12)
            If (Type .EQ. 60) JALABEL  = denstypel_bb(12)
            If (Type .EQ. 61) JALABEL  = denstype_tl (12)

         ELSE IF (Type .EQ. 62 .OR. Type .EQ. 63 .OR. Type .EQ. 64 .OR.
     +            Type .EQ. 65 .OR. Type .EQ. 66) THEN
            If (Type .EQ. 62) JALABEL  = denstyper_aa(13)
            If (Type .EQ. 63) JALABEL  = denstyper_bb(13)
            If (Type .EQ. 64) JALABEL  = denstypel_aa(13)
            If (Type .EQ. 65) JALABEL  = denstypel_bb(13)
            If (Type .EQ. 66) JALABEL  = denstype_tl (13)

         ELSE IF (Type .EQ. 67 .OR. Type .EQ. 68 .OR. Type .EQ. 69 .OR.
     +            Type .EQ. 70 .OR. Type .EQ. 71) THEN
            If (Type .EQ. 67) JALABEL  = denstyper_aa(14)
            If (Type .EQ. 68) JALABEL  = denstyper_bb(14)
            If (Type .EQ. 69) JALABEL  = denstypel_aa(14)
            If (Type .EQ. 70) JALABEL  = denstypel_bb(14)
            If (Type .EQ. 71) JALABEL  = denstype_tl (14)
         ELSE IF (Type .EQ. 72 .OR. Type .EQ. 73 .OR. Type .EQ. 74 .OR.
     +            Type .EQ. 75 .OR. Type .EQ. 76) THEN
            If (Type .EQ. 72) JALABEL  = denstyper_aa(15)
            If (Type .EQ. 73) JALABEL  = denstyper_bb(15)
            If (Type .EQ. 74) JALABEL  = denstypel_aa(15)
            If (Type .EQ. 75) JALABEL  = denstypel_bb(15)
            If (Type .EQ. 76) JALABEL  = denstype_tl (15)

         ELSE IF (Type .EQ. 77 .OR. Type .EQ. 78 .OR. Type .EQ. 79 .OR.
     +            Type .EQ. 80 .OR. Type .EQ. 81) THEN
            If (Type .EQ. 77) JALABEL  = denstyper_aa(16)
            If (Type .EQ. 78) JALABEL  = denstyper_bb(16)
            If (Type .EQ. 79) JALABEL  = denstypel_aa(16)
            If (Type .EQ. 80) JALABEL  = denstypel_bb(16)
            If (Type .EQ. 81) JALABEL  = denstype_tl (16)

         ELSE IF (Type .EQ. 82 .OR. Type .EQ. 83 .OR. Type .EQ. 84 .OR.
     +            Type .EQ. 85 .OR. Type .EQ. 86) THEN
            If (Type .EQ. 82) JALABEL  = denstyper_aa(17)
            If (Type .EQ. 83) JALABEL  = denstyper_bb(17)
            If (Type .EQ. 84) JALABEL  = denstypel_aa(17)
            If (Type .EQ. 85) JALABEL  = denstypel_bb(17)
            If (Type .EQ. 86) JALABEL  = denstype_tl (17)

         ELSE IF (Type .EQ. 87 .OR. Type .EQ. 88 .OR. Type .EQ. 89 .OR.
     +            Type .EQ. 90 .OR. Type .EQ. 91) THEN
            If (Type .EQ. 87) JALABEL  = denstyper_aa(18)
            If (Type .EQ. 88) JALABEL  = denstyper_bb(18)
            If (Type .EQ. 89) JALABEL  = denstypel_aa(18)
            If (Type .EQ. 90) JALABEL  = denstypel_bb(18)
            If (Type .EQ. 91) JALABEL  = denstype_tl (18)

         ELSE IF (Type .EQ. 92 .OR. Type .EQ. 93 .OR. Type .EQ. 94 .OR.
     +            Type .EQ. 95 .OR. Type .EQ. 96) THEN
            If (Type .EQ. 92) JALABEL  = denstyper_aa(19)
            If (Type .EQ. 93) JALABEL  = denstyper_bb(19)
            If (Type .EQ. 94) JALABEL  = denstypel_aa(19)
            If (Type .EQ. 95) JALABEL  = denstypel_bb(19)
            If (Type .EQ. 96) JALABEL  = denstype_tl (19)
         ELSE IF (Type .EQ. 97 .OR. Type .EQ. 98 .OR. Type .EQ. 99 .OR.
     +            Type .EQ. 100 .OR. Type .EQ. 101) THEN
            If (Type .EQ. 97)  JALABEL  = denstyper_aa(20)
            If (Type .EQ. 98)  JALABEL  = denstyper_bb(20)
            If (Type .EQ. 99)  JALABEL  = denstypel_aa(20)
            If (Type .EQ. 100)  JALABEL  = denstypel_bb(20)
            If (Type .EQ. 101) JALABEL  = denstype_tl (20)
         ELSE IF (Type .EQ. 102 .OR. Type .EQ. 103 .OR. Type .EQ. 104
     +            .OR. Type .EQ. 105 .OR. Type .EQ. 106) THEN
            If (Type .EQ. 102)  JALABEL  = denstyper_aa(21)
            If (Type .EQ. 103)  JALABEL  = denstyper_bb(21)
            If (Type .EQ. 104)  JALABEL  = denstypel_aa(21)
            If (Type .EQ. 105)  JALABEL  = denstypel_bb(21)
            If (Type .EQ. 106)  JALABEL  = denstype_tl (21)
         ELSE IF (Type .EQ. 107 .OR. Type .EQ. 108 .OR. Type .EQ. 109
     +            .OR. Type .EQ. 110 .OR. Type .EQ. 111) THEN
            If (Type .EQ. 107)  JALABEL  = denstyper_aa(22)
            If (Type .EQ. 108)  JALABEL  = denstyper_bb(22)
            If (Type .EQ. 109)  JALABEL  = denstypel_aa(22)
            If (Type .EQ. 110)  JALABEL  = denstypel_bb(22)
            If (Type .EQ. 111)  JALABEL  = denstype_tl (22)

         ELSE IF (Type .EQ. 112 .OR. Type .EQ. 113 .OR. Type .EQ. 114
     +            .OR. Type .EQ. 115 .OR. Type .EQ. 116) THEN
            If (Type .EQ. 112)  JALABEL  = denstyper_aa(23)
            If (Type .EQ. 113)  JALABEL  = denstyper_bb(23)
            If (Type .EQ. 114)  JALABEL  = denstypel_aa(23)
            If (Type .EQ. 115)  JALABEL  = denstypel_bb(23)
            If (Type .EQ. 116)  JALABEL  = denstype_tl (23)
         ELSE IF (Type .EQ. 117 .OR. Type .EQ. 118 .OR. Type .EQ. 119
     +            .OR. Type .EQ. 120 .OR. Type .EQ. 121) THEN
            If (Type .EQ. 117)  JALABEL  = denstyper_aa(24)
            If (Type .EQ. 118)  JALABEL  = denstyper_bb(24)
            If (Type .EQ. 119)  JALABEL  = denstypel_aa(24)
            If (Type .EQ. 120)  JALABEL  = denstypel_bb(24)
            If (Type .EQ. 121)  JALABEL  = denstype_tl (24)
         ELSE IF (type .EQ. -3) THEN ! Special for now 5-2011
            FLAG = 1 ! anything but density
            JALABEL='AOOVRLAP'
            WRITE (*,*)
            WRITE (*,*) ' Writing atomic overlap integrals to JOBARC '
            WRITE (*,*)
        ELSE IF (type .EQ.-5) THEN 
C
C This block is reserved to dump the record that flag the states with
C spin multplicity. The triplets carry 2 and singlet carry 1.
C
         IEOM_ROOTS = EOM_ROOTS*2 
         CALL  DPUTREC (1, 'JOBARC', "SPINMULT", IEOM_ROOTS,DENS)

         ELSE
            IF (.NOT.(TYPE .EQ. -1 .OR. TYPE .EQ. -2 .OR. TYPE .EQ. -5))
     +      THEN
                WRITE (*,*)
                WRITE (*,*) ' No density for this integer label! '
                RETURN
            ENDIF 

         ENDIF

         IF (Type .NE. -1) Then
             NDIM = NDIM1 * NDIM2
             CALL  DPUTREC (1, 'JOBARC', JALABEL, NDIM,DENS)
         ENDIF 

C The following block of code is no longer necessary. All the
C reordering and rescaling is transfered to a2proc routines. 
C From here write everyting as they come out from ACES III.
C Ajith Perera. 01/2014.
C
CSSS         IF (intpkg .eq. flocke_package) THEN
CSSS
CSSS             do nu = 1, ndim2
CSSS             do mu = 1, ndim1
CSSS                scr(mu,nu) = dens(mu,nu)
CSSS             enddo
CSSS             enddo
CSSS            
CSSS             IF (FLAG .EQ. 0) THEN
CSSS
CSSS                do nu = 1, ndim2
CSSS                do mu = 1, ndim1
CSSS                   ierd1 = ibuf (ierdind+mu-1)
CSSS                   ierd2 = ibuf (ierdind+nu-1)
CSSS                   dscale = dbuf(iscale_fac+ierd1-1) * 
CSSS     +                      dbuf(iscale_fac+ierd2-1)
CSSS                   dens(ierd1,ierd2) = scr(mu,nu) / dscale
CSSS                enddo
CSSS                enddo
CSSS
CSSS             ELSE IF (FLAG .EQ. 1) THEN
CSSS
CSSS                do nu = 1, ndim2
CSSS                do mu = 1, ndim1
CSSS                   ierd1 = ibuf (ierdind+mu-1)
CSSS                   ierd2 = ibuf (ierdind+nu-1)
CSSS                   dscale = dbuf(iscale_fac+ierd1-1) * 
CSSS     +                      dbuf(iscale_fac+ierd2-1)
CSSS                   dens(ierd1,ierd2) = scr(mu,nu) * dscale
CSSS                enddo
CSSS                enddo
CSSS
CSSS             ENDIF
CSSS
CSSS            do nu = 1, ndim2
CSSS
CSSS               do mu = 1, ndim1
CSSS                  scr (mu) = dens (mu,nu)
CSSS               enddo
CSSS
CSSS               do mu = 1, ndim1
CSSS                  ierd = ibuf (ierdind+mu-1)
CSSS                  dens (ierd,nu) = scr(mu) / dbuf (iscale_fac+ierd-1)
CSSS               enddo
CSSS
CSSS            enddo
CSSS
CSSS            do mu = 1, ndim1
CSSS
CSSS               do nu = 1, ndim2
CSSS                  scr(nu) = dens (mu,nu)
CSSS               enddo
CSSS
CSSS               do nu = 1, ndim2
CSSS                  ierd = ibuf (ierdind+nu-1)
CSSS                  dens (mu,ierd) = scr (nu) / dbuf (iscale_fac+ierd-1)
CSSS               enddo
CSSS
CSSS            enddo
CSSS
CSSS         ENDIF
CSSS         IF (FLAG .EQ. 0 .OR. FLAG .EQ. 1) 
CSSS     +   CALL  DESCRAMBLE_MATRIX
CSSS     +
CSSS     +               (DENS,NDIM1,NFPS,IXSHELLS,NSHELLS)
CSSS     +
CSSS     +
CSSS         Write(6,*) "Inside dens2jobarc, ndim1, ndim2", ndim1, ndim2
C

      ENDIF
C
C Block for the RHF wave functions. Write both Alpha and Beta densities
C
      ELSE 

      FLAG = -1 
      IF (me .EQ. 0) THEN
 
         IF (type .eq. -1) CALL DPUTREC (1, 'JOBARC',"CIS     ",
     &                                   1,SCR)
         IF (type .eq. -2) THEN
             NROOTS = EOM_ROOTS
             CALL DGETREC(1, 'JOBAC', 'EOM_EXES', NROOTS, EOM_EXE_S)
             CALL ACES_JA_TRUNCATE('MARKER  ',1)
             CALL DPUTREC (1, 'JOBARC',"EOM     ", 1, SCR)
             CALL DpUTREC(1, 'JOBAC', 'EOM_EXES', NROOTS, EOM_EXE_S)
             RETURN
         ENDIF 

         IF      (type .EQ.  0) THEN
            JALABEL='HFDENSTY'
         ELSE IF (type .EQ.  1) THEN
            JALABEL='TDENSITY'
         ELSE IF (Type .EQ. 2 .OR. Type .EQ. 3 .OR. Type .EQ. 4) THEN
            If (Type .EQ. 2) JALABEL  = denstyper(1)
            If (Type .EQ. 3) JALABEL  = denstypel(1)
            If (Type .EQ. 4) JALABEL  = denstype_tl (1)
         ELSE IF (Type .EQ. 5 .OR. Type .EQ. 6 .OR. Type .EQ. 7) THEN
            If (Type .EQ. 5) JALABEL  = denstyper(2)
            If (Type .EQ. 6) JALABEL  = denstypel(2)
            If (Type .EQ. 7) JALABEL  = denstype_tl (2)
         ELSE IF (Type .EQ. 8 .OR. Type .EQ. 9 .OR. Type .EQ. 10) THEN
            If (Type .EQ. 8)  JALABEL  = denstyper(3)
            If (Type .EQ. 9)  JALABEL  = denstypel(3)
            If (Type .EQ. 10) JALABEL  = denstype_tl (3)
         ELSE IF (Type .EQ. 11 .OR. Type .EQ. 12 .OR. Type .EQ. 13) THEN
            If (Type .EQ. 11) JALABEL  = denstyper(4)
            If (Type .EQ. 12) JALABEL  = denstypel(4)
            If (Type .EQ. 13) JALABEL  = denstype_tl (4)
         ELSE IF (Type .EQ. 14 .OR. Type .EQ. 15 .OR. Type .EQ. 16) THEN
            If (Type .EQ. 14) JALABEL  = denstyper(5)
            If (Type .EQ. 15) JALABEL  = denstypel(5)
            If (Type .EQ. 16) JALABEL  = denstype_tl (5)
         ELSE IF (Type .EQ. 17 .OR. Type .EQ. 18 .OR. Type .EQ. 19) THEN
            If (Type .EQ. 17) JALABEL  = denstyper(6)
            If (Type .EQ. 18) JALABEL  = denstypel(6)
            If (Type .EQ. 19) JALABEL  = denstype_tl (6)
         ELSE IF (Type .EQ. 20 .OR. Type .EQ. 21 .OR. Type .EQ. 22) THEN
            If (Type .EQ. 20) JALABEL  = denstyper(7)
            If (Type .EQ. 21) JALABEL  = denstypel(7)
            If (Type .EQ. 22) JALABEL  = denstype_tl (7)
         ELSE IF (Type .EQ. 23 .OR. Type .EQ. 24 .OR. Type .EQ. 25) THEN
            If (Type .EQ. 23) JALABEL  = denstyper(8)
            If (Type .EQ. 24) JALABEL  = denstypel(8)
            If (Type .EQ. 25) JALABEL  = denstype_tl (8)
         ELSE IF (Type .EQ. 26 .OR. Type .EQ. 27 .OR. Type .EQ. 28) THEN
            If (Type .EQ. 26) JALABEL  = denstyper(9)
            If (Type .EQ. 27) JALABEL  = denstypel(9)
            If (Type .EQ. 28) JALABEL  = denstype_tl (9)
         ELSE IF (Type .EQ. 29 .OR. Type .EQ. 30 .OR. Type .EQ. 31) THEN
            If (Type .EQ. 29) JALABEL  = denstyper(10)
            If (Type .EQ. 30) JALABEL  = denstypel(10)
            If (Type .EQ. 31) JALABEL  = denstype_tl (10)
         ELSE IF (Type .EQ. 32 .OR. Type .EQ. 33 .OR. Type .EQ. 34) THEN
            If (Type .EQ. 32) JALABEL  = denstyper(11)
            If (Type .EQ. 33) JALABEL  = denstypel(11)
            If (Type .EQ. 34) JALABEL  = denstype_tl (11)
         ELSE IF (Type .EQ. 35 .OR. Type .EQ. 36 .OR. Type .EQ. 37) THEN
            If (Type .EQ. 35) JALABEL  = denstyper(12)
            If (Type .EQ. 36) JALABEL  = denstypel(12)
            If (Type .EQ. 37) JALABEL  = denstype_tl (12)
         ELSE IF (Type .EQ. 38 .OR. Type .EQ. 39 .OR. Type .EQ. 40) THEN
            If (Type .EQ. 38) JALABEL  = denstyper(13)
            If (Type .EQ. 39) JALABEL  = denstypel(13)
            If (Type .EQ. 40) JALABEL  = denstype_tl (13)
         ELSE IF (Type .EQ. 41 .OR. Type .EQ. 42 .OR. Type .EQ. 43) THEN
            If (Type .EQ. 41) JALABEL  = denstyper(14)
            If (Type .EQ. 42) JALABEL  = denstypel(14)
            If (Type .EQ. 43) JALABEL  = denstype_tl (14)
         ELSE IF (Type .EQ. 44 .OR. Type .EQ. 45 .OR. Type .EQ. 46) THEN
            If (Type .EQ. 44) JALABEL  = denstyper(15)
            If (Type .EQ. 45) JALABEL  = denstypel(15)
            If (Type .EQ. 46) JALABEL  = denstype_tl (15)
         ELSE IF (Type .EQ. 47 .OR. Type .EQ. 48 .OR. Type .EQ. 49) THEN
            If (Type .EQ. 47) JALABEL  = denstyper(16)
            If (Type .EQ. 48) JALABEL  = denstypel(16)
            If (Type .EQ. 49) JALABEL  = denstype_tl (16)
         ELSE IF (Type .EQ. 50 .OR. Type .EQ. 51 .OR. Type .EQ. 52) THEN
            If (Type .EQ. 50) JALABEL  = denstyper(17)
            If (Type .EQ. 51) JALABEL  = denstypel(17)
            If (Type .EQ. 52) JALABEL  = denstype_tl (17)
         ELSE IF (Type .EQ. 53 .OR. Type .EQ. 54 .OR. Type .EQ. 55) THEN
            If (Type .EQ. 53) JALABEL  = denstyper(18)
            If (Type .EQ. 54) JALABEL  = denstypel(18)
            If (Type .EQ. 55) JALABEL  = denstype_tl (18)
         ELSE IF (Type .EQ. 56 .OR. Type .EQ. 57 .OR. Type .EQ. 58) THEN
            If (Type .EQ. 56) JALABEL  = denstyper(19)
            If (Type .EQ. 57) JALABEL  = denstypel(19)
            If (Type .EQ. 58) JALABEL  = denstype_tl (19)
         ELSE IF (Type .EQ. 59 .OR. Type .EQ. 60 .OR. Type .EQ. 61) THEN
            If (Type .EQ. 59) JALABEL  = denstyper(20)
            If (Type .EQ. 60) JALABEL  = denstypel(20)
            If (Type .EQ. 61) JALABEL  = denstype_tl (20)
         ELSE IF (Type .EQ. 62 .OR. Type .EQ. 63 .OR. Type .EQ. 64) THEN
            If (Type .EQ. 62) JALABEL  = denstyper(21)
            If (Type .EQ. 63) JALABEL  = denstypel(21)
            If (Type .EQ. 64) JALABEL  = denstype_tl (21)
         ELSE IF (Type .EQ. 65 .OR. Type .EQ. 66 .OR. Type .EQ. 67) THEN
            If (Type .EQ. 65) JALABEL  = denstyper(22)
            If (Type .EQ. 66) JALABEL  = denstypel(22)
            If (Type .EQ. 67) JALABEL  = denstype_tl (22)
         ELSE IF (Type .EQ. 68 .OR. Type .EQ. 69 .OR. Type .EQ. 70) THEN
            If (Type .EQ. 68) JALABEL  = denstyper(23)
            If (Type .EQ. 69) JALABEL  = denstypel(23)
            If (Type .EQ. 70) JALABEL  = denstype_tl (23)
         ELSE IF (Type .EQ. 71 .OR. Type .EQ. 72 .OR. Type .EQ. 73) THEN
            If (Type .EQ. 71) JALABEL  = denstyper(24)
            If (Type .EQ. 72) JALABEL  = denstypel(24)
            If (Type .EQ. 73) JALABEL  = denstype_tl (24)

         ELSE IF (type .EQ. -3) THEN ! Special for now 5-2011
            FLAG = 1 ! anything but density
            JALABEL='AOOVRLAP'
            WRITE (*,*)
            WRITE (*,*) ' Writing atomic overlap integrals to JOBARC '
            WRITE (*,*)

        ELSE IF (type .EQ.-5) THEN 
C
C This block is reserved to dump the record that flag the states with
C spin multplicity. The triplets carry 2 and singlet carry 1.
C
         IEOM_ROOTS = EOM_ROOTS*2 
         CALL  DPUTREC (1, 'JOBARC', "SPINMULT", IEOM_ROOTS,DENS)

         ELSE 
            IF (.NOT.(TYPE .EQ. -1 .OR. TYPE .EQ. -2 .OR. TYPE .EQ. -5))
     +      THEN 
            WRITE (*,*)
            WRITE (*,*) ' No density for this integer label! '
            ENDIF 
           
            RETURN

         ENDIF 

         IF (TYPE .NE. -1) THEN 
            NDIM = NDIM1 * NDIM2
            CALL  DPUTREC (1, 'JOBARC', JALABEL, NDIM,DENS)
         ENDIF 

       ENDIF 
      ENDIF 

      RETURN
      END
