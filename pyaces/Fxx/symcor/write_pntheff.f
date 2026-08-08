










      subroutine write_pntheff(scr)
c     
c     In this subroutine the PNTHEFF file is written for each electronic structure method.
c     This file can be processed further by the program vibron (using option refine_coupling)
c     
      implicit none








c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler























cYAU - ACES3 stuff . . . we hope - #include <aces.par>




c
c global variables in the program vibron are set.
c
c maximum dimensions
c
      integer max_mode, max_state, max_quanta,max_irrep,
     $    max_tile
      parameter (max_mode=50)
      parameter (max_state=30)
      parameter (max_quanta=10000)
      parameter (max_irrep=20)
      parameter (max_tile=30)
c
      double precision linear_couplings(max_mode, max_state, max_state)
      double precision quadratic_couplings(max_mode, max_mode,
     $    max_state, max_state)
      double precision cubic_couplings(max_mode, max_state, max_state)
      double precision quartic_couplings(max_mode, max_state, max_state)
      double precision frequencies(max_mode, max_state)
      double precision e_electronic(max_state, max_state)
      double precision e_weight(max_state)
      logical diag_quadratic, off_quadratic, linear, auto, vecfol,
     $    use_rfa, do_franck_condon, interpret, tmom_derivative,
     $    unitary, sequential, symmetric, complex, cubic, quartic,
     $    lcdspectrum
c
c The number of irreps in the full symmetry group and the number of
c modes per block.
c
      integer nirrepf, numirpf(max_irrep)
c
c The number of normal modes considered
c
      integer nmode
c
c The number of electronic states
c
      integer nstate, nstate_orig
c
      integer nsym
c
c The number of symmetric modes
c
      integer ncoupling
c
c The number of coupling modes
c
c     The number of quanta per mode.
c
      integer nquanta(max_mode)
c
c The number of electronic absorbing states
c
      integer nstate_abs
c
c The number of quanta in the absorbing states
c
      integer nquanta_abs(max_mode)
c
c the maximum number of setup state in diagonalization algorithms.
c
      integer max_setup
c
c transition moments, and their gradients along normal modes. 
c They couple initial and final states.
c
      double precision
     $    tmom0(3, max_state, max_state),
     $    grad_tmom(3, max_mode, max_state, max_state),
     $    tmag0(3, max_state, max_state),
     $    grad_tmag(3, max_mode, max_state, max_state)
c
c indicates either real or complex solutions in transformation
c and subsequent lanczos calculation
c
      integer n_complex
c
c To perform more effective H C multiplications and to reduce the 
c dimensions of the vector, we partition the nomal modes into tiles. For each tile
c we assign a maximum excitation level
c
      integer
     $    dim_tile(max_tile), level_tile(max_tile), tile_low(max_tile),
     $    tile_high(max_tile), ntiles, off_annihilation(max_tile),
     $    off_creation(max_tile), off_quanta(max_tile), n_vibron_stub,
     $    n_graph_stub
c
      character*30 vibron_filename, vibron_stub,
     $    graph_stub, graph_filename
c
      common /vibron_int/ nmode, nstate, nsym, ncoupling, nquanta,
     $     nirrepf, numirpf, nstate_abs, nquanta_abs, max_setup,
     $    n_complex, nstate_orig, n_vibron_stub, n_graph_stub
      common /vibron_log/ diag_quadratic, off_quadratic, linear, auto,
     $    vecfol, use_rfa, do_franck_condon, interpret, tmom_derivative,
     $    unitary, sequential, symmetric, complex, cubic, quartic,
     $    lcdspectrum
      common /vibron_real/ linear_couplings, frequencies, e_electronic,
     $    quadratic_couplings, tmom0, grad_tmom, e_weight,
     $    cubic_couplings, quartic_couplings, tmag0, grad_tmag
      common /vibron_tile/
     $    dim_tile, level_tile, tile_low,
     $    tile_high, ntiles, off_annihilation,
     $    off_creation, off_quanta
      common /vibron_char/ vibron_stub, graph_stub
c
c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End: c


c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end



      integer ncases
      parameter (ncases=10)
      character*8 ndimheff(ncases)
      character*8 nameheff(ncases)
      character*8 namenvec(ncases)
      character*8 namehvec(ncases)
      character*8 namehdiab(ncases)
      character*8 nametran(ncases)
      character*8 nametmom(ncases)
      character*8 refheff(ncases)
      character*8 reftmom(ncases)
      character*8 pntheff(ncases)
      character*8 pnttmom(ncases)
c
      data ndimheff /'NEFF_IP2', 'NEFF_EA2', 'NEFF_EE1', 'NEFF_EE3',
     $                'NEFFDIP1', 'NEFFDIP3', 'NEFFDEA1', 'NEFFDEA3',
     $                'NEFFDEE1', 'NEFFDEE3' /
      data nameheff /'HEFF_IP2', 'HEFF_EA2', 'HEFF_EE1', 'HEFF_EE3',
     $                'HEFFDIP1', 'HEFFDIP3', 'HEFFDEA1', 'HEFFDEA3',
     $                'HEFFDEE1', 'HEFFDEE3' /
      data namenvec /'NVEC_IP2', 'NVEC_EA2', 'NVEC_EE1', 'NVEC_EE3',
     $     'NVECDIP1', 'NVECIP3', 'NVECDEA1', 'NVECDEA3',
     $     'NVECDEE1', 'NVECDEE3' /
      data namehvec /'HVEC_IP2', 'HVEC_EA2', 'HVEC_EE1', 'HVEC_EE3',
     $     'HVECDIP1', 'HVECDIP3', 'HVECDEA1', 'HVECDEA3',
     $     'HVECDEE1', 'HVECDEE3' /
c     
      data namehdiab /'HDIA_IP2', 'HDIA_EA2', 'HDIA_EE1', 'HDIA_EE3',
     $     'HDIADIP1', 'HDIADIP3', 'HDIADEA1', 'HDIADEA3',
     $     'HDIADEE1', 'HDIADEE3' /
      data nametran /'TRAN_IP2', 'TRAN_EA2', 'TRAN_EE1', 'TRAN_EE3',
     $     'TRANDIP1', 'TRANDIP3', 'TRANDEA1', 'TRANDEA3',
     $     'TRANDEE1', 'TRANDEE3' /
      data refheff /'HREF_IP2', 'HREF_EA2', 'HREF_EE1', 'HREF_EE3',
     $                'HREFDIP1', 'HREFDIP3', 'HREFDEA1', 'HREFDEA3',
     $                'HREFDEE1', 'HREFDEE3' /
      data pntheff /'HPNT_IP2', 'HPNT_EA2', 'HPNT_EE1', 'HPNT_EE3',
     $                'HPNTDIP1', 'HPNTDIP3', 'HPNTDEA1', 'HPNTDEA3',
     $                'HPNTDEE1', 'HPNTDEE3' /
      data nametmom /'TMOM_IP2', 'TMOM_EA2', 'TMOM_EE1', 'TMOM_EE3',
     $                'TMOMDIP1', 'TMOMDIP3', 'TMOMDEA1', 'TMOMDEA3',
     $                'TMOMDEE1', 'TMOMDEE3' /
      data reftmom /'TREF_IP2', 'TREF_EA2', 'TREF_EE1', 'TREF_EE3',
     $                'TREFDIP1', 'TREFDIP3', 'TREFDEA1', 'TREFDEA3',
     $                'TREFDEE1', 'TREFDEE3' /
      data pnttmom /'TPNT_IP2', 'TPNT_EA2', 'TPNT_EE1', 'TPNT_EE3',
     $                'TPNTDIP1', 'TPNTDIP3', 'TPNTDEA1', 'TPNTDEA3',
     $                'TPNTDEE1', 'TPNTDEE3' /
c
c Effective hamiltonians could be stored on JOBARC for a variety of methods.
c IP / EA / EE / DIP / DEA /DEE indicate the IP/EA-eomcc methods or the STEOM variants. DEE is eomee
c The last digit 1, 2 or 3 indicates the multiplicity of the manifold.
c
c The mrcc program writes the # of states included (ndimheff), and it writes the
c records 'nameheff' (ndim*ndim) dimensional and a transition moment record nametmom,
c of dimension 3*ndim. they contain the effective hamiltonian and the tmoms in 
c the diabatic basis.
c
c The rest of these records is internal to symcor. The reference values (refheff, reftmom)
c need to be
c stored and the finite difference records for all geometry points (pntheff, pnttmom)
c
c finally we will calculate heffgrd, tmomgrd, and possibly hefffcm and tmomfcm later on.
c These will be written to files HGRD_method and TGRD_method etc.
c






































































































































































































c This common block contains the IFLAGS and IFLAGS2 arrays for JODA ROUTINES
c ONLY! The reason is that it contains both arrays back-to-back. If the
c preprocessor define MONSTER_FLAGS is set, then the arrays are compressed
c into one large (currently) 600 element long array; otherwise, they are
c split into IFLAGS(100) and IFLAGS2(500).

c iflags(100)  ASVs reserved for Stanton, Gauss, and Co.
c              (Our code is already irrevocably split, why bother anymore?)
c iflags2(500) ASVs for everyone else

      integer        iflags(100), iflags2(500)
      common /flags/ iflags,      iflags2
      save   /flags/




c     
      integer i_sided, npoint, npoint2, nmode2, ijunk, ione,
     $     irp, irpi, irpj, iunit2, ndim, ndim2,
     $     iref, ipnt, ipnt2, i, nspecial, nrest, ioff, ioff2,
     $     icount, icase
      double precision stpsiz, scr(*)
      character*4 doit
      character*6 ndimheff0, nameheff0
c     
 888  format(A50,A20)
c     
c     get some of the basic info from Jobarc / Flags
c     
      ione = 1
      i_sided = iflags2(157)
      CALL GETREC(20,'JOBARC','NUMPOINT',IONE,npoint)
c     
      CALL GETREC(-1,'JOBARC','DANGERUS',IONE,IJUNK)
      if(iflags(79).eq.0.and.ijunk.eq.0)doit='FULL'
      if(iflags(79).eq.1.or.ijunk.eq.1)doit='COMP'
      CALL GETREC(20,'JOBARC',DOIT//'NIRX',IONE,nirrepf)
      call getrec(20, 'JOBARC', 'NMODES_F', nirrepf, numirpf)
c
      if (IFLAGS(57) .lt. 0) then
         write(6,*)  ' We used displacements along ',
     $        'dimensionless reduced coordinates.'
         stpsiz = DFLOAT(IFLAGS(57)) / 1000.0d0
         write(6,511) stpsiz
 511     format(T3,' Elementary stepsize ', F8.3, ' Dimensionless ')
      else
         STPSIZ=DFLOAT(IFLAGS(57))*10.0D-5
         write(6,*) ' Use Normal mode displacements '
         WRITE(6,500)STPSIZ
 500     FORMAT(T3,'Step size ',F8.5,' amu**(1/2) * bohr.')
      endif
c     
      write(6,*) ' @write_pntheff : normal mode information'
      write(6,800) nirrepf
 800  format('  Number of symmetries : ', I8)
      write(6,*) ' Number of normal modes per symmetry block '
      write(6,*) '   IRREP       # of Modes '
      do irp = 1, nirrepf
         write(6,801) irp, numirpf(irp)
      enddo
 801  format(2I8)
c     
      write(6,*) 
c     
      nmode = 0
      do irp = 1, nirrepf
         nmode = nmode + numirpf(irp)
      enddo
      write(6,802) nmode
 802  format( '  Total number of normal modes ', I8)
c     
      nmode2 = 0
      if (i_sided .ge. 3 .and. i_sided .ne. 5) then
         do irp = 1, nirrepf
            nmode2 = nmode2 + numirpf(irp)*(numirpf(irp)-1) / 2
         enddo
      endif
      if (i_sided .eq. 4 .or. i_sided .eq. 7) then
         do irpi = 1, nirrepf
            do irpj = 1, irpi - 1
               nmode2 = nmode2 + numirpf(irpi) * numirpf(irpj)
            enddo
         enddo
      endif
c     
c     calculate full set of points, including points that we know by symmetry.
c     
      npoint2 = nmode
      if (i_sided .ge. 2) then
         npoint2 = npoint2 + nmode
      endif
      if (i_sided .ge. 3 .and. i_sided .ne. 5) then
         npoint2 = npoint2 + 2*nmode2
      endif
      if (i_sided .ge. 5) then
         npoint2 = npoint2 + 2 * nmode
      endif
c
      if (i_sided .eq. 8) then
         npoint2 = npoint
         nmode2 = npoint - 2 * nmode
      elseif (i_sided .ge. 9) then
         npoint2 = npoint
         nmode2 = npoint - 4 * nmode
      endif
c
      write(6,804) nmode2
 804  format(' Total number of quadratic displacements ', I8)
c     
 806     format(' Unique number of points: ', I8)
         write(6,806) npoint
         write(6,807) npoint2
 807     format( ' Expanded number of points: ', I8)

      nspecial = nmode
      nrest = nmode - nspecial

      iunit2 = 185
      open (unit=iunit2,file='PNTHEFF',form='formatted')
      rewind (iunit2)
c     
c write parent energy array in usual format
c
c
         ndimheff0 = 'NEFF_0'
         nameheff0 = 'HEFF_0'
         ndim = 1
         write(6,805) ndimheff0, ndim
c     
         ndim2 = ndim * ndim
         iref = 1
         call getrec(20,'JOBARC','REFPAR_E',
     $        iintfp, scr(iref))
c     
         write(iunit2,888) ' Heff points for method ', nameheff0
         write(iunit2,888) ' Reference hamiltonian'
         do i = 1, ndim2
            write(iunit2,887) scr(iref+i-1)
         enddo
         write(iunit2,889) i_sided, ndim2, nmode2, npoint2, stpsiz
c     
         ipnt = 1
         ipnt2 = ipnt + ndim2*npoint
         call getrec(20,'JOBARC','PNTENERG',
     $        npoint*iintfp, scr(ipnt))
c     
            do i = 1, ndim2*npoint
               write(iunit2,887) scr(ipnt+i-1)
            enddo
c
c write fake transition moments for convenience in VIBRON
c
c     
c     process transition moments
c     
         iref = 1
c
         call dzero(scr(iref), 3*ndim)
c         call getrec(20,'JOBARC',reftmom(icase),
c     $        ndim*3*iintfp, scr(iref))
c     
         write(iunit2,888) ' TMOM points for method ', nameheff0
         write(iunit2,888) ' Reference transition moments'
         do i = 1, 3*ndim
            write(iunit2,887) scr(iref+i-1)
         enddo
c     
         write(iunit2,889) i_sided, ndim*3, nmode2, npoint2, stpsiz

         ipnt = 1
         ipnt2 = ipnt + ndim*3*npoint
         call dzero(scr(ipnt), 3*ndim*npoint)
c         call getrec(20,'JOBARC',pnttmom(icase),
c     $        ndim*3*npoint*iintfp, scr(ipnt))
c     
         do i = 1, ndim*3*npoint
            write(iunit2,887) scr(ipnt+i-1)
         enddo
c
         write(iunit2,888) ' End point info ', nameheff0
c     
         write(6,808) nameheff0
c
c     done with parent energy
c
      do 10 icase = 1, ncases
c     
         ndim = 0
         call getrec(-1, 'JOBARC', ndimheff(icase), ione, ndim)
         write(6,805) ndimheff(icase), ndim
 805     format(' Process vibronic_coupling ', A12, /,
     $        '   # of electronic states: ', I10)
         if (ndim .eq. 0) goto 10
c     
         ndim2 = ndim * ndim
         iref = 1
         call getrec(20,'JOBARC',refheff(icase),
     $        ndim2*iintfp, scr(iref))
c     
         write(iunit2,888) ' Heff points for method ', nameheff(icase)
         write(iunit2,888) ' Reference hamiltonian'
         do i = 1, ndim2
            write(iunit2,887) scr(iref+i-1)
         enddo
         write(iunit2,889) i_sided, ndim2, nmode2, npoint2, stpsiz
 889     format(4I12, E12.4)
 887     format(F25.12)
c     
         ipnt = 1
         ipnt2 = ipnt + ndim2*npoint
         call getrec(20,'JOBARC',pntheff(icase),
     $        ndim2*npoint*iintfp, scr(ipnt))
c     
         if (npoint .ne. npoint2) then
            ioff = ipnt
            ioff2 = ipnt2
            icount = 1
            write(6,*) ' copy from point , # ', icount, nmode
            call dcopy(nmode*ndim2, scr(ioff), 1,
     $           scr(ioff2), 1)
            ioff2 = ioff2 + nmode*ndim2
            ioff = ioff + nmode*ndim2
            icount = icount + nmode
            if (i_sided .ge. 2) then
               write(6,*) ' copy from point , # ', icount, nspecial
               call dcopy(nspecial*ndim2, scr(ioff), 1,
     $              scr(ioff2), 1)            
               ioff2 = ioff2 + nspecial * ndim2
               ioff = ioff + (nspecial-nmode) * ndim2
               icount = icount + (nspecial -nmode)
               write(6,*) ' copy from point , # ', icount, nrest
               call dcopy(nrest*ndim2, scr(ioff), 1,
     $              scr(ioff2), 1)            
               ioff2 = ioff2 + nrest * ndim2
               ioff = ioff + nmode * ndim2
               icount = icount + nmode
            endif
            if (i_sided .ge. 3 .and. i_sided .ne. 5) then
               write(6,*) ' copy from point , # ', icount, 2*nmode2
               call dcopy(2*nmode2*ndim2, scr(ioff), 1,
     $              scr(ioff2), 1)            
               ioff2 = ioff2 + 2*nmode2 * ndim2
               ioff = ioff + 2*nmode2 * ndim2
               icount = icount + 2*nmode2
            endif
            if (i_sided .ge. 5) then
               write(6,*) ' copy from point , # ', icount, nmode
               call dcopy(nmode*ndim2, scr(ioff), 1,
     $              scr(ioff2), 1)            
               ioff2 = ioff2 + nmode * ndim2
               ioff = ioff + nmode*ndim2
               icount = icount + nmode
               write(6,*) ' copy from point , # ', icount, nspecial
               call dcopy(nspecial*ndim2, scr(ioff), 1,
     $              scr(ioff2), 1)            
               ioff2 = ioff2 + nspecial * ndim2
               ioff = ioff + (nspecial-nmode) * ndim2
               icount = icount + (nspecial - nmode)
               write(6,*) ' copy from point , # ', icount, nrest
               call dcopy(nrest*ndim2, scr(ioff), 1,
     $              scr(ioff2), 1)            
               ioff2 = ioff2 + nrest * ndim2
               ioff = ioff + nmode * ndim2
               icount = icount + nmode
            endif            
            write(6,*) ' final icount, npoint ', icount, npoint
            if (ioff .ne. ipnt + npoint*ndim2) then
               write(6,*) ' Something wrong in calc_couplings a'
               call aces_exit(1)
            endif
            if (ioff2 .ne. ipnt2 + npoint2*ndim2) then
               write(6,*) ' Something wrong in calc_couplings b'
               call aces_exit(1)
            endif
            do i = 1, ndim2*npoint2
               write(iunit2,887) scr(ipnt2+i-1)
            enddo
         else 
            do i = 1, ndim2*npoint
               write(iunit2,887) scr(ipnt+i-1)
            enddo
         endif
c     
c     process transition moments
c     
         iref = 1
c     
         call getrec(20,'JOBARC',reftmom(icase),
     $        ndim*3*iintfp, scr(iref))
c     
         write(iunit2,888) ' TMOM points for method ', nameheff(icase)
         write(iunit2,888) ' Reference transition moments'
         do i = 1, 3*ndim
            write(iunit2,887) scr(iref+i-1)
         enddo
c     
         write(iunit2,889) i_sided, ndim*3, nmode2, npoint2, stpsiz

         ipnt = 1
         ipnt2 = ipnt + ndim*3*npoint
         call getrec(20,'JOBARC',pnttmom(icase),
     $        ndim*3*npoint*iintfp, scr(ipnt))
c     
         if (npoint .ne. npoint2) then
            ioff = ipnt
            ioff2 = ipnt2
            call dcopy(nmode*ndim*3, scr(ioff), 1,
     $           scr(ioff2), 1)
            ioff2 = ioff2 + nmode*ndim*3
            ioff = ioff + nmode*ndim*3
            if (i_sided .ge. 2) then
               call dcopy(nspecial*ndim*3, scr(ioff), 1,
     $              scr(ioff2), 1)            
               ioff2 = ioff2 + nspecial * ndim*3
               ioff = ioff + (nspecial-nmode) * ndim*3
               call dcopy(nrest*ndim*3, scr(ioff), 1,
     $              scr(ioff2), 1)            
               ioff2 = ioff2 + nrest * ndim*3
               ioff = ioff + nmode * ndim*3
            endif
            if (i_sided .ge. 3 .and. i_sided .ne. 5) then
               call dcopy(2*nmode2*ndim*3, scr(ioff), 1,
     $              scr(ioff2), 1)            
               ioff2 = ioff2 + 2*nmode2 * ndim*3
               ioff = ioff + 2*nmode2 * ndim*3
            endif
            if (i_sided .ge. 5) then
               call dcopy(nmode*ndim*3, scr(ioff), 1,
     $              scr(ioff2), 1)            
               ioff2 = ioff2 + nmode * ndim*3
               ioff = ioff + nmode*ndim*3
               call dcopy(nspecial*ndim*3, scr(ioff), 1,
     $              scr(ioff2), 1)            
               ioff2 = ioff2 + nspecial * ndim*3
               ioff = ioff + (nspecial-nmode) * ndim*3
               call dcopy(nrest*ndim*3, scr(ioff), 1,
     $              scr(ioff2), 1)            
               ioff2 = ioff2 + nrest * ndim*3
               ioff = ioff + nmode * ndim*3
            endif            
            if (ioff .ne. ipnt + npoint*ndim*3) then
               write(6,*) ' Something wrong in calc_couplings a'
               call aces_exit(1)
            endif
            if (ioff2 .ne. ipnt2 + npoint2*ndim*3) then
               write(6,*) ' Something wrong in calc_couplings b'
               call aces_exit(1)
            endif
            do i = 1, ndim*3*npoint2
               write(iunit2,887) scr(ipnt2+i-1)
            enddo
         else 
            do i = 1, ndim*3*npoint
               write(iunit2,887) scr(ipnt+i-1)
            enddo
         endif
c     
         write(iunit2,888) ' End point info ', nameheff(icase)
c     
         write(6,808) nameheff(icase)
 808     format(' Finished processing case :', A20)
c     
 10   continue
c     
      write(6,*) ' All done write_pntheff before close'
      close(unit=iunit2, status = 'KEEP')
      write(6,*) ' All done write_pntheff '
c     
      return
      end
