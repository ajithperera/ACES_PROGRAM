










      subroutine write_vibron_input(freq, numvib, scr)
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
      integer ione, i_sided, ijunk, irp, nmode2, npoint, irpi, irpj,
     $     i, j, k, nirpf, iunit, perline, ioff, ndim, ndim2,
     $     iref, ipnt, igrd, ifcm, ifreq, itop,
     $     ioff1, ioff2, ioff3, ioff4, numvib, icase,
     $     ncouple, nsymm
      character*4 doit, fnameirp
      character*8 labelf(14)
      character*50 string80, string
      double precision
     $     scr(*), stpsiz, freq(numvib), fact, e0, fred, two, x

      two = 2.0d0
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
      stpsiz=dfloat(iflags(57))*10.0d-5
      call getrec(1, 'JOBARC', 'FREQS_0 ', numvib*iintfp,
     $     freq)
c     
      nmode = 0
      do irp = 1, nirrepf
         nmode = nmode + numirpf(irp)
      enddo
c     
c     All results will be written to file 'Vibronic_Input'
c     
      iunit = 172
      perline = 8
      open(UNIT=iunit, FILE='Vibron_Input',
     $     FORM='FORMATTED')
      rewind(iunit)
      nsymm = numirpf(1)
      ncouple = nmode - nsymm
c     
c     write common info for this calculation
c     
 922  format(i6, t20, ' ! units in eV')
 923  format('xxx', t20, ' ! Number of electronic states')
 924  format(i6, t20, ' ! Number of symmetric normal modes ')
 925  format(i6, t20, ' ! Number of non-A1 normal modes ')
 926  format(i6, t20, ' ! Number of Lanczos iterations ')
 927  format(50i4)
 928  format(F8.3)
 930  format(t2,A50)
c     
      i = 0
      write(iunit,922) i
      write(iunit, 923) 
      write(iunit,924) nsymm
      write(iunit,925) ncouple
      i = 500
      write(iunit,926) i
      do i = 1, nmode
         nquanta(i) = 0
      enddo
      write(iunit, 927) (nquanta(i), i = 1, nsymm)
      write(iunit,704)
     $     (freq(i) /8065.5d0, i=1, nsymm)
 704  format(50F12.6)
c     
      write(iunit,927) (nquanta(i), i = 1, ncouple)
      write(iunit,704) (freq(i) /8065.5d0, i=nsymm+1, nmode)
c
 660  format(20A)
      x = 0.1
      write(iunit,928) x
      x = 3.0
      write(iunit,928) x
      string='auto'
      write(iunit,*) string(1:4)
      write(iunit,*)
      string = '*vibron'
      write(iunit,660) string(1:7)
      string = 'short_title=cp'
      write(iunit,*) string(1:14)
      string = 'refine_coupling=on'
      write(iunit,*) string(1:18)
      string = 'cubic=on'
      write(iunit,*) string(1:8)
      string =  'quartic=on'
      write(iunit,*) string(1:10)
      string =  '*solve_coupling'
      write(iunit,660) string(1:16)
      string = 'purify_couplings=on'
      write(iunit,*) string(1:19)
      string = '*end'
      write(iunit,660) string(1:4)
      write(iunit,*)
c     
      close(unit=iunit, status='KEEP')
c     
      write(6,*) ' All done write_vibron_input '
c     
      return
      end
