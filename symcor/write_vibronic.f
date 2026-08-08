










      subroutine write_vibronic(freq, numvib, scr)
c     
c     Write an initial vibronic coupling file, which can be processed by VIBRON
c     and turned into a full Vibronic_Coupling file using info on PNTHEFF
c     which is written by write_pntheff
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
     $     ioff1, ioff2, ioff3, ioff4, numvib, icase
      character*4 doit, fnameirp
      character*8 labelf(14)
      character*80 string80
      character*6 nameheff0, ndimheff0
      double precision scr(*), stpsiz, freq(numvib), fact, e0, fred, two
      logical reduced

      two = 2.0d0

c     
c     Documentation on conversion factors:
c     
c     q_i = sqrt(M omega / hbar) * X_i = sqrt(omega * m_u / hbar) * a_0 * Q_i
c     = sqrt(2 pi c nu_i^tilde * m_u / hbar) * a_0 * Q_i
c     = f * sqrt(nu_i^tilde) * Q_i , where nu_i^tilde is the frequency of the
c     parent state normal mode in cm-1, Q_i is the mass weighted normal coordinate, 
c     while q_i = 1/sqrt(2) [ a_i^+ + a_i] is the dimensionless normal coordinate.
c     The conversion factor f evaluates to f=0.0911355
c     
      fred = 0.0911355d0
      fact = 27.2113956324672d0 / 0.0911355d0

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
      if (iflags(57) .lt. 0) then
         reduced = .true.
         stpsiz = dfloat(iflags(57)) /1000.0d0
      else
         reduced = .false.
         stpsiz=dfloat(iflags(57))*10.0d-5
      endif
c     
      write(6,*) ' @write_vibronic : normal mode information'
      write(6,800) nirrepf
 800  format('  Number of symmetries : ', I8)
      write(6,*) ' Number of normal modes per symmetry block '
      write(6,*) '   IRREP       # of Modes '
      do irp = 1, nirrepf
         write(6,801) irp, numirpf(irp)
      enddo
 801  format(2I8)
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
      write(6,804) nmode2
 804  format(' Total number of quadratic displacements ', I8)
c     
c     Up to now the same as in write_pntheff 
c     


      do i = 1, nirrepf
         do k = 1, 8
            labelf(i)(k:k) = ' '
         enddo
      enddo
c     
      open(unit=70, file='IRPNAMES',status='OLD')
      rewind(70)
      read(70,*) string80
      read(70,*) nirpf
      do i = 1, nirpf
         read(70,*) labelf(i)
      enddo
      close(70)
c     
      if (nirpf .ne. nirrepf) then
         write(6,*) ' @write_vibronic, Something suspicious : '
         write(6,*) ' nirpf, nirrepf ', nirpf, nirrepf
c         call aces_exit(1)
      endif
c     
c     All results will be written to file 'Vibronic_Coupling'
c     
      iunit = 172
      perline = 8
      open(UNIT=iunit, FILE='Vibronic_Coupling',
     $     FORM='FORMATTED')
      rewind(iunit)
c     
      write(iunit,888) '  Vibronic Coupling Constants '
      write(iunit,888) '  All units are in eV or dimensionless '
      write(iunit,901)
      write(iunit, 888) ' '
      if (i_sided .eq. 1) then
         write(iunit, 888) ' Single-sided differentiation has been used'
         write(iunit, 888) ' Gradients (linear couplings) only'
      elseif (i_sided .eq. 2) then
         write(iunit,888) ' Double-sided differentiation has been used'
         write(iunit,888) ' Gradients and diagonal elements of Hessian.'
      elseif (i_sided .eq. 3) then
         write(iunit,888) ' Double-sided differentiation has been used'
         write(iunit,888) ' Gradients and totally symmetric Hessian '
      elseif (i_sided .eq. 4) then
         write(iunit,888) ' Double-sided differentiation has been used'
         write(iunit,888) ' Gradients and Complete Hessian '
      elseif (i_sided .eq. 5) then
         write(iunit,888) ' Double-sided differentiation has been used'
         write(iunit,888) ' Gradients and diagonal elements of Hessian',
     $        ' + Cubic/Quartic.'
      elseif (i_sided .eq. 6) then
         write(iunit,888) ' Double-sided differentiation has been used'
         write(iunit,888) ' Gradients and totally symmetric Hessian ',
     $        '+ Cubic/Quartic.'
      elseif (i_sided .eq. 7) then
         write(iunit,888) ' Double-sided differentiation has been used'
         write(iunit,888) ' Gradients and Complete Hessian ',
     $        '+ Cubic/Quartic.'
      endif
      write(iunit,888) ' '
c     
      call getrec(1, 'JOBARC', 'FREQS_0 ', numvib*iintfp,
     $     freq)
      write(iunit,888) ' Parent state Symmetric modes in eV '
      ioff = 1
      do irp = 1, 1
         write(iunit,704)
     $        (freq(ioff+i) /8065.5d0, i=0,numirpf(irp)-1)
         ioff = ioff + numirpf(irp)
      enddo
 704  format(20F12.8)
c     
      write(iunit,888) ' Parent state Non-Symmetric modes in eV '
      write(iunit,704) (freq(i) /8065.5d0, i=ioff, nmode)
      write(iunit,888) ' '
c     
c process ground state (a bit stupid for historical reasons ...
c
c     
         ndim = 1
         ndimheff0 = 'Neff_0'
         nameheff0 = 'Heff_0'
c         call getrec(-1, 'JOBARC', ndimheff(icase), ione, ndim)
         write(6,805) ndimheff0, ndim
c     
         write(iunit,888) ' Vibronic Coupling elements for case ',
     $        nameheff0
         write(iunit,901)
         write(iunit,901)
         write(iunit, 888)
c     
c     write common info for this calculation
c     
         write(iunit,811) i_sided
         write(iunit,812) nirrepf
         write(iunit,813) (numirpf(irp), irp=1,nirrepf)
c     
         call getrec(1, 'JOBARC', 'FREQS_0 ', numvib*iintfp,
     $        freq)
         ioff = 1
         do irp = 1, nirrepf
            if (numirpf(irp) .ne. 0) then
               do i = 1, 4
                  fnameirp(i:i) = labelf(irp)(i:i)
               enddo
               write(iunit,814) irp, numirpf(irp), fnameirp,
     $              (freq(ioff+i), i=0,numirpf(irp)-1)
            endif
            ioff = ioff + numirpf(irp)
         enddo
         write(iunit,888)
c     
c     process heff part
c     
         ndim2 = ndim*ndim
c     
         iref = 1
         ipnt = iref + 1
         igrd = ipnt + npoint
         ifcm = igrd + nmode
         ifreq = ifcm + nmode+nmode2
         itop = ifreq+nmode
c     
         do k = 1, nmode
            scr(ifreq+k-1) = fact / sqrt(abs(freq(k)))
         enddo
c     
c     first treat parent energy
c     
         call getrec(20,'JOBARC','REFPAR_E',
     $        iintfp, scr(iref))
         write(6,*) ' Ground state energy at reference point ',
     $        scr(iref)
c     
         call getrec(20,'JOBARC','PNTENERG',
     $        npoint*iintfp, scr(ipnt))
c     
         if (i_sided .eq. 1) then
            do k = 1, nmode
               scr(igrd+k-1) = 
     $              (scr(ipnt+k-1) - scr(iref)) / stpsiz
            enddo
         else
c     
c     two-sided gradients
c     
            ioff1 = ipnt-1
            ioff2 = ipnt + nmode -1
            do k = 1, nmode
               scr(igrd+k-1) = 
     $              (scr(ioff1+k) -
     $              scr(ioff2+k)) / (2.0d0*stpsiz)
            enddo
c     
c     diagonal elements of Hessian
c     
            do k = 1, nmode
               scr(ifcm+k-1) =
     $              (scr(ioff1+k) +
     $              scr(ioff2+k) - 2.0d0 * scr(iref)) / (stpsiz**2)
            enddo
c     
         endif
c     
         if (i_sided .ge. 3 .and. i_sided .ne. 5) then
c     
c     off-diagonal elements of Hessian (totally symmetric only)
c     
            ioff1 = ipnt - 1
            ioff2 = ioff1 + nmode
            ioff3 = ioff2 + nmode
            ioff4 = ioff3 + nmode2
            k=0
            do irp = 1, nirrepf
               do i = 1, numirpf(irp)
                  do j = i+1, numirpf(irp)
                     k=k+1
                     write(6,*)
     $                    scr(ioff3+k), scr(ioff4+k),
     $                    scr(ioff1+i), scr(ioff1+j),
     $                    scr(ioff2+i), scr(ioff2+j),
     $                    two*scr(iref)
c     
                     scr(ifcm+nmode+k-1) =
     $                    (scr(ioff3+k) + scr(ioff4+k)
     $                    - scr(ioff1+i) - scr(ioff1+j)
     $                    - scr(ioff2+i) - scr(ioff2+j)
     $                    + two*scr(iref) )
                     scr(ifcm+nmode+k-1) = scr(ifcm+nmode+k-1)
     $                    / (stpsiz**2)
                  enddo
               enddo
               ioff1 = ioff1 + numirpf(irp)
               ioff2 = ioff2 + numirpf(irp)
               ioff3 = ioff3 + numirpf(irp)*(numirpf(irp)-1) / 2
               ioff4 = ioff4 + numirpf(irp)*(numirpf(irp)-1) / 2
            enddo
c     
         endif
c     
         write(iunit,888)
         write(iunit,888) ' Gradients parent energy along normal modes'
         write(iunit, 901)
         write(iunit,888)
         do k = 1, nmode
            write(iunit,912) k, scr(igrd+k-1)
         enddo
c
         if (reduced) then
            fact = 27.2114 * 8065.5 
         if (i_sided .ne. 1) then
            write(iunit,888)
            write(iunit,888) ' Diagonal Hessian parent energy '
            write(iunit, 901)
            write(iunit,888)
            do k = 1, nmode
               write(iunit,914) k, 
     $              abs(scr(ifcm+k-1))*fact, 
     $              abs(scr(ifcm+k-1)) * fact / 8065.5
            enddo
         endif
c     
         if (.false. .and. (i_sided .ge. 3 .and. i_sided .ne. 5)) then
            write(6,*) ' off-diagonal elements hessian parent E'
c     
            k=0
            ioff = ifcm+nmode
            do irp = 1, nirrepf
               do i = 1, numirpf(irp)
                  do j = 1, i-1
                     write(iunit,916) irp, i, j,
     $                    abs(scr(ioff+k))*fact,
     $                    abs(scr(ioff+k)) * fact / 8065.5
                     k=k+1
                  enddo
               enddo
            enddo
         endif
         write(iunit,888)

         else
         if (i_sided .ne. 1) then
            write(iunit,888)
            write(iunit,888) ' Diagonal Hessian parent energy '
            write(iunit, 901)
            write(iunit,888)
            do k = 1, nmode
               write(iunit,914) k, 
     $              sqrt(abs(scr(ifcm+k-1)))*5.14048D03,
     $              sqrt(abs(scr(ifcm+k-1))) *5.14048D03 / 8065.5
            enddo
         endif
c     
         if (.false. .and. (i_sided .ge. 3 .and. i_sided .ne. 5)) then
            write(6,*) ' off-diagonal elements hessian parent E'
c     
            k=0
            ioff = ifcm+nmode
            do irp = 1, nirrepf
               do i = 1, numirpf(irp)
                  do j = 1, i-1
                     write(iunit,916) irp, i, j,
     $                    sqrt(abs(scr(ioff+k)))*5.14048D03,
     $                    sqrt(abs(scr(ioff+k))) *5.14048D03 / 8065.5
                     k=k+1
                  enddo
               enddo
            enddo
         endif
         write(iunit,888)
         endif
c     
c     This was known information concerning the Parent State.
c     
c         call getrec(20,'JOBARC',refheff(icase),
c     $        ndim2*iintfp, scr(iref))
c     
         call getrec(20,'JOBARC','REFPAR_E',
     $        iintfp, e0)
         scr(iref) = e0
c     
c     subract reference energy from diagonal and convert to EV's
c     
         j = iref
         do i = 1, ndim
            scr(j) = scr(j) - e0
            j = j + ndim+1
         enddo
         call SSCAL(ndim2, 27.2113956324672d0, scr(iref), 1)
c     
         write(iunit,888) ' Reference Hamiltonian '
         call put_matrix(scr(iref), ndim, ndim, perline, iunit,.false.)
c     
         write(iunit,888)
         write(iunit,888) ' All information concerning method ',
     $        nameheff0
         write(iunit,888)
         write(6,*) ' @write_vibronic: done with ', nameheff0
c     
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
         write(iunit,888) ' Vibronic Coupling elements for case ',
     $        nameheff(icase)
         write(iunit,901)
         write(iunit,901)
 901     format(70('-'))
         write(iunit, 888)
c     
c     write common info for this calculation
c     
         write(iunit,811) i_sided
 811     format(t3,'Vibronic grid in ACESII : ', i4)
         write(iunit,812) nirrepf
 812     format(t3,'Total number of vibrational irreps', i4)
         write(iunit,813) (numirpf(irp), irp=1,nirrepf)
 813     format(t3,'Number of modes per symmetry', 20i4)
c     
         call getrec(1, 'JOBARC', 'FREQS_0 ', numvib*iintfp,
     $        freq)
         ioff = 1
         do irp = 1, nirrepf
            if (numirpf(irp) .ne. 0) then
               do i = 1, 4
                  fnameirp(i:i) = labelf(irp)(i:i)
               enddo
               write(iunit,814) irp, numirpf(irp), fnameirp,
     $              (freq(ioff+i), i=0,numirpf(irp)-1)
            endif
            ioff = ioff + numirpf(irp)
         enddo
 814     format(t3,'[',i2,']', t9,i4, t15, a4, t25, 20F8.2)
         write(iunit,888)
c     
c     process heff part
c     
         ndim2 = ndim*ndim
c     
         iref = 1
         ipnt = iref + 1
         igrd = ipnt + npoint
         ifcm = igrd + nmode
         ifreq = ifcm + nmode+nmode2
         itop = ifreq+nmode
c     
         do k = 1, nmode
            scr(ifreq+k-1) = fact / sqrt(abs(freq(k)))
         enddo
c     
c     first treat parent energy
c     
         call getrec(20,'JOBARC','REFPAR_E',
     $        iintfp, scr(iref))
         write(6,*) ' Ground state energy at reference point ',
     $        scr(iref)
c     
         call getrec(20,'JOBARC','PNTENERG',
     $        npoint*iintfp, scr(ipnt))
c     
         if (i_sided .eq. 1) then
            do k = 1, nmode
               scr(igrd+k-1) = 
     $              (scr(ipnt+k-1) - scr(iref)) / stpsiz
            enddo
         else
c     
c     two-sided gradients
c     
            ioff1 = ipnt-1
            ioff2 = ipnt + nmode -1
            do k = 1, nmode
               scr(igrd+k-1) = 
     $              (scr(ioff1+k) -
     $              scr(ioff2+k)) / (2.0d0*stpsiz)
            enddo
c     
c     diagonal elements of Hessian
c     
            do k = 1, nmode
               scr(ifcm+k-1) =
     $              (scr(ioff1+k) +
     $              scr(ioff2+k) - 2.0d0 * scr(iref)) / (stpsiz**2)
            enddo
c     
         endif
c     
         if (i_sided .ge. 3 .and. i_sided .ne. 5) then
c     
c     off-diagonal elements of Hessian (totally symmetric only)
c     
            ioff1 = ipnt - 1
            ioff2 = ioff1 + nmode
            ioff3 = ioff2 + nmode
            ioff4 = ioff3 + nmode2
            k=0
            do irp = 1, nirrepf
               do i = 1, numirpf(irp)
                  do j = i+1, numirpf(irp)
                     k=k+1
                     write(6,*)
     $                    scr(ioff3+k), scr(ioff4+k),
     $                    scr(ioff1+i), scr(ioff1+j),
     $                    scr(ioff2+i), scr(ioff2+j),
     $                    two*scr(iref)
c     
                     scr(ifcm+nmode+k-1) =
     $                    (scr(ioff3+k) + scr(ioff4+k)
     $                    - scr(ioff1+i) - scr(ioff1+j)
     $                    - scr(ioff2+i) - scr(ioff2+j)
     $                    + two*scr(iref) )
                     scr(ifcm+nmode+k-1) = scr(ifcm+nmode+k-1)
     $                    / (stpsiz**2)
                  enddo
               enddo
               ioff1 = ioff1 + numirpf(irp)
               ioff2 = ioff2 + numirpf(irp)
               ioff3 = ioff3 + numirpf(irp)*(numirpf(irp)-1) / 2
               ioff4 = ioff4 + numirpf(irp)*(numirpf(irp)-1) / 2
            enddo
c     
         endif
c     
         write(iunit,888)
         write(iunit,888) ' Gradients parent energy along normal modes'
         write(iunit, 901)
         write(iunit,888)
         do k = 1, nmode
            write(iunit,912) k, scr(igrd+k-1)
         enddo
 912     format(t3, 'Normal mode ', i4, t22, 'Gradient :', E16.10)
c
         if (reduced) then
            fact = 27.2114 * 8065.5 
         if (i_sided .ne. 1) then
            write(iunit,888)
            write(iunit,888) ' Diagonal Hessian parent energy '
            write(iunit, 901)
            write(iunit,888)
            do k = 1, nmode
               write(iunit,914) k, 
     $              abs(scr(ifcm+k-1))*fact, 
     $              abs(scr(ifcm+k-1)) * fact / 8065.5
            enddo
         endif
c     
         if (.false. .and. (i_sided .ge. 3 .and. i_sided .ne. 5)) then
            write(6,*) ' off-diagonal elements hessian parent E'
c     
            k=0
            ioff = ifcm+nmode
            do irp = 1, nirrepf
               do i = 1, numirpf(irp)
                  do j = 1, i-1
                     write(iunit,916) irp, i, j,
     $                    abs(scr(ioff+k))*fact,
     $                    abs(scr(ioff+k)) * fact / 8065.5
                     k=k+1
                  enddo
               enddo
            enddo
         endif
         write(iunit,888)

         else
         if (i_sided .ne. 1) then
            write(iunit,888)
            write(iunit,888) ' Diagonal Hessian parent energy '
            write(iunit, 901)
            write(iunit,888)
            do k = 1, nmode
               write(iunit,914) k, 
     $              sqrt(abs(scr(ifcm+k-1)))*5.14048D03,
     $              sqrt(abs(scr(ifcm+k-1))) *5.14048D03 / 8065.5
            enddo
         endif
 914     format(t3, 'Normal mode ', i4, t22, 
     $        ' Frequency : ', F8.2, ' cm-1', F12.6, ' eV')
c     
         if (.false. .and. (i_sided .ge. 3 .and. i_sided .ne. 5)) then
            write(6,*) ' off-diagonal elements hessian parent E'
c     
            k=0
            ioff = ifcm+nmode
            do irp = 1, nirrepf
               do i = 1, numirpf(irp)
                  do j = 1, i-1
                     write(iunit,916) irp, i, j,
     $                    sqrt(abs(scr(ioff+k)))*5.14048D03,
     $                    sqrt(abs(scr(ioff+k))) *5.14048D03 / 8065.5
                     k=k+1
                  enddo
               enddo
            enddo
         endif
 916     format(t3, '[', i1, ']', t8, i4, t15, i4, t23, F10.2,
     $        t35, 'cm-1', t40, F10.6, t52, 'eV')
         write(iunit,888)
         endif
c     
c     This was known information concerning the Parent State.
c     
         call getrec(20,'JOBARC',refheff(icase),
     $        ndim2*iintfp, scr(iref))
c     
         call getrec(20,'JOBARC','REFPAR_E',
     $        iintfp, e0)
         write(6,*) ' Do not subtract e0 '
         e0 = 0.0d0
c     
c     subract reference energy from diagonal and convert to EV's
c     
         j = iref
         do i = 1, ndim
            scr(j) = scr(j) - e0
            j = j + ndim+1
         enddo
         call SSCAL(ndim2, 27.2113956324672d0, scr(iref), 1)
c     
         write(iunit,888) ' Reference Hamiltonian '
         call put_matrix(scr(iref), ndim, ndim, perline, iunit,.false.)
c     
         write(iunit,888)
         write(iunit,888) ' All information concerning method ',
     $        nameheff(icase)
         write(iunit,888)
         write(6,*) ' @write_vibronic: done with ', nameheff(icase)
c     
 10   continue
c     
      close(unit=iunit, status='KEEP')
c     
      write(6,*) ' All done write_vibronic '
c     
      return
      end
