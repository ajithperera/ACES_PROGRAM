















































































































































































































      subroutine init
      implicit none



      double precision     dTimeStart
      common /aces2_times/ dTimeStart
      save   /aces2_times/

      integer iflags(100),iflags2(500)
      common /flags/ iflags
      common /flags2/ iflags2

      logical
     &        mrcc,
     &        geom_opt, analytical_gradient,numerical_gradient,
     &        raman, vib_specs,Cmpt_props_numrcl,
     &        fno, hf_scf, ks_scf, sewswitch, plain_scf, direct,
     &        nddo_guess, dirmp2, mopac_guess, bexport, 
     &        Single_point, Single_point_energy,	
     &        Single_point_gradient, First_order_props,
     &        Higher_order_props, NMR_SHIFTS, NMR_SPNSPN,
     &        NLO_PROPS, TDHF, bHyper, bMolden,blCCSDT,
     &        BEextrap, BGextrap, BCextrap, Analytical_hessian,
     &        hf_dft,sl_oep,Manual_FD,dkh,contract 
      common /sv_bool/
     &        mrcc,
     &        geom_opt, analytical_gradient,numerical_gradient,
     &        raman, vib_specs,Cmpt_props_numrcl,
     &        fno, hf_scf, ks_scf, sewswitch, plain_scf, direct,
     &        nddo_guess, dirmp2, mopac_guess, bexport,
     &        Single_point, Single_point_energy,	
     &        Single_point_gradient, First_order_props,
     &        Higher_order_props, NMR_SHIFTS, NMR_SPNSPN,
     &        NLO_PROPS, TDHF, bHyper, bMolden,blCCSDT,
     &        BEextrap, BGextrap, BCextrap, Analytical_hessian,
     &        hf_dft,sl_oep,Manual_FD,dkh,contract
      save   /sv_bool/

      character
     &          xjoda*79,
     &          integral_package*79,
     &          der_integral_package*79
      common /sv_char/
     &          xjoda,
     &          integral_package,
     &          der_integral_package
      save   /sv_char/


      integer linblnk, iTmp
      character arg*20, alaska_inp*11
      integer i, f_iargc, num_args,length
      logical file_exists, gen_quartic



c ---------------------------------------------------------------------------
c Analyze the IFLAGS arrays first -------------------------------------------
c ---------------------------------------------------------------------------
c                  = 
c .TRUE. if INIT_HESSIAN=MOPAC
      mopac_guess = (iflags2(8).eq.2)

c .TRUE. if DIRECT=ON and INTEGRALS!=GAMESS
      direct=((iflags2(154).eq.1).and.
     &        (iflags(56).ne.5))

c .TRUE. if CALC=SCF and EXCITE=NONE
      plain_scf=((iflags(2).eq.0).and.(iflags(87).lt.3))

c .TRUE. if GUESS=NDDO
      nddo_guess = (iflags(45).eq.2)

c .TRUE. if SCF_TYPE=HF
      hf_scf=(iflags2(153).eq.0)
      ks_scf=(iflags2(153).eq.1)
      hf_dft=(iflags2(153).eq.2)
      sl_oep=(iflags2(153).ge.3)
      dkh   =(iflags2(167).gt.0)
      contract=(iflags2(168).eq.0)

c .TRUE. if DIRECT=ON, INTEGRALS=GAMESS, and CALC=MBPT(2)
      dirmp2=((iflags2(154).eq.1).and.
     &        (iflags(56).eq.5).and.
     &        (iflags(2).eq.1))

      fno=(iflags2(144).gt.0)

      geom_opt=(iflags2(5).ne.0)

      raman=(iflags2(151).eq.1)
      vib_specs=(iflags(54).ne.0)

c Analytical_gradient=.TRUE. if GRAD_CALC=ANALYTICAL
c Numerical_gradient =.TRUE. if GRAD_CALC=NUMERICAL

      analytical_gradient=(iflags2(138).eq.1)
      numerical_gradient =(iflags2(138).eq.2) 
c
c .TRUE. if the DERIVATIVE_LEVEL=SECOND
      Analytical_hessian = (iflags(3) .eq. 2) .and.
     &                     .NOT. hf_dft
c
c In general, this is wrong since ACES3 is not just mrcc code.
      mrcc=(iflags2(132).eq.3)
c
cmn include procedure to do multiple points on a surface
cmn  for the purpose of calculating quartic force constants.
cmn  In MRCC branch a number of such options are available. Here only one.
c
      gen_quartic = iflags2(165).eq.1
     $     .and. .not. mrcc
c
c Single Point calculations: Energy;gradients;hessians
c There is no direct keyword that we can acces to choose single point.
c The single point is true if not a geometry optimization, finite difference 
c frequency or finite difference property (energy gradient included) calculation. 
c
      Single_point          = .NOT. (Geom_opt .OR. Vib_specs .OR. 
     &                               Numerical_gradient .OR.
     $                                  gen_quartic) 
      Single_point_energy   = (Single_point .AND. (.NOT. 
     &                         (Analytical_gradient .OR. 
     &                          Numerical_gradient)))
c                
c First_order_props; .TRUE. if PROPS = FIRST_ORDER (1)
c 
      First_order_props = (Iflags(18) .EQ. 1)
c
c Higher_order_props; .TRUE. if PROPS = NMR, EOM_NLO, J_FC,....(> 1)
c
      Higher_order_props = (Iflags(18) .GT. 1)
      print *, 'Higher_order_props = ', Higher_order_props
c
      Single_point_gradient = (Analytical_gradient .OR.
     &                         Numerical_gradient) .AND. 
     &                                             .NOT. 
     &                        (First_order_props   .OR. 
     &                         Higher_order_props  .OR.
     &                         Geom_opt            .OR. 
     &                         Vib_specs)
c
c NMR_SHIFTS; .TRUE. if PROPS = NMR
c
      NMR_SHIFTS = ((Iflags(18) .EQ. 3) .OR.
     &             (Iflags(18) .EQ.  4) .OR.
     &             (Iflags(18) .EQ.  5) .OR.
     &             (Iflags(18) .EQ.  6) .OR.
     &             (Iflags(18) .EQ. 12))

c
c NMR_SPNSPN; .TRUE. if PROPS = J_FC,J_SD,J_PSO,J_DS), JSC_ALL
c
      NMR_SPNSPN = ((Iflags(18) .EQ. 8) .OR.
     &             (Iflags(18) .EQ.  9) .OR.
     &             (Iflags(18) .EQ. 10) .OR.
     &             (Iflags(18) .EQ. 13))
c
c NLO_PROPS; .TRUE. if PROPS = EOM_NLO
c
      NLO_PROPS = (Iflags(18) .EQ. 11) .OR. 
     &            (Higher_order_props)
c
c TDHF; .TRUE. if PROPS = TDHF
c
      TDHF = (Iflags2(103) .NE. 0)
      print *, 'TDHF = ', TDHF
      Cmpt_props_numrcl = (Analytical_gradient .AND. Raman 
     &                                        .AND. Vib_specs)
     &                .OR. (.NOT. Analytical_gradient   .AND. 
     &                           (Single_point_gradient .OR.
     &                            Higher_order_props))
c
c bHyper;bMolden; TRUE if EXTERNAL = HYPERCHEM, MOLDEN
c
      bHyper  = (iflags2(6) .eq. 1)
      bMolden = (iflags2(6) .eq. 2)
c
c blCCSDT; TRUE if CALCLEVEL = ACCSD(T)
c
      blCCSDT = (iflags(2).eq.42)
c
c bGExtrap;bEExtrap;bCExtrap TRUE if EXTRAPOLATE = GRADIENT, ENERGY, COMBO
c
      bGExtrap = (iflags2(9) .eq .1)
      bEExtrap = (iflags2(9) .eq .2)
      bCExtrap = (iflags2(9) .eq .3)
c ---------------------------------------------------------------------------
c Analyze the command-line arguments next -----------------------------------
c ---------------------------------------------------------------------------

c set the defaults and override them with the CL.
      sewswitch=.FALSE.
      integral_package='xvmol'
      if (hf_scf) then
         if ((iflags2(154).eq.1).and.
     &       (iflags(56).eq.5).and.
     &       (iflags(2).eq.0)) then
            der_integral_package='xscfgrd'
         else
            der_integral_package='xvdint'
         end if
      else if (ks_scf .or. hf_dft) then
         der_integral_package='xvdint'
      end if

      num_args = f_iargc()

c BEGIN GETARG LOOP
      do i = 1, num_args
         call f_getarg(i,arg)

c SEWARD
         if (arg(1:1).eq.'s') then
            sewswitch=.TRUE.
            integral_package= 'molcas run seward MOLCAS.INP'
         end if

c ALASKA
         if (arg(1:1).eq.'a') then
            inquire(file='ALASKA.INP',exist=file_exists)
            if (file_exists) then
               alaska_inp='ALASKA.INP'
            else
               alaska_inp='MOLCAS.INP'
            endif
            der_integral_package= 'molcas run alaska '//alaska_inp
         end if

c END GETARG LOOP
      end do

c verify consistency
      if ((.not.sewswitch).and.(iflags(56).eq.4)) then
         write(*,*) '@ACES2: VMOL cannot generate MOLCAS files.'
         write(*,*) '        Try running `xaces2 s` instead.'
         call c_exit(1)
      end if
      if (dirmp2.and.analytical_gradient.and.(geom_opt.or.vib_specs))
     &   then
         write(*,*) '@ACES2: There are no gradients for direct MBPT(2).'
         call c_exit(1)
      end if
      if (dirmp2.and.iflags(74).ne.0) then
         write(*,*) '@ACES2: HFSTAB does not work for direct MBPT(2).'
         call c_exit(1)
      end if
      if ( ((.not.plain_scf).and.(.not.hf_scf)).and.
     &     (geom_opt.or.vib_specs.or.raman)             ) then
         write(*,*)
         write(*,*) '@ACES2: WARNING! post-scf energy gradients ',
     &              'using KS-SCF orbitals are not fully tested.'
 1       write(*,*)
      end if



      return
      end

