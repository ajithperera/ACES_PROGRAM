










      Subroutine Scf_higher_order_prps(DO_TDHF, DO_NMR_SHIFTS, 
     &                                 DO_NMR_SPNSPN, DO_NLO_PROPS)
c  
      Logical DO_TDHF, DO_NMR_SHIFTS, DO_NMR_SPNSPN,
     &        DO_NLO_PROPS
c


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


c
      If (DO_TDHF) Then
         Call Runit('xvprops')
         Call Runit('xvtran')
         Call Runit('xtdhf')
      Else if (DO_NMR_SHIFTS .OR. DO_NMR_SPNSPN) Then
         Call Runit('xvtran')
         Call Runit('xintprc')
         Call Runit('xvdint')
         Call Runit('xcphf')
         Call Runit('xjoda')
      Else if (DO_NLO_PROPS) Then
         Call Runit('xvtran')
         Call Runit('xintprc')
         Call Runit('xvdint')
         Call Runit('xcphf')
      Endif
c
      Return
      End
