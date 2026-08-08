










      Subroutine Post_scf_higher_order_prps(DO_NMR_SHIFTS,
     &                                      DO_NMR_SPNSPN,
     &                                      DO_NLO_PROPS)
c
       Logical DO_NMR_SHIFTS, DO_NMR_SPNSPN, DO_NLO_PROPS 
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
      IF (DO_NMR_SHIFTS) Then
c
         Call Runit('xlambda')
         Call Runit('xdens')
         Call Runit('xvdint')
         Call Runit('xcphf')
         Call Runit('xnmr')
c
        If (iflags(90).ne.0) then
c
c o sequential treatment of B-field components requires
c   reevaluation of GIAO integrals
c
            call runit('rm -f I*X')
            call runit('xvdint')
            call runit('xnmr')
            call runit('rm -f I*Y')
            call runit('xvdint')
            call runit('xnmr')
            call runit('rm -f I*Z')
         End if
         Call Runit('xjoda')
c
      Else if (DO_NMR_SPNSPN .OR. DO_NLO_PROPS) Then
c
         Call Runit('xlambda')
         Call Runit('xvprops')
         Call Runit('xvcceh')
c
      Endif
c
      Return
      End
