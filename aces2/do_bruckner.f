













































































































































































































      Subroutine do_brueckner

      Logical Converged
      Logical Rot_grad
      Integer Icycle 
C


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



      Write(6,"(a)"), "@-DO_POST_Brueckner Entering Brueckner block"
      Rot_grad = Iflags2(178) .Eq. 4
      Converged = .FALSE.
      Icycle = 1
      Do while (.NOT. Converged)
          If (Icycle .EQ. 1) Then
              Call Scf_eneg
          Else
              If (dkh .and. contract) Call Runit("xdkh")
              Call Runit("xvscf")
          Endif
          Call Change_orbital_space
          Call Prep4_post_scf
          Call Post_scf_eneg_4brueck(Rot_grad)
          Call A2getrec(20, 'JOBARC', 'BRUKTEST', 1, Itest)
          IF (Itest .EQ. 1) Converged = .TRUE.
          Icycle = Icycle + 1
      Enddo
c
      if (First_order_props) Then
          Call Runit("xlambda")     
          Call Runit("xdens")     
          Call Runit("xvprops")     
      Endif
C
C Excited state calculations with Brueckner orbitals. At convergence
C we have T2 amplitudes to form Hbar.
C
      If (Iflags(87) .GT. 2 .AND.
     &         Iflags2(174) .EQ. 1) Then
         Call Runit("xhbar")
         Call Runit("xfsip")
      Else If (Iflags(87) .GT. 2 .AND. (.NOT.
     &         Iflags2(174) .EQ. 1)) Then
         Call Runit("xhbar")
         Call Runit("xvee")
      Endif 
      
      Return
      End
