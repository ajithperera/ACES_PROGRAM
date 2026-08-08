










      Subroutine Scf_anlytcl_grad(Do_derint)



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



      integer icore, is, ius
      logical Do_derint


CSSS      IF (Do_derint) call Runit(der_integral_package)
C
      If (Hf_Dft .and. Do_derint) Then
C
         Call Runit('xvtran')
         Call Runit('xintprc')
         Call Runit(der_integral_package)
         Call Runit('xvksdint')
C       
      Else If (Do_derint .and. ks_scf) Then 
C
         Call Runit(der_integral_package)
         Call Runit('xvksdint')
C
      Else If (Do_derint) Then
C            
         Call Runit(der_integral_package)
C
      End If
     
      call c_gtod(is,ius)
      print '(a,f10.1,a)', 'ACES2: Total elapsed time is ',
     &                      is+1.d-6*ius-dTimeStart,' seconds'
         
         
      Return
      End


