











      subroutine runvcc



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



c KJW 9/29/99
c support for Stan Kucharski's' ccsdtq program
      if (iflags2(152).eq.1) then
         call runit('xrdycom')
         call runit('xccsdtq < COM')
         return
      end if

1     call runit('xvcc')

c JDW 3/31/94
c noniterative 5th-order triples and quadruples
      if (iflags(2).eq.12.or.(iflags(2).ge.26.and.iflags(2).le.30)) then
         call runit('xvcc5t')
         call runit('xvcc5q')
      end if

c JDW 3/17/95
c CC5SD[T]
      if (iflags(2).eq.31) call runit('xvcc5t')

c Extrapolate energy (only for correlated methods), 01/2006, Ajith Perera
      if (bEExtrap.or.bCExtrap) then
         call runit('xa2proc extrap energy')
      end if

      if (iflags(22).eq.1) then
         call a2getrec(20,'JOBARC','BRUKTEST',1,itest)
         if(itest.eq.1)return
         if (nddo_guess) call runit('xnddo')
         If (dkh .and. contract) Call Runit("xdkh")
         call runit('xvscf')
         call runit('xvtran')
         call runit('xintprc')
         goto 1
      end if

      end

