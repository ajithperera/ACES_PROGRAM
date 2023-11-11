











      subroutine runoombpt
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


      integer ndrgeo,itest
      logical bDoneOO

      bDoneOO = .false.
      ndrgeo = 0
      itest = 0

      call a2getrec(20,'JOBARC','IFLAGS',100,iflags)
      do while (.not. bDoneOO)

        call runit('xoombpt2')
        call a2getrec(20,'JOBARC','IFLAGS',100,iflags)
        call a2getrec(20,'JOBARC','OOMBPTST',1,itest)

C If itest == 2 then error
        if (itest .eq. 2) stop
C If itest == 0 or 1 run xvscf and precc
        call a2putrec(20,'JOBARC','NOTRANS ',1,0)
        If (dkh .and. contract) Call Runit("xdkh")
        call runit('xvscf')
        call precc(ndrgeo)

        if (itest .eq. 1) bDoneOO = .true.

      end do

      iflags(2) = 1
      call a2putrec(20,'JOBARC','IFLAGS',100,iflags)

      call runvcc

      return
      end
