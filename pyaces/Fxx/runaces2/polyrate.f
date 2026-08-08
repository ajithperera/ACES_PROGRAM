











C JDW 5/12/95 Force calculation in POLYRATE jobs
C     5/25/95 Only let this happen if we have cartesian coordinates.
C             Note that we can in principle do a force with dropped
C             core, we cannot do hessian, hence dropped core addition
C             here is currently of minor value.

      subroutine polyrate



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




      ndrgeo=0
      call a2getrec(20,'JOBARC','NDROPGEO',1,ndrgeo)

      if (.not.direct) call runit(integral_package)
      call runit('xvmol2ja')
      if (nddo_guess) call runit('xnddo')
      If (dkh .and. contract) Call Runit("xdkh")
      call runit('xvscf')
      if (iflags(2).gt.0) then
         if (ndrgeo.ne.0) then
            call runit('cp JOBARC JOBARC_AM')
            call runit('cp JAINDX JAINDX_AM')
         end if
         call runit('xvtran')
         call runit('xintprc')
         call runit('xvcc')
         call runit('xlambda')
         if (ndrgeo.ne.0) call runchg
         call runit('xdens')
         call runit('xanti')
         call runit('xbcktrn')
      end if
      call runit('xvdint')

      return
      end

