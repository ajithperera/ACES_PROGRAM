








































































































































































































      subroutine precc(ndrgeo)
      implicit none

      integer ndrgeo, ndropa
      integer itest



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



c   o vtran/intprc were already run if HFSTAB=(ON|FOLLOW)
c     we have to re-run them if (HFSTAB=FOLLOW && 'SCFKICK'!=0) || DROPMO!=0
      if (iflags(74).eq.2) then
         call a2getrec(-1,'JOBARC','SCFKICK',1,itest)
      else
         itest=0
      end if

c   o use NUMDROPA since ndrgeo is 0 for single-point calcs
      call a2getrec(-1,'JOBARC','NUMDROPA',1,ndropa)

      if (iflags(74).eq.0.or.ndropa.ne.0.or.
     &    (iflags(74).eq.2.and.itest.ne.0)  ) then
         if (ndrgeo.ne.0)  then
            call runit('cp JOBARC JOBARC_AM')
            call runit('cp JAINDX JAINDX_AM')
         end if
         call runit('xvtran')
         if (iflags2(152).eq.1) then
c         o KJW 9/29/99 for Stan Kucharski's ccsdtq program
            call runit('cp HF2 HF2N')
            call runit('xintprc')
            call runit('mv HF2N HF2')
         else
            call runit('xintprc')
         end if
      end if

      return
      end

