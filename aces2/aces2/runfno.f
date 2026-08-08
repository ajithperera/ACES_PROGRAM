











      subroutine runfno
      integer ndropa,icalc,iabcd,igamma_abcd
      logical IntsExist


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





































































































































































































      
      ndrgeo1=0
      call aces_ja_init
      call getrec(20,'JOBARC','IFLAGS  ',100,iflags)
      icalc=iflags(2)
      iflags(2)=1
      iabcd=iflags(93)
      iflags(93)=0
      ivtran=iflags(83)
      iflags(83)=0
      igamma_abcd=iflags(100)
      iflags(100)=0
      call putrec(20,'JOBARC','IFLAGS  ',100,iflags)
      call getrec(20,'JOBARC','NUMDROPA',1,ndropa)
      call putrec(20,'JOBARC','FNOFREEZ',1,ndropa)
      call putrec(20,'JOBARC','NUMDROPA',1,0)
      if (iflags(11).ne.0) then
         call putrec(20,'JOBARC','NUMDROPB',1,0)
      end if
      if ((geom_opt.or.vib_specs).and.analytical_gradient) then
         print *,'@RUNFNO FNO gradients not yet supported.'
         call aces_exit(1)
         call putrec(20,'JOBARC','NDROPGEO',1,0)
      endif
      call aces_ja_fin
      call runit('xvtran')
      call runit('xintprc')
      call runit('xvcc')
      call runit('xfno')

      call aces_ja_init
      call getrec(20,'JOBARC','IFLAGS  ',100,iflags)
      call getrec(20,'JOBARC','IFLAGS2 ',500,iflags2)
      iflags(2)=icalc
      iflags(93)=iabcd
      iflags(83)=ivtran
      iflags(100)=igamma_abcd
      call putrec(20,'JOBARC','IFLAGS  ',100,iflags)
      call aces_ja_fin
      inquire(file='IIII',exist=IntsExist)
      if (.not.IntsExist) then
         call runit(integral_package)
         call runit('xvmol2ja')         
      endif
      return
      end

