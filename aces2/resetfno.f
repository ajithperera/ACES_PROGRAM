











      subroutine resetfno(eval,nbas,file)
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


C     Common blocks
      integer iintln,iifltln,iintfp 
      common/machsp/iintln,iifltln,iintfp 
C     Input Variables
      integer nbas
      character*2 file
C     Pre-allocated Local Variables
      double precision eval(nbas)      
C     Local variables
      integer nirrep,nocco(2),vrtdrop(8,2),vrt(8,2),uhf,virtd,virtk,
     &   virt,irrep,ii,spin,ieval,dropocc
      character*8 recname(2)
      data recname /'SCFEVLA0','SCFEVLB0'/
C - - - - - - -- -- - - - - - - - - - - - - - - - - - - - - - - -
      uhf = 0
      if (iflags(11).ne.0) uhf = 1
      call aces_ja_init
      call getrec(20,'JOBARC','NOCCORB',2,nocco)
      call getrec(20,'JOBARC','COMPNIRR',1,nirrep)
      call getrec(20,'JOBARC','SYMPOPVA',nirrep,vrt(1,1))
      call getrec(20,'JOBARC','SYMPOPVB',nirrep,vrt(1,2))
      call getrec(20,'JOBARC','FNODROP ',16,vrtdrop)
      call getrec(20,'JOBARC','FNOFREEZ',1,dropocc)
      do spin=1,uhf+1
         call getrec(20,'JOBARC',recname(spin),nbas*iintfp,eval)
         if (file.eq.'DM') then
            ieval=nocco(spin)+dropocc
         else
            ieval=nocco(spin)
         endif
         do irrep=1,nirrep
            virt=vrt(irrep,spin)
            virtd=vrtdrop(irrep,spin)
            if (file.eq.'DM') then
               virtk=virt
            else
               virtk=virt-virtd
            endif
            ieval=ieval+virtk
            do ii=1,virtd
               eval(ieval+ii)=eval(ieval+ii)-999.0d0
            end do
            ieval=ieval+virtd
         end do
         call putrec(20,'JOBARC',recname(spin),nbas*iintfp,eval)
      end do
      call aces_ja_fin
      end
