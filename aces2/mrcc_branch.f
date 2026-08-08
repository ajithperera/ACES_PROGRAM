











      subroutine mrcc_branch



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






































































































































































































      integer istat, ishell
      character*79 szHyper, szMolden

      write(6,*) ' Entered mrcc_branch '
      write(6,*) ' geom_opt : ', geom_opt
      write(6,*) ' vib_specs: ', vib_specs
      write(6,*)
c
      if (geom_opt .and. vib_specs) then
         write(6,*) ' @MRCC_new: optimize geometries + frequencies'
         if (analytical_gradient) then
            istat=0
            do while (istat.eq.0)
               call mrcc_grad
               call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
               call runit('xjoda')
               call a2getrec(1,'JOBARC','JODADONE',1,istat)
            end do
         else
            call opt_numrcl_4mrcc
         end if
         Call Prep4_post_opt_freq 
         call mrcc_vib
         return
      end if


      if (geom_opt) then
         if (analytical_gradient) then
            istat=0
            do while (istat.eq.0)
               call mrcc_grad
               call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
               call runit('xjoda')
               call a2getrec(1,'JOBARC','JODADONE',1,istat)
            end do
         else
            call opt_numrcl_4mrcc
         end if
         return
      end if

      if (vib_specs) then
         call mrcc_vib
         return
      end if
c
      if (iflags2(3) .eq. 1) then
cmn
c Resonance Raman calculation
cmn
         write(6,*) ' @mrcc_branch: enter mrcc_grad'
          call mrcc_grad
C
c 
       else if (iflags2(156).eq.1 .or.
     $         iflags2(160).eq.1 .or.
     $         iflags2(165).eq.1 .or.
     $          iflags2(158).eq.1 ) then
c
         write(6,*) ' @mrcc_branch: enter mrcc_heff'
C
          call mrcc_heff

       else
         write(6,*) ' @mrcc_branch: enter energy_4mrcc', mrcc
          call energy_4mrcc
          szMolden  = 'xa2proc molden'
          call runit(szMolden)
      end if
c
      return
      end

