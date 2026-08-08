











      subroutine opt_numrcl(Calc_level)



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


c parallel_aces.com : begin

c This common block contains the MPI statistics for each MPI process. The values
c are initialized in the acescore library.

      external aces_bd_parallel_aces




      integer                nprocs, irank, icpuname

      character*(256) szcpuname

      common /parallel_aces/ nprocs, irank, icpuname,
     &                       szcpuname
      save   /parallel_aces/

c parallel_aces.com : end

      Character*4 Calc_level
      Logical Converged

c   o generate an initial approximate Hessian. 12/04, Ajith Perera
      if (mopac_guess) call runit('xmopac')

      istat = 0
      do while (istat.eq.0)
c
         If (Calc_level .Eq. "SCF ") Then
            Call Scf_eneg
         Else
c
c At the moment I do not know whether HFSTABILITY=FOLLOW should be allowed.
c Questions are do we do it only in the first cycle, or every cycle or 
c intermittently.

            If (h_IFLAGS_brueckner .EQ. 1) Then
               Do while (.NOT. Converged)
                  Call Scf_eneg
                  Call Change_orbital_space
                  Call Prep4_post_scf
                  Call Post_scf_eneg 
                  Call A2getrec(20, 'JOBARC', 'BRUKTEST', 1, Itest)
                  IF (Itest .EQ. 1) Converged = .TRUE.
               Enddo
            Else
               Call Scf_eneg
               Call Change_orbital_space
               Call Prep4_post_scf                
               Call Post_scf_eneg
            Endif
         Endif
c
CSSS         call c_gtod(is,ius)
CSSS         print '(a,f10.1,a)', 'ACES2: Total elapsed time is ',
CSSS     &                        is+1.d-6*ius-dTimeStart,' seconds'
         call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
         call runit(xjoda)
         call a2getrec(1,'JOBARC','JODADONE',1,istat)
      end do

      return
      end

