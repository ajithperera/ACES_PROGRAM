










      Subroutine Freq_numrcl(Calc_level, Analytical_grads)
c
      Implicit None
c


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





































































































































































































c
      Character*4 Calc_level
      Character*79 szGExtrap  
      Logical Analytical_grads
      Integer Istat, Is, Ius, Irank
c
      szGExtrap = 'xa2proc grad_extrap'
c
c This shell script (with minor mods) can perform a restartable SCF vib freq
c rmjunk; restore; xjoda || exit 1
c analytical_gradient=$(some_test)
c while true
c do xvmol; xvmol2ja; xvscf 
c    test analytical_gradient && xvdint
c    xa2proc clrdirty
c    backup; restore; xjoda || exit 1
c    test $(xa2proc jareq FNDFDONE i 1 | tail -1) -ne 0 && break
c done
c clean; rmbackup
c exit
c
c The goal is to compartmentalize the computation of the "targets" (energy, polarizability..)
c that we do the finite differencing. Need to be very carefull about the overlaoded variable
c Analytical_grads. 

      If (Calc_level .Eq. "SCF ") Then
         iRank=0
         istat=0
         Do while (istat.eq.0)
c
c The energy block, this is the entry point for any future one particle energy methods.
c
            call scf_eneg
c
            If (Analytical_grads) then
               Call Scf_anlytcl_grad(.TRUE.)
            End if
c
            call a2getrec(1,'JOBARC','LASTGEOM',1,istat)
c
            if (istat.eq.0) then
               call c_gtod(is,ius)
               print '(a,f10.1,a)',
     &               'ACES2: Total elapsed time is ',
     &                is+1.d-6*ius-dTimeStart,' seconds'
               call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
               call runit(xjoda)
            else

c
               if (irank.eq.0) then
                  call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
                  If (.NOT. Manual_FD) call runit(xjoda)
               end if
            Endif 
c
            call rmfiles
         End do
c
c Else for the Calc_level
c
      Else
c         
         iRank=0
         istat=0
         do while (istat.eq.0)
c
            Call Scf_eneg
c
c Post-SCF energy generating memeber executables are "xdirmp2" (AO-direct MBP2(22)
c "xvcc" (ACES II CC program) and "xvee" is ACES2 EOM-CC program. This is the 
c entry point for all the future energy generating methods.
c
            If (dirmp2) Then
               Call runit('xdirmp2')
            Else
               Call Change_orbital_space
               Call Prep4_post_scf   
               Call Post_scf_eneg
c
               If (iflags(87).gt.1) Then
                  Call Runit('xlambda')
                  Call Runit('xvee')
               End if
            Endif
c
c Any of these energy generating methods has analytical gradients? 
c
            If (Analytical_grads) then
c
               Call Post_scf_anlytcl_grad(.TRUE.)
c
c If the gradients need to modified externally (ie. extrapolation)
c do it here (note that no extrapolation for HF-SCF calculations),
c 01/2006, Ajith Perera. When POST-HF exact frequencies are availble
c the extrpolation has to be done differently. 
c
               If (bGExtrap .or. bCExtrap) then
                  call runit(szGExtrap)
               endif
            Endif
c
            call a2getrec(1,'JOBARC','LASTGEOM',1,istat)
c            
            if (istat.eq.0) then
CSSS               call c_gtod(is,ius)
CSSS               print '(a,f10.1,a)',
CSSS     &               'ACES2: Total elapsed time is ',
CSSS     &                is+1.d-6*ius-dTimeStart,' seconds'
               call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
               call runit(xjoda)
            else
c
c
               if (irank.eq.0) then
                  call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
                  If (.NOT. Manual_FD) call runit(xjoda)
               end if
c
            end if
            call rmfiles
c
         End do
c
      Endif

      Return
      End
