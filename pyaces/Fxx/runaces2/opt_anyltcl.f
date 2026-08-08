










      Subroutine Opt_anlytcl(Calc_level)



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

      integer icore, is, ius
      logical Do_derint
      Character*4 Calc_level
      character*79 szGExtrap 


      ndrgeo=0
      szGExtrap = 'xa2proc grad_extrp'
      call a2getrec(20,'JOBARC','NDROPGEO',1,NDRGEO)

c   o generate an initial approximate Hessian. 12/04, Ajith Perera
      if (mopac_guess) call runit('xmopac')

      icycle=0
c
 1000 continue
c
      icycle=icycle+1

      If (Calc_level .Eq. "SCF " ) Then
c
c Let's not run the integral derivetive during this call to scf_anlytcl_grad
c to avoid duplicate runs if and when analytic Hessian is calculated, but set
c do_derint to TRUE
c
         Call Scf_eneg
         Call Scf_anlytcl_grad (.FALSE.)
         Do_derint = .TRUE.
c 
c Since analytic Hessian is available for HF-SCF, we can use it guide the optimization
c by recomputing the Hessian by interval (specified by the user).
c
         If (iflags(55).gt.0 .and. .not. (ks_scf .or. hf_dft)) then
c
            ijunk=mod(icycle-1,iflags(55))
c
            if (ijunk.eq.0) then
               Do_derint=.FALSE.
               itmp=iflags(54)
               itmp2=iflags(3)
               iflags(3)=2
               iflags(54)=1
               call a2putrec(20,'JOBARC','IFLAGS  ',100,iflags)
               call runit('xvtran')
               call runit('xintprc')
c
c AP - Whenever it is done, McKinley from MOLCAS program should go here
c
               call runit('xvdint')
               call runit('xcphf')
C--DEBUG ONLY---
CSSS               call runit('xprep')
C--DEBUG ONLY---

               iflags(3)=itmp2
               iflags(54)=itmp
               call a2putrec(20,'JOBARC','IFLAGS  ',100,iflags)
            end if
c
c When EVAL_HESS=0, generate a Hessian only at the start. Equivalent
c to using FCMINT file from VIB=EXACT calculation.
c
         Else if (iflags(55) .eq. 0 .and. .not. ks_scf) Then
            if (icycle-1 .eq. 0) then
               Do_derint=.FALSE.
               itmp=iflags(54)
               itmp2=iflags(3)
               iflags(3)=2
               iflags(54)=1
               call a2putrec(20,'JOBARC','IFLAGS  ',100,iflags)
               call runit('xvtran')
               call runit('xintprc')
c
c AP - Whenever it is done, McKinley from MOLCAS program should go here
c
               call runit('xvdint')
               call runit('xcphf')
C--DEBUG ONLY---
CSSS               call runit('xprep')
C--DEBUG ONLY---

               iflags(3)=itmp2
               iflags(54)=itmp
               call a2putrec(20,'JOBARC','IFLAGS  ',100,iflags)
            end if
         Endif

         if (Do_derint .and. ks_scf) then
             Call runit(der_integral_package)
             Call Runit('xvksdint')
         else if (Do_derint .and. hf_dft) then 
             Call runit('xvtran')
             Call runit('xintprc')
             Call runit(der_integral_package)
             Call Runit('xvksdint')
         else if (Do_derint) then
             call runit(der_integral_package)
         endif
C
         call c_gtod(is,ius)
         print '(a,f10.1,a)', 'ACES2: Total elapsed time is ',
     &                        is+1.d-6*ius-dTimeStart,' seconds'
c
         call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
         call runit('xjoda')
         call a2getrec(1,'JOBARC','JODADONE',1,istat)
         if (istat.ne.0) return
         call rmfiles
         goto 1000
c
c The Else block for the Post SCF optimizations
c
      Else if (Calc_level .EQ. "PSCF") Then
c
         Call Scf_eneg
         If (Iflags(87) .GT. 2) Call Runit
     &                                       ("xvprops")
         Call Change_orbital_space
         Call Prep4_post_scf
CSSS         Call Post_scf_eneg
CSSS         Call Post_scf_anlytcl_grad(.TRUE.)
         Do_derint = .FALSE. 
c
c At the moment we do not have analytic Hessians, even when we have analytic
c Hessian, I don't think it is wise to compute analytic Hessian to guide
c optimization, instead we can use the HF-SCF Hessian.
c
         if (iflags(55) .gt. 0) then
c
            ijunk=mod(icycle-1,iflags(55))
c
            if (ijunk.eq.0) then
               itmp0=iflags(2)
               itmp=iflags(54)
               itmp2=iflags(3)
               iflags(2)=0
               iflags(3)=2
               iflags(54)=1
               call a2putrec(20,'JOBARC','IFLAGS  ',100,iflags)
c
c AP - Whenever it is done, McKinley from MOLCAS program should go here
c
               call runit('xvdint')
               call runit('xcphf')
               iflags(2)=itmp0
               iflags(3)=itmp2
               iflags(54)=itmp
               call a2putrec(20,'JOBARC','IFLAGS  ',100,iflags)
            end if
c
c When EVAL_HESS=0, generate a Hessian only at the start. Equivalent
c to using FCMINT file from VIB=EXACT calculation.
c
         Else if (iflags(55) .eq. 0) Then 
             If (Icycle-1 .eq. 0) then
                itmp0=iflags(2)
                itmp=iflags(54)
                itmp2=iflags(3)
                iflags(2)=0
                iflags(3)=2
                iflags(54)=1
                call a2putrec(20,'JOBARC','IFLAGS  ',100,iflags)
c
c AP - Whenever it is done, McKinley from MOLCAS program should go here
c
                call runit('xvdint')
                call runit('xcphf')
                iflags(2)=itmp0
                iflags(3)=itmp2
                iflags(54)=itmp
                call a2putrec(20,'JOBARC','IFLAGS  ',100,iflags)
             end if
         end if
c
         Call Post_scf_eneg
         Call Post_scf_anlytcl_grad(.TRUE.)
c
c We have no choice but run integral derivative code twice if we ues analytic SCF
c Hessian to guide the optimization.
c
         IF (Do_derint) Call Runit(der_integral_package)
c
c If the gradients need to modified externally (ie. extrapolation) 
c do it here (note that no extrapolation for HF-SCF calculations),
c 01/2006, Ajith Perera. 
c
c         If (bGExtrap .or. bCExtrap) then
c            call runit(szGExtrap) 
c         endif
c
CSSS         call c_gtod(is,ius)
CSSS         print '(a,f10.1,a)', 'ACES2: Total elapsed time is ',
CSSS     &                        is+1.d-6*ius-dTimeStart,' seconds'
c
         call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
         call runit('xjoda')
         call a2getrec(1,'JOBARC','JODADONE',1,istat)
         if (istat.ne.0) return
         call rmfiles
         goto 1000
c
c Endif for the Calc_level
c
      End if

      Return
      End

