










      Subroutine Post_scf_eneg_4brueck(Rot_grad)
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
      Character*79 szEExtrap 
      Logical Rot_grad
      szEExtrap = 'xa2proc ener_extrp'
      Call a2getrec(20,'JOBARC','NDROPGEO',1,NDRGEO)
C
C The regular post_scr_energy routine has excited states logic. But
C we are now trying to develop Breuckner based excited state
C calculations. So, it become necessary first determine Brueckner
C and then proceed to exciated energies, 04/2015, Ajith Perera.
c
c This Saving of All MO JOBARC/JAINDX as JOBARC_AM/JAINDX_AM is only relevent
c for DROPMO gradients. It is done here even for energy only as a organizational
c convenince.
c
CSSS      If (Ndrgeo.ne.0)  then
CSSS         Call Runit('cp JOBARC JOBARC_AM')
CSSS         Call Runit('cp JAINDX JAINDX_AM')
CSSS      End if
C
c The Musial & Kucharski need to have HF2 file around which normaly
c get deleted after processing. I sure hope we can change this key-word!
c 
      
      If (Iflags2(152) .EQ. 1)  Then 
c
         Call Runit('xvtran')
         Call runit('cp HF2 HF2N')
         Call Runit('xintprc')
         Call runit('mv HF2N HF2')
         Call Runit('xccsdtq < COM')
c
      Else If (Dirmp2) Then
c
         Call Runit('xdirmp2')
c
C Finite-order Brueckner work for COPT, 11/2016 (Ajith Perera for
C Y. Park). 
C
      Else if (Iflags2(174) .EQ. 1) then
 
         Call Runit('xcopt')

      Else
C
         If (Iflags(87) .NE. 2 .AND. .NOT. Rot_grad) 
     &       Call Runit('xvcc')
         If (Iflags(87) .NE. 2 .AND. Rot_grad) 
     &       Call Runit('xpccd')
c
c Noniterative 5th-order triples and quadruples on top of CCSD
c
         If (Iflags(2) .EQ. 12 .OR.
     &      (Iflags(2) .GE. 26 .AND. 
     &       Iflags(2) .LE. 30)) Then
             Call Runit('xvcc5t')
             Call Runit('xvcc5q')
         Endif
C
         If (BlCCSDT) Then
           Call Runit('xlambda')
           Call Runit('xlcct')
         Endif
c
         If (Iflags(2) .EQ. 31) Call Runit('xvcc5t')
c
c Extrapolate energy (only for post-SCF methods), 01/2006, Ajith Perera
c
         If (bEextrap .OR. bCExtrap) Then
            Call Runit(szEExtrap)
         Endif
c
       Endif
c 
      Return
      End
