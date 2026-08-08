













































































































































































































      Subroutine Prep4_post_opt_freq
      
      Implicit none
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


c     Maximum string length of terminal lines
      INTEGER LINELEN
      PARAMETER (LINELEN=80)

      Integer Iflag_store1(100), Length, I_havegeom, iOPTARC,
     &        Iflag_store2(500), Itmp, I_lastgeom
      Character*(linelen) szOPTARC, szTmp
      Logical OPTARC_presnt 
c
      Call A2getrec(10, 'JOBARC', 'IFLAGS  ', 100, Iflag_store1)
      Call A2getrec(10, 'JOBARC', 'IFLAGS2 ', 500, Iflag_store2)
      Call A2getrec(10,  'JOBARC', 'HAVEGEOM', 1, I_havegeom)
c
      Call gfname('OPTARC',szOPTARC,iOPTARC)
      Inquire(file=szOPTARC(1:iOPTARC),exist=OPTARC_presnt)
c

      If (I_havegeom .GT. 0 .AND. OPTARC_presnt .AND. 
     &    Iflags(54) .GT. 0)
     &    Call A2putrec(10, 'JOBARC', 'VIB_POPT', 1, 1)
c
      Iflag_store2(5) = 0
      Call A2getrec(10, 'JOBARC', 'NFDIRREP', 1, ITmp)
      If (Itmp .NE. 0) Itmp = 0
      Call A2putrec(10, 'JOBARC', 'NFDIRREP', 1, Itmp)
      Call A2getrec(10, 'JOBARC', 'FNDFDONE', 1, iTmp)
      If (Itmp .NE. 0) Itmp = 0
      Call A2putrec(10, 'JOBARC', 'FNDFDONE', 1, iTmp)
c
c 
c   o Prepare the derivative level for vibrational frequency.
c
      If (Iflags(54) .GE. 3) 
     &    Iflag_store1(3) = 1
c
c   o VIB=EXACT need especial care
c
      If (Iflags(2) .EQ. 1 .AND. Iflags(54) 
     &   .EQ. 1) Iflag_store1(3) = 2
c
      If ((Iflags(54) .EQ. 1)        .AND. 
     &    (Iflags2(5) .NE. 0) .AND. 
     &     I_havegeom .GT. 0) 
     &     Iflag_store1(3) = 2
c
c   o Set the step size for the numerical vibrational frequency calculation
c
      If (Iflags(54) .GE. 3) Then
         If (Iflags2(138).eq.2) then
c          o numerical gradients
                Iflag_store1(57) = 200
         Else
                Iflag_store1(57) = 50
         Endif
      Endif

       Print*, "The fd_stepsize is geting updated (new value)"
       Print*, Iflag_store1(57)
c
      Call A2putrec(10, 'JOBARC', 'IFLAGS  ', 100, Iflag_store1)
      Call A2putrec(10, 'JOBARC', 'IFLAGS2 ', 500, Iflag_store2)
c
c    o Truncate the JOBARC from JODAOUT for so that lower symmetry
c      points can pproperly run
c
      Call Aces_ja_init
      Call Aces_ja_truncate('JODAOUT ',1)
      Call Aces_ja_fin

      If ((Iflags(54) .GE. 3) .AND. 
     &     Iflags2(138) .eq. 2) Then
           Call A2putrec(1, 'JOBARC', 'FNDFDONE', 1, 0)
           Call A2putrec(1, 'JOBARC', 'POSTOPFD', 1, 1)
      Endif
c
      Call Runit("xjoda") 
c
      If ((Iflags(54) .GE. 3) .AND. 
     &     Iflags2(138) .eq. 2) Then
           Call A2putrec(1, 'JOBARC', 'POSTOPFD', 1, 0)
      Endif
c
c    o Restore the records that got deleted for future use.
c
      Call A2putrec(10, 'JOBARC', 'LASTGEOM', 1, 0)
      Call A2putrec(10, 'JOBARC', 'HAVEGEOM', 1, I_havegeom)
c
      Call Runit("rm -f OPTARC")
c
      Return
      End

