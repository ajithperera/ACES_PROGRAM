










      Subroutine Get_oneh(Oneh1,Oneh2,Oneh_aa,Oneh_bb,Scr1,
     +                    Scr2,E1_aa,E1_bb,Nbasis)

      Implicit Double Precision(A-H,O-Z)

      Dimension Oneh1(Nbasis*(Nbasis+1)/2)
      Dimension Oneh2(Nbasis,Nbasis) 
      Dimension Oneh_aa(Nbasis,Nbasis)
      Dimension Oneh_bb(Nbasis,Nbasis)
 
      Dimension Scr1(Nbasis*Nbasis),Scr2(Nbasis*Nbasis)

c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end

      Integer cc_maxcyc
      Integer Act_min_a,Act_min_b,Act_max_a,Act_max_b
      Logical Ring_cc,Brueck,Active_space
      Double Precision ocn_oa,Ocn_ob,Ocn_va,Ocn_vb
      Double Precision Denom_tol,Brueck_tol
      Dimension E_corr(0:500)

      Common /ccsdlight_vars/Ring_cc,Brueck,cc_conv,cc_maxcyc,
     +                       ocn_oa(Maxbasfn),ocn_ob(Maxbasfn),
     +                       ocn_va(Maxbasfn),ocn_vb(Maxbasfn),
     +                       E_corr,Denom_tol,Brueck_tol,
     +                       Act_min_a,Act_min_b,Act_max_a,
     +                       Act_max_b,Active_space




      Ldim = Nbasis*(Nbasis+1)/2
      Call Getrec(20,"JOBARC","ONEHAO  ",Ldim,Oneh1)
      Call Expnd2(Oneh1,Oneh2,Nbasis)


      Nbasis2 = Nbasis*Nbasis 
      If (Iend .GT. Maxcor) Call Insmem("prep_fock",Iend,Maxcor)

      Call Form_mo_fock(Oneh2,Oneh_aa,Scr1,Scr2,Nbasis,Nbasis,1)
      Call Form_mo_fock(Oneh2,Oneh_bb,Scr1,Scr2,Nbasis,Nbasis,2)

      E1_aa = 0.0D0
      DO I = 1, Nbasis 
         C = Ocn_oa(I)
         E1_aa = E1_aa + Oneh_aa(I,I)*C
      Enddo

      E1_bb = 0.0D0
      DO I = 1, Nbasis 
         C = Ocn_ob(i)
         E1_bb = E1_bb + oneh_bb(i,i)*C
      Enddo


      Return
      End
