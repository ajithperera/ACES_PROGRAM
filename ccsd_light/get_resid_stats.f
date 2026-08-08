










      Subroutine Get_resid_stats(T1resid_aa,T1resid_bb,T2resid_aa,
     +                           T2resid_bb,T2resid_ab,Tresid_rms,
     +                           Tresid_max,Nocc_a,Nocc_b,Nvrt_a,
     +                           Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Data Izero /0/

c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end

      Integer cc_maxcyc
      Integer Act_min_a,Act_min_b,Act_max_a,Act_max_b
      Integer Lineq_mxcyc
      Logical Ring_cc,Brueck,Active_space,Regular
      Double Precision ocn_oa,Ocn_ob,Ocn_va,Ocn_vb
      Double Precision Denom_tol,Brueck_tol,Lineq_tol
      Double Precision Rfac
      Dimension E_corr(0:500)

      Common /ccsdlight_vars/Ring_cc,Brueck,cc_conv,cc_maxcyc,
     +                       ocn_oa(Maxbasfn),ocn_ob(Maxbasfn),
     +                       ocn_va(Maxbasfn),ocn_vb(Maxbasfn),
     +                       E_corr,Denom_tol,Brueck_tol,Lineq_tol,
     +                       Act_min_a,Act_min_b,Act_max_a,
     +                       Act_max_b,Active_space,Lineq_mxcyc,
     +                       Regular,Rfac
     +                       




      Dimension T1resid_aa(Nvrt_a*Nocc_a)
      Dimension T1resid_bb(Nvrt_b*Nocc_b)
      Dimension T2resid_aa(Nvrt_a*Nvrt_a*Nocc_a*Nocc_a)
      Dimension T2resid_bb(Nvrt_b*Nvrt_b*Nocc_b*Nocc_b)
      Dimension T2resid_ab(Nvrt_a*Nvrt_b*Nocc_a*Nocc_b)

      Laa   = Nvrt_a*Nocc_a
      Lbb   = Nvrt_b*Nocc_b
      Laaaa = Nvrt_a*Nvrt_a*Nocc_a*Nocc_a
      Lbbbb = Nvrt_b*Nvrt_b*Nocc_b*Nocc_b
      Labab = Nvrt_a*Nvrt_b*Nocc_a*Nocc_b

      Dt1_res_aa = Ddot(Laa,T1resid_aa,1,T1resid_aa,1)
      Dt1_res_bb = Ddot(Lbb,T1resid_bb,1,T1resid_bb,1)
      Dt2_res_aa = Ddot(Laaaa,T2resid_aa,1,T2resid_aa,1)
      Dt2_res_bb = Ddot(Lbbbb,T2resid_bb,1,T2resid_bb,1)
      Dt2_res_ab = Ddot(Labab,T2resid_ab,1,T2resid_ab,1)
      
      If (Laa .Gt. Izero) Dt1_res_aa = Dsqrt(Dt1_res_aa/Laa)
      If (Lbb .Gt. Izero) Dt1_res_bb = Dsqrt(Dt1_res_bb/Lbb)
      If (Laaaa .Gt. Izero) Dt2_res_aa = Dsqrt(Dt2_res_aa/Laaaa)
      If (Lbbbb .Gt. Izero) Dt2_res_bb = Dsqrt(Dt2_res_bb/Lbbbb)
      If (Labab .Gt. Izero) Dt2_res_ab = Dsqrt(Dt2_res_ab/Labab)

      Tresid_rms = Max(Dt1_res_aa,Dt1_res_bb,Dt2_res_aa,Dt2_res_bb,
     +             Dt2_res_ab)

      Write(6,"(a,F15.10)") "The root mean square of the residual:",
     +                       Tresid_rms

      Dt1_res_aa = T1resid_aa(Isamax(Laa,T1resid_aa,1))
      Dt1_res_bb = T1resid_bb(Isamax(Lbb,T1resid_bb,1))
      Dt2_res_aa = T2resid_aa(Isamax(Laaaa,T2resid_aa,1))
      Dt2_res_bb = T2resid_ab(Isamax(Lbbbb,T2resid_bb,1))
      Dt2_res_aa = T2resid_ab(Isamax(Labab,T2resid_ab,1))
    
      
      Tresid_max = Max(Dt1_res_aa,Dt1_res_bb,Dt2_res_aa,Dt2_res_bb,
     +                 Dt2_res_ab)
      Write(6,"(a,F15.10)") "The maximum element of the residual :",
     +                       Tresid_max

      Return
      End 
