










      subroutine Set_active_space(Work,Maxcor,Nocc_a,Nocc_b,Nvrt_a,
     +                            Nvrt_b,Length)
     
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Integer Act_start

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




      L1_aa  = Nocc_a*Nvrt_a
      L1_bb  = Nocc_b*Nvrt_b
      L2_aa  = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L2_bb  = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L2_ab  = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      Length = L1_aa+L1_bb+L2_aa+L2_bb+L2_ab

      Ibgn  = 1
      Iend  = Ibgn + Length 
   
      If (Iend .Gt. Maxcor)
     +   Call Insmem("@-set_active_space",Iend,Maxcor)
     
      Call Form_exc_mask(Work(Ibgn),Length,Nocc_a,Nocc_b,Nvrt_a,
     +                   Nvrt_b)

      Sum = 0.0D0
      Do i=1,Length
         Sum = Sum + Work(Ibgn+i-1)*Work(Ibgn+i-1)
      Enddo
      Write(6,"(a,1x,2(1x,ES8.2E2))") " @-set_active_space,P*P     = ",
     +                                  Sum
      Write(6,*) 

      Return
      End
