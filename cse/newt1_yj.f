










      Subroutine Newt1_yj(T1aa_new,T1bb_new,T1resid_aa,T1resid_bb,
     +                    H1bar_aa,H1bar_bb,Dt1_aa,Dt1_bb,Dt1,Nocc_a,
     +                    Nocc_b,Nvrt_a,Nvrt_b,Nbasis)

      Implicit Double Precision(A-H,O-Z)
      
      Dimension T1aa_new(Nvrt_a*Nocc_a)
      Dimension T1bb_new(Nvrt_b*Nocc_b)

      Dimension Dt1_aa(Nvrt_a*Nocc_a)
      Dimension Dt1_bb(Nvrt_b*Nocc_b)

      Dimension T1resid_aa(Nvrt_a,Nocc_a)
      Dimension T1resid_bb(Nvrt_b,Nocc_b)

      Dimension H1bar_aa(Nvrt_a*Nocc_a,Nvrt_a*Nocc_a)
      Dimension H1bar_bb(Nvrt_b*Nocc_b,Nvrt_b*Nocc_b)

      Data Onem,Onep,Dzro/-1.0D0,1.0D0,0.0D0/

      Integer A, B, I, J

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




      L_aa = Nocc_a*Nvrt_a
      L_bb = Nocc_b*Nvrt_b

      Call Dscal(L_aa,Onem,T1resid_aa,1)
      Call Dscal(L_bb,Onem,T1resid_bb,1)

      Call Dgemm("N","N",L_aa,1,L_aa,Onep,H1bar_aa,L_aa,T1resid_aa,
     +            L_aa,Dzro,Dt1_aa,L_aa)
      Call Dgemm("N","N",L_bb,1,L_bb,Onep,H1bar_bb,L_bb,T1resid_bb,
     +            L_bb,Dzro,Dt1_bb,L_bb)

C Construct t1(new) = t1(old) + dt1

      Call Daxpy(L_aa,Onep,DT1_aa,1,T1aa_new,1)
      Call Daxpy(L_bb,Onep,DT1_bb,1,T1bb_new,1)

      L_aa = Nocc_a*NvRT_a
      L_bb = Nocc_b*Nvrt_b
      call checksum("DT1aa     :",DT1_aa,L_aa)
      call checksum("DT1bb     :",DT1_bb,L_bb)
      Trace = 0.0D0
      do i=1,Nocc_a
      Index = (I-1)*Nvrt_a + I
      Trace = Trace + (T1aa_new(Index) + ocn_oa(i))
      enddo
      Write(6,"(a,F12.8)") "Trace of the density :", Trace   
      DT1_alpha = Ddot(L_aa,DT1_aa,1,DT1_aa,1)
      DT1_Beta  = Ddot(L_bb,DT1_bb,1,DT1_bb,1)

      DT1 = Max(DT1_alpha,DT1_Beta)
      DT1 = Dsqrt(Dabs(DT1))

      Return
      End
