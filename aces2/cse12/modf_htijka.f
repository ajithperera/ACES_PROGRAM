










      Subroutine Modf_htijka(T1aa,T1bb,W5_aa,W5_bb,W5_ab,W5_ba,
     +                       W_aa,W_bb,W_ab,Nocc_a,Nocc_b,Nvrt_a,
     +                       Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension W5_aa(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension W5_bb(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension W5_ab(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension W5_ba(Nocc_a,Nocc_b,Nvrt_a,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Integer M,N,I,E,F

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



 
      Do F = 1, Nvrt_a
      Do I = 1, Nocc_a
      Do N = 1, Nocc_a
      Do M = 1, Nocc_a 
         T = 0.0D0
      DO E = 1, Nvrt_a
         If (Ocn_oa(I) .NE. 0.0D0) Then
            C = 1.0D0/Ocn_oa(I)
         Else
            C = 1.0D0
         Endif
         T = T + T1aa(E,I)*W_aa(M,N,E,F)*C
      Enddo 
         W5_aa(M,N,I,F) = W5_aa(M,N,I,F) - T
      Enddo 
      Enddo
      Enddo 
      Enddo

      Do F = 1, Nvrt_b
      Do I = 1, Nocc_b
      Do N = 1, Nocc_b
      Do M = 1, Nocc_b
         T = 0.0D0
      DO E = 1, Nvrt_b
         If (Ocn_ob(i) .NE. 0.0D0) Then
            C = 1.0D0/Ocn_ob(i)
         Else
            C = 1.0D0
         Endif
         T = T + T1bb(e,i)*W_bb(m,n,e,f)*C 
      Enddo  
         W5_bb(m,n,i,f) = W5_bb(m,n,i,f) - T
      Enddo
      Enddo
      Enddo
      Enddo

      Do F = 1, Nvrt_b
      Do I = 1, Nocc_a
      Do N = 1, Nocc_b
      Do M = 1, Nocc_a
         T = 0.0D0
      DO E = 1, Nvrt_a
         If (Ocn_oa(I) .NE. 0.0D0) Then
            C = 1.0D0/Ocn_oa(I)
         Else
            C = 1.0D0
         Endif
         T = T + T1aa(E,I)*W_ab(M,n,E,f)*C 
      Enddo  
         W5_ab(M,n,I,f) = W5_ab(M,n,I,f) - T
      Enddo
      Enddo
      Enddo
      Enddo

      Do I = 1, Nocc_b
      Do F = 1, Nvrt_a
      Do N = 1, Nocc_b
      Do M = 1, Nocc_a
         T = 0.0D0
      DO E = 1, Nvrt_b
         If (Ocn_ob(i) .NE. 0.0D0) Then
            C = 1.0D0/Ocn_ob(i)
         Else
            C = 1.0D0
         Endif
         T = T + T1bb(e,i)*W_ab(M,n,F,e)*C 
      Enddo  
         W5_ba(M,n,F,i) = W5_ba(M,n,F,i) - T
      Enddo
      Enddo
      Enddo
      Enddo

      Return
      End 

