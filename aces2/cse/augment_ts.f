










      Subroutine Augment_ts(T1aa,T1bb,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,
     +                      Daa,Dbb)

      Implicit Double Precision(A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension Daa(Nvrt_a,Nocc_a)
      Dimension Dbb(Nvrt_b,Nocc_b)

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




      Integer I,J,A,B

      Laa = Nvrt_a*Nocc_a
      LBb = Nvrt_b*Nocc_b
      Call Dzero(Daa, Laa)
      Call Dzero(Daa, Lbb)

      Index = 0
      Do I = 1, Nocc_a
         Index = Index + 1
      Do A = 1, Nvrt_a
         Daa(I,I) = Ocn_oa(index)
      Enddo
      Enddo

      Index = 0
      Do I = 1, Nocc_b
         Index = Index + 1
      Do A = 1, Nvrt_b
         Dbb(I,I) = Ocn_ob(index)
      Enddo
      Enddo

      Write(6,"(a)") "The Alpha and Beta density matrices" 
      Call output(Daa,1,Nocc_a,1,Nvrt_a,Nocc_a,Nvrt_a,1)
      Write(6,*)
      Call output(Dbb,1,Nocc_a,1,Nvrt_a,Nocc_a,Nvrt_a,1)

      Return
      End
