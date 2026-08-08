










      Subroutine Form_t2abef(T2aa,T2bb,T2ab,W_aa,W_bb,W_ab,
     +                       T2abef_1,T2abef_2,T2abef_3,
     +                       Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension T2abef_1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension T2abef_2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension T2abef_3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Integer M,N,A,B,E,F

      Data Dzero,Done,Quart,Half/0.0D0,1.0D0,0.250D0,0.50D0/

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




C T2abef_1(AB,EF) = 1/4T2(AB,MN)*W(MN,EF) 

      Scale = Dzero
      Fact1 = Quart 
      Fact2 = Half 
 
      Do F = 1, Nvrt_a
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
CSSS      T2abef_1(A,B,E,F) = Scale*T2abef_1(A,B,E,F)
      T2abef_1(A,B,E,F) = Dzero
      Do M = 1, Nocc_a
      Do N = 1, Nocc_a
      T2abef_1(A,B,E,F) = T2abef_1(A,B,E,F) + 
     +                    T2aa(A,B,M,N)*W_aa(M,N,E,F)*Fact1
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C T2abef_2(ab,ef) = 1/4T2(ab,mn)*W(mn,ef) 

      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
CSSS      T2abef_2(a,b,e,f) = Scale*T2abef_2(a,b,e,f)
      T2abef_2(a,b,e,f) = Dzero 
      Do M = 1, Nocc_b
      Do N = 1, Nocc_b
      T2abef_2(a,b,e,f) = T2abef_2(a,b,e,f) + 
     +                    T2bb(a,b,m,n)*W_bb(m,n,e,f)*Fact1
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C T2abef_3(Ab,Ef) = 1/2T2(Ab,Mn)*W(Mn,Ef)

      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
CSSS      T2abef_3(A,b,E,F) = Scale*T2abef_3(A,b,E,F)
      T2abef_3(A,b,E,F) = Dzero 
      Do M = 1, Nocc_a
      Do N = 1, Nocc_b
      T2abef_3(A,b,E,f)= T2abef_3(A,b,E,f) + 
     +                   T2ab(A,b,M,n)*W_ab(M,n,E,f)*Fact2
      Enddo    
      Enddo 
      Enddo
      Enddo
      Enddo
      Enddo

      Return
      End
