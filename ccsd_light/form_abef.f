










      Subroutine Form_abef(T2aa,T2bb,T2ab,T1aa,T1bb,W2_aa,W2_bb,w2_ab,
     +                     W_aa,W_bb,W_ab,W4_aa,W4_bb,W4_ab,W4_ba,
     +                     Wabef_1,Wabef_2,Wabef_3,Nocc_a,Nocc_b,
     +                     Nvrt_a,Nvrt_b,NR)

      Implicit Double Precision (A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension W2_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension W2_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension W2_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Dimension W4_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W4_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W4_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W4_ba(Nvrt_a,Nvrt_b,Nocc_a,Nvrt_a)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension Wabef_1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension Wabef_2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension Wabef_3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Integer M,N,A,B,E,F
      Logical NR

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




      IF (NR) then
         Fact_aa = 0.50D0
         Fact_bb = 0.50D0
         Fact_ab = 1.00D0
      Else
         Fact_aa = 0.250D0
         Fact_bb = 0.250D0
         Fact_ab = 0.500D0
      Endif
      
C Wabef_1(AB,EF) = W2(AB,EF) + 1/4T2(AB,MN)*W(MN,EF) 
 
      Do F = 1, Nvrt_a
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         C = (1.0D0-Ocn_va(A))*(1.0D0-Ocn_va(B))
         Wabef_1(A,B,E,F) =  W2_aa(A,B,E,F)*C
      Do M = 1, Nocc_a
      Do N = 1, Nocc_a
CSSS      C = Ocn_oa(M)*Ocn_oa(N)
      C = 1.0D0
      Wabef_1(A,B,E,F) = Wabef_1(A,B,E,F) + 
     +                   Fact_aa*T2aa(A,B,M,N)*W_aa(M,N,E,F)*C 
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_1(A,B,E,F) = -P_(A,B)T1(B,M)W(E,F,A,M)

      Do F = 1, Nvrt_a
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
CSSS         C =  Ocn_oa(M)
         C = (1.0D0-Ocn_va(A))
         T = T + T1aa(B,M)*W4_aa(E,F,A,M)*C
      Enddo
         Wabef_1(A,B,E,F) =  Wabef_1(A,B,E,F) - T
         Wabef_1(B,A,E,F) =  Wabef_1(B,A,E,F) + T
      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_2(ab,ef) = W2(ab,ef) + 1/4T2(ab,mn)*W(mn,ef) 

      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         C = (1.0D0-Ocn_vb(a))*(1.0D0-Ocn_vb(b))
         Wabef_2(a,b,e,f) = W2_bb(a,b,e,f)*C
      Do M = 1, Nocc_b
      Do N = 1, Nocc_b
CSSS      C = Ocn_ob(m)*Ocn_ob(n)
      C = 1.0D0
      Wabef_2(a,b,e,f) = Wabef_2(a,b,e,f) + 
     +                   Fact_bb*T2bb(a,b,m,n)*W_bb(m,n,e,f)*C 
      Enddo
      Enddo

      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_2(a,b,e,f) = -P_(a,b)T1(b,m)*W(e,f,a,m)

      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T = 0.0D0
      Do M = 1, Nocc_b
CSSS         C =  Ocn_ob(m)
         C = (1.0D0-Ocn_vb(a))
         T = T + T1bb(b,m)*W4_bb(e,f,a,m)*C
      Enddo
         Wabef_2(a,b,e,f) =  Wabef_2(a,b,e,f) - T
         Wabef_2(b,a,e,f) =  Wabef_2(b,a,e,f) + T
      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_3(Ab,Ef) = W2(Ab,Ef) + T2(Ab,Mn)*W(Mn,Ef)

      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         C = (1.0D0-Ocn_va(A))*(1.0D0-Ocn_vb(b))
         Wabef_3(A,b,E,F) = W2_ab(A,b,E,f)*C

      Do M = 1, Nocc_a
      Do N = 1, Nocc_b
      Wabef_3(A,b,E,f)= Wabef_3(A,b,E,f) + 
     +                  Fact_ab*T2ab(A,b,M,n)*W_ab(M,n,E,f)
      Enddo    
      Enddo 

      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_3(Ab,Ef) = -T1(b,m)*W(E,f,A,m) 
      
      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a  
         T = 0.0D0
      Do M = 1, Nocc_b
CSSS         C = Ocn_ob(m)
         C = (1.0D0-Ocn_va(A))
         T = T + T1bb(b,m)*W4_ab(E,f,A,m)*C
      Enddo
         Wabef_3(A,b,E,F) = Wabef_3(A,b,E,F) - T
      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_3(Ab,Ef) = -T1(A,M)*W(E,f,b,M)->T1(A,M)*W(E,f,M,b)

      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
CSSS         C = Ocn_oa(M)
         C = (1.0D0-Ocn_vb(b))
         T = T + T1aa(A,M)*W4_ba(E,f,M,b)*C
      Enddo
         Wabef_3(A,b,E,F) = Wabef_3(A,b,E,F) - T
      Enddo
      Enddo
      Enddo
      Enddo

      L_aaaa = Nvrt_a*Nvrt_a*Nvrt_a*Nvrt_a
      L_bbbb = Nvrt_b*Nvrt_b*Nvrt_b*Nvrt_b
      L_abab = Nvrt_a*Nvrt_b*Nvrt_a*Nvrt_b
      call checksum("Form_wABEF:",Wabef_1,L_aaaa)
      call checksum("Form_wabef:",Wabef_2,L_bbbb)
      call checksum("Form_wAbEf:",Wabef_3,L_abab)

      Return
      End
