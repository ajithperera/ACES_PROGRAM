










      Subroutine T2_inmbej(T2aa,T2bb,T2ab,W_aa,W_bb,W_ab,Wmbej_1,
     +                     Wmbej_2,Wmbej_3,Wmbej_4,Wmbej_5,Wmbej_6,
     +                     Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_b)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension Wmbej_1(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Wmbej_2(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Wmbej_3(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Wmbej_4(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension Wmbej_5(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension Wmbej_6(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Integer M,N,A,F,E,B,J

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




C Wmbej_1(MB,EJ) <- - 1/2T2(FB,JN)*W2(MN,EF) - 1/2T2(fB,Jn)*W2(Mn,Ef)
 
      Do J = 1, Nocc_a
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do M = 1, Nocc_a
         T = 0.0
      Do F = 1, Nvrt_a
      Do N = 1, Nocc_a
         T = T - 0.50D0*T2aa(F,B,J,N)*W_aa(M,N,E,F)
      Enddo
      Enddo
         Wmbej_1(M,B,E,J) = T
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_a
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do M = 1, Nocc_a
      Do F = 1, Nvrt_b
      Do N = 1, Nocc_b
         Wmbej_1(M,B,E,J)= Wmbej_1(M,B,E,J) + 
     +                     0.50D0*T2ab(B,f,J,n)*W_ab(M,n,E,f)
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C Wmbej_2(mb,ej) <- - 1/2T2(fb,jn)*W2(mn,ef) + 1/2T2(Fb,jN)*W2(mN,eF)

      Do J = 1, Nocc_b
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_b
         T = 0.0D0
      Do F = 1, Nvrt_b
      Do N = 1, Nocc_b
         T = T - 0.50D0*T2bb(f,b,j,n)*W_bb(m,n,e,f)
      Enddo
      Enddo
         Wmbej_2(m,b,e,j) = T
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_b
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_b
      Do F = 1, Nvrt_a
      do N = 1, Nocc_a
         Wmbej_2(m,b,e,j)= Wmbej_2(m,b,e,j) + 
     +                     0.50D0*T2ab(F,b,N,j)*W_ab(N,m,F,e)
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C Wmbej_3(Mb,Ej) <- - 1/2T2(fb,jn)*W(Mn,Ef) + 1/2T2(Fb,jN)*W(MN,EF)

      Do J = 1, Nocc_b 
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_a
         T = 0.0D0
      Do F = 1, Nvrt_b
      Do N = 1, Nocc_b
         T = T  + 0.50D0*T2bb(b,f,j,n)*W_ab(M,n,E,f)
      Enddo
      Enddo
        Wmbej_3(M,b,E,j) = T
      Enddo
      Enddo
      Enddo
      Enddo
  
      Do J = 1, Nocc_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_a
      Do F = 1, Nvrt_a
      Do N = 1, Nocc_a
         Wmbej_3(M,b,E,j)= Wmbej_3(M,b,E,j) +  
     +                     0.50D0*T2ab(F,b,N,j)*W_aa(M,N,E,F)
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C Wmbej_4(mB,eJ) <- - 1/2T2(BF,JN)*W2(mN,eF) + 1/2T2(Bf,Jn)*W2(mn,ef)

      Do J = 1, Nocc_a
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_a
      Do M = 1, Nocc_b
         T = 0.0D0
      Do F = 1, Nvrt_a
      Do N = 1, Nocc_a
         T = T + 0.50D0*T2aa(B,F,J,N)*W_ab(N,m,F,e)
      Enddo
      Enddo
         Wmbej_4(m,B,e,J) = T
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_a
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_a
      Do M = 1, Nocc_b
      Do F = 1, Nvrt_b
      Do N = 1, Nocc_b
         Wmbej_4(m,B,e,J) = Wmbej_4(m,B,e,J) + 
     +                      0.50D0*T2ab(B,f,J,n)*W_bb(m,n,e,f)
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C Wmbej_5(Mb,Je) <-  + 1/2T2(bF,JN)*W2(Mn,Fe) 

      Do J = 1, Nocc_a
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_a
         T = 0.0D0
      Do F = 1, Nvrt_a
      Do N = 1, Nocc_b
         T = T - 0.50D0*T2ab(F,b,J,n)*W_ab(M,n,F,e)
      Enddo
      Enddo
         Wmbej_5(M,b,J,e) = T
      Enddo
      Enddo
      Enddo
      Enddo

C Wmbej_6(mB,jE) <-  + 1/2T2(fB,jN)*W2(mN,fE) 

      Do J = 1, Nocc_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a 
      Do M = 1, Nocc_b
         T = 0.0D0
      Do F = 1, Nvrt_b
      Do N = 1, Nocc_a
         T = T - 0.50D0*T2ab(B,f,N,j)*W_ab(N,m,E,f)
      Enddo
      Enddo
         Wmbej_6(m,B,j,E) = T 
      Enddo
      Enddo
      Enddo
      Enddo

      Write(6,"(a)") " @-T2_inmbej"
      L_aaaa = Nocc_a*Nvrt_a*Nvrt_a*Nocc_a
      L_bbbb = Nocc_b*Nvrt_b*Nvrt_b*Nocc_b
      L_abab = Nocc_a*Nvrt_b*Nvrt_a*Nocc_b
      L_baba = Nocc_b*Nvrt_a*Nvrt_b*Nocc_a
      L_abba = Nocc_a*Nvrt_b*Nocc_a*Nvrt_b
      L_baab = Nocc_b*Nvrt_a*Nocc_b*Nvrt_a
      call checksum("T2->MBEJ  :",Wmbej_1,L_aaaa)
      call checksum("T2->mbej  :",Wmbej_2,L_bbbb)
      call checksum("T2->MbEj  :",Wmbej_3,L_abab)
      call checksum("T2->mBeJ  :",Wmbej_4,L_baba)
      call checksum("T2->mBEj  :",Wmbej_5,L_abba)
      call checksum("T2->MbeJ  :",Wmbej_6,L_baab)
      Write(6,*)

      Return
      End
