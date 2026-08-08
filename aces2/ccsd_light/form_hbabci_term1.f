










      Subroutine Form_hbabci_term1(Hbar_abci1,Hbar_abci2,Hbar_abci3,
     +                             Hbar_abci4,W4_aa,W4_bb,W4_ab,W4_ba,
     +                             Hbar_aibc1,Hbar_aibc2,Hbar_aibc3,
     +                             Hbar_aibc4,T2aa,T2bb,T2ab,Nocc_a,
     +                             Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension W4_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W4_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W4_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W4_ba(Nvrt_a,Nvrt_b,Nocc_a,Nvrt_b)

      Dimension Hbar_aibc1(Nvrt_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension Hbar_aibc2(Nvrt_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension Hbar_aibc3(Nvrt_a,Nocc_b,Nvrt_a,Nvrt_b)
      Dimension Hbar_aibc4(Nvrt_b,Nocc_a,Nvrt_b,Nvrt_a)

      Dimension Hbar_abci1(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Hbar_abci2(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Hbar_abci3(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Hbar_abci4(Nvrt_b,Nvrt_a,Nvrt_b,Nocc_a)

      Integer A,B,C,I,M,E

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




C AAAA contribution
      
      Do I = 1, Nocc_a
      Do C = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0
         D = (1.0D0-Ocn_va(A))
      Do E = 1, Nvrt_a
      Do M = 1, Nocc_a
         T = T + T2aa(B,E,I,M)*Hbar_aibc1(A,M,C,E)
      Enddo
      Enddo
         Hbar_abci1(A,B,C,I) = Hbar_abci1(A,B,C,I) + T*D
         Hbar_abci1(B,A,C,I) = Hbar_abci1(B,A,C,I) - T*D
      Enddo
      Enddo
      Enddo
      Enddo

      Do I = 1, Nocc_a
      Do C = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0
         D = (1.0D0-Ocn_va(A))
      Do E = 1, Nvrt_b
      Do M = 1, Nocc_b
         T = T + T2ab(B,e,I,m)*Hbar_aibc3(A,m,C,e)
      Enddo
      Enddo
         Hbar_abci1(A,B,C,I) = Hbar_abci1(A,B,C,I) + T*D
         Hbar_abci1(B,A,C,I) = Hbar_abci1(B,A,C,I) - T*D
      Enddo
      Enddo
      Enddo
      Enddo

C BBBB contribution

      Do I = 1, Nocc_b
      Do C = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T = 0.0D0
         D = (1.0D0-Ocn_vb(a))
      Do E = 1, Nvrt_b
      Do M = 1, Nocc_b
         T = T + T2bb(b,e,i,m)*Hbar_aibc2(a,m,c,e)
      Enddo
      Enddo
         Hbar_abci2(a,b,c,i) = Hbar_abci2(a,b,c,i) + T*D
         Hbar_abci2(b,a,c,i) = Hbar_abci2(b,a,c,i) - T*D
      Enddo
      Enddo
      Enddo
      Enddo

      Do I = 1, Nocc_b
      Do C = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T = 0.0D0
         D = (1.0D0-Ocn_vb(a))
      Do E = 1, Nvrt_a
      Do M = 1, Nocc_a
         T = T + T2ab(E,b,M,i)*Hbar_aibc4(a,M,c,E)
      Enddo
      Enddo
         Hbar_abci2(a,b,c,i) = Hbar_abci2(a,b,c,i) + T*D
         Hbar_abci2(b,a,c,i) = Hbar_abci2(b,a,c,i) - T*D
      Enddo
      Enddo
      Enddo
      Enddo

C ABAB contribution

      Do I = 1, Nocc_b
      Do C = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
         D = (1.0D0-Ocn_va(A))
      Do E = 1, Nvrt_b
      Do M = 1, Nocc_b
         T = T + T2bb(b,e,i,m)*Hbar_aibc3(A,m,C,e) 
      Enddo
      Enddo
        Hbar_abci3(A,b,C,i) = Hbar_abci3(A,b,C,i) + T*D
      Enddo
      Enddo
      Enddo
      Enddo

      Do I = 1, Nocc_b
      Do C = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
         D = (1.0D0-Ocn_va(A))
      Do E = 1, Nvrt_a
      Do M = 1, Nocc_a
         T = T + T2ab(E,b,M,i)*Hbar_aibc1(A,M,C,E)
      Enddo
      Enddo
         Hbar_abci3(A,b,C,i) = Hbar_abci3(A,b,C,i) + T*D
      Enddo
      Enddo
      Enddo
      Enddo

      Do I = 1, Nocc_b
      Do C = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
         D = (1.0D0-Ocn_vb(b))
      Do E = 1, Nvrt_b
      Do M = 1, Nocc_a
         T = T + T2ab(A,e,M,i)*Hbar_aibc4(b,M,e,C)
      Enddo
      Enddo
         Hbar_abci3(A,b,C,i) = Hbar_abci3(A,b,C,i) - T*D
      Enddo
      Enddo
      Enddo
      Enddo

C BABA contribution

      Do I = 1, Nocc_a
      Do C = 1, Nvrt_b
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_b
         T = 0.0D0
         D = (1.0D0-Ocn_vb(a))
      Do E = 1, Nvrt_a
      Do M = 1, Nocc_a
         T = T + T2aa(B,E,I,M)*Hbar_aibc4(a,M,c,E)
      Enddo
      Enddo
         Hbar_abci4(a,B,c,I) = Hbar_abci4(a,B,c,I) + T*D
      Enddo
      Enddo
      Enddo
      Enddo

      Do I = 1, Nocc_a
      Do C = 1, Nvrt_b
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_b
         T = 0.0D0
         D = (1.0D0-Ocn_vb(a))
      Do E = 1, Nvrt_b
      Do M = 1, Nocc_b
         T = T + T2ab(B,e,I,m)*Hbar_aibc2(a,m,c,e)
      Enddo
      Enddo
         Hbar_abci4(a,B,c,I) = Hbar_abci4(a,B,c,I) + T*D
      Enddo
      Enddo
      Enddo
      Enddo

      Do I = 1, Nocc_a
      Do C = 1, Nvrt_b
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_b
         T = 0.0D0
         D = (1.0D0-Ocn_va(B))
      Do E = 1, Nvrt_a
      Do M = 1, Nocc_b
         T = T + T2ab(E,a,I,m)*Hbar_aibc3(B,m,E,c)
      Enddo
      Enddo
         Hbar_abci4(a,B,c,I) = Hbar_abci4(a,B,c,I) - T*D 
      Enddo
      Enddo
      Enddo
      Enddo

      Return
      End
