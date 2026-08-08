










      Subroutine Form_hbabci_term5(Hbar_abci1,Hbar_abci2,Hbar_abci3,
     +                             Hbar_abci4,Wtt_mbej1,Wtt_mbej2,
     +                             Wtt_mbej3,Wtt_mbej4,Wtt_mbej5,
     +                             Wtt_mbej6,Taa,Tbb,Nocc_a,Nocc_b,
     +                             Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension Taa(Nvrt_a,Nocc_a)
      Dimension Tbb(Nvrt_b,Nocc_b)

      Dimension Hbar_abci1(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Hbar_abci2(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Hbar_abci3(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Hbar_abci4(Nvrt_b,Nvrt_a,Nvrt_b,Nocc_a)

      Dimension Wtt_mbej1(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Wtt_mbej2(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Wtt_mbej3(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Wtt_mbej4(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension Wtt_mbej5(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension Wtt_mbej6(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Integer A,B,C,I,M,N

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
      Do M = 1, Nocc_a
         T = T + Taa(A,M)*Wtt_mbej1(M,B,C,I)
      Enddo
         Hbar_abci1(A,B,C,I) = Hbar_abci1(A,B,C,I) - T
         Hbar_abci1(B,A,C,I) = Hbar_abci1(B,A,C,I) + T
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
      Do M = 1, Nocc_b
         T = T + Tbb(a,m)*Wtt_mbej2(m,b,c,i)
      Enddo
         Hbar_abci2(a,b,c,i) = Hbar_abci2(a,b,c,i) - T
         Hbar_abci2(b,a,c,i) = Hbar_abci2(b,a,c,i) + T
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
      Do M = 1, Nocc_a
         T = T + Taa(A,M)*Wtt_mbej3(M,b,C,i)
      Enddo
         Hbar_abci3(A,b,C,i) = Hbar_abci3(a,b,c,i) - T
      Enddo
      Enddo
      Enddo
      Enddo

      Do I = 1, Nocc_b
      Do C = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_b
         T = T + Tbb(b,m)*Wtt_mbej6(m,A,i,C)
      Enddo
         Hbar_abci3(A,b,C,i) = Hbar_abci3(A,b,C,i) - T
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
      Do M = 1, Nocc_b
         T = T + Tbb(a,m)*Wtt_mbej4(m,B,c,I)
      Enddo
         Hbar_abci4(a,B,c,I) = Hbar_abci4(a,B,c,I) - T
      Enddo
      Enddo
      Enddo
      Enddo

      Do I = 1, Nocc_a
      Do C = 1, Nvrt_b
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_b
         T = 0.0D0
      Do M = 1, Nocc_a
         T = T + Taa(B,M)*Wtt_mbej5(M,a,I,c)
      Enddo
         Hbar_abci4(a,B,c,I) = Hbar_abci4(a,B,c,I) - T
      Enddo
      Enddo
      Enddo
      Enddo

      Return
      End
