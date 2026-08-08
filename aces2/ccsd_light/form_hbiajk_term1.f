










      Subroutine Form_hbiajk_term1(Hbar_iajk1,Hbar_iajk2,Hbar_iajk3,
     +                             Hbar_iajk4,W5_aa,W5_bb,W5_ab,W5_ba,
     +                             T2aa,T2bb,T2ab,Nocc_a,Nocc_b,Nvrt_a,
     +                             Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension W5_aa(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension W5_bb(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension W5_ab(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension W5_ba(Nocc_a,Nocc_b,Nvrt_a,Nocc_b)

      Dimension Hbar_iajk1(Nocc_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Hbar_iajk2(Nocc_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Hbar_iajk3(Nocc_a,Nvrt_b,Nocc_a,Nocc_b)
      Dimension Hbar_iajk4(Nocc_b,Nvrt_a,Nocc_b,Nocc_a)

      Integer I,J,K,M,E,A

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
      
      Do K = 1, Nocc_a
      Do J = 1, Nocc_a
      Do A = 1, Nvrt_a
      Do I = 1, Nocc_a
         T = 0.0D0
         C = Ocn_oa(J)
      Do E = 1, Nvrt_a
      Do M = 1, Nocc_a
         T = T + T2aa(A,E,K,M)*W5_aa(I,M,J,E)*C
      Enddo
      Enddo
         Hbar_iajk1(I,A,J,K) = Hbar_iajk1(I,A,J,K) + T
         Hbar_iajk1(I,A,K,J) = Hbar_iajk1(I,A,K,J) - T
      Enddo
      Enddo
      Enddo
      Enddo

      Do K = 1, Nocc_a
      Do J = 1, Nocc_a
      Do A = 1, Nvrt_a
      Do I = 1, Nocc_a
         T = 0.0D0
         C = Ocn_oa(J)
      Do E = 1, Nvrt_b
      Do M = 1, Nocc_b
         T = T + T2ab(A,e,K,m)*W5_ab(I,m,J,e)*C
      Enddo
      Enddo
         Hbar_iajk1(I,A,J,K) = Hbar_iajk1(I,A,J,K) + T
         Hbar_iajk1(I,A,K,J) = Hbar_iajk1(I,A,K,J) - T
      Enddo
      Enddo
      Enddo
      Enddo

C BBBB contribution

      Do K = 1, Nocc_b
      Do J = 1, Nocc_b
      Do A = 1, Nvrt_b
      Do I = 1, Nocc_b
         T = 0.0D0
         C = Ocn_ob(j)
      Do E = 1, Nvrt_b
      Do M = 1, Nocc_b
         T = T + T2bb(a,e,k,m)*W5_bb(i,m,j,e)*C
      Enddo
      Enddo
         Hbar_iajk2(i,a,j,k) = Hbar_iajk2(i,a,j,k) + T
         Hbar_iajk2(i,a,k,j) = Hbar_iajk2(i,a,k,j) - T
      Enddo
      Enddo
      Enddo
      Enddo

      Do K = 1, Nocc_b
      Do J = 1, Nocc_b
      Do A = 1, Nvrt_b
      Do I = 1, Nocc_b
         T = 0.0D0
         C = Ocn_ob(j)
      Do E = 1, Nvrt_a
      Do M = 1, Nocc_a
         T = T + T2ab(E,a,M,k)*W5_ba(M,i,E,j)*C
      Enddo
      Enddo
         Hbar_iajk2(i,a,j,k) = Hbar_iajk2(i,a,j,k) + T
         Hbar_iajk2(i,a,k,j) = Hbar_iajk2(i,a,k,j) - T
      Enddo
      Enddo
      Enddo
      Enddo

C ABAB contribution

      Do K = 1, Nocc_b
      Do J = 1, Nocc_a
      Do A = 1, Nvrt_b
      Do I = 1, Nocc_a
         T = 0.0D0
         C = Ocn_oa(J)
      Do E = 1, Nvrt_b
      Do M = 1, Nocc_b
         T = T + T2bb(a,e,k,m)*W5_ab(I,m,J,e)*C
      Enddo
      Enddo
         Hbar_iajk3(I,a,J,k) = Hbar_iajk3(I,a,J,k) + T
      Enddo
      Enddo
      Enddo
      Enddo

      Do K = 1, Nocc_b
      Do J = 1, Nocc_a
      Do A = 1, Nvrt_b
      Do I = 1, Nocc_a
         T = 0.0D0
         C = Ocn_oa(J)
      Do E = 1, Nvrt_a
      Do M = 1, Nocc_a
         T = T + T2ab(E,a,M,k)*W5_aa(I,M,J,E)*C
      Enddo
      Enddo
         Hbar_iajk3(I,a,J,k) = Hbar_iajk3(I,a,J,k) + T
      Enddo
      Enddo
      Enddo
      Enddo

      Do K = 1, Nocc_b
      Do J = 1, Nocc_a
      Do A = 1, Nvrt_b
      Do I = 1, Nocc_a
         T = 0.0D0
         C = Ocn_ob(k)
      Do E = 1, Nvrt_a
      Do M = 1, Nocc_b
         T = T + T2ab(E,a,J,m)*W5_ba(I,m,E,k)*C
      Enddo
      Enddo
         Hbar_iajk3(I,a,J,k) = Hbar_iajk3(I,a,J,k) - T
      Enddo
      Enddo
      Enddo
      Enddo

C BABA contribution

      Do K = 1, Nocc_a
      Do J = 1, Nocc_b
      Do A = 1, Nvrt_a
      Do I = 1, Nocc_b
         T = 0.0D0
         C = Ocn_ob(j)
      Do E = 1, Nvrt_a
      Do M = 1, Nocc_a
         T = T + T2aa(A,E,K,M)*W5_ba(M,i,E,j)*C
      Enddo
      Enddo
         Hbar_iajk4(i,A,j,K) = Hbar_iajk4(i,A,j,K) + T
      Enddo
      Enddo
      Enddo
      Enddo

      Do K = 1, Nocc_a
      Do J = 1, Nocc_b
      Do A = 1, Nvrt_a
      Do I = 1, Nocc_b
         T = 0.0D0
         C = Ocn_ob(j)
      Do E = 1, Nvrt_b
      Do M = 1, Nocc_b
         T = T + T2ab(A,e,K,m)*W5_bb(i,m,j,e)*C
      Enddo
      Enddo
         Hbar_iajk4(i,A,j,K) = Hbar_iajk4(i,A,j,K) + T
      Enddo
      Enddo
      Enddo
      Enddo

      Do K = 1, Nocc_a
      Do J = 1, Nocc_b
      Do A = 1, Nvrt_a
      Do I = 1, Nocc_b
         T = 0.0D0
         C = Ocn_oa(K)
      Do E = 1, Nvrt_b
      Do M = 1, Nocc_a
         T = T + T2ab(A,e,M,j)*W5_ab(M,i,K,e)*C
      Enddo
      Enddo
         Hbar_iajk4(i,A,j,K) = Hbar_iajk4(i,A,j,K) - T
      Enddo
      Enddo
      Enddo
      Enddo



      Return
      End
