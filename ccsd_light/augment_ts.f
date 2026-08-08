










      Subroutine Augment_ts(T1aa,T1bb,T2aa,T2bb,T2ab,Nocc_a,Nocc_b,
     +                      Nvrt_a,Nvrt_b,Daa,Dbb)

      Implicit Double Precision(A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

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

      C = 1.0D00
      D = 1.0D00
      Do I = 1, Nocc_a
      Do A = 1, Nvrt_a
         T1aa(A,I) = C*T1aa(A,I) + Daa(A,I)*D
      Enddo
      Enddo

      Do I = 1, Nocc_b
      Do A = 1, Nvrt_b
         T1bb(a,i) = C*T1bb(a,i) + Dbb(a,i)*D
      Enddo
      Enddo

      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T2aa(A,B,I,J) = C*T2aa(A,B,I,J) + Daa(A,I)*Daa(B,J)*D
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T2bb(a,b,i,j) = C*T2bb(a,b,i,j) + Dbb(a,i)*Dbb(b,j)*D
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T2ab(A,b,I,j) = C*T2ab(A,b,I,j) + Daa(A,I)*Dbb(b,j)*D
      Enddo
      Enddo
      Enddo
      Enddo

      Return
      End
