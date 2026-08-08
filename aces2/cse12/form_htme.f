










      Subroutine Form_htme(T1aa,T1bb,Fme_a,Fme_b,W_aa,W_bb,W_ab,
     +                     Fockov_a,Fockov_b,Nocc_a,Nocc_b,Nvrt_a,
     +                     Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension Fme_a(Nocc_a,Nvrt_a),Fme_b(Nocc_b,Nvrt_b)
      Dimension Fockov_a(Nocc_a,Nvrt_a),Fockov_b(Nocc_b,Nvrt_b)

      Integer M,E,N,F

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



      
C Fme(M,E) = T1(N,F)*W(MN,EF) + T1(n,f)*W(M,n,E,f) + Fockov(M,E)

      DO E = 1, Nvrt_a
      Do M = 1, Nocc_a
         T = 0.0D0
         Fme_a(M,E) = Fockov_a(M,E)
      Do N = 1, Nocc_a
      Do F = 1, Nvrt_a 
         T = T + T1aa(F,N)*W_aa(M,N,E,F)
      Enddo
      Enddo
         Fme_a(M,E) = Fme_a(M,E) + T
      Enddo
      Enddo

      DO E = 1, Nvrt_a
      Do M = 1, Nocc_a
         T = 0.0D0
      Do N = 1, Nocc_b
      Do F = 1, Nvrt_b
         T = T + T1bb(f,n)*W_ab(M,n,E,f)
      Enddo
      Enddo
         Fme_a(M,E) = Fme_a(M,E) + T
      Enddo
      Enddo

C Fme(m,e) = T1(n,f)*W(m,n,e,f) + T1(N,F)*W(m,N,e,F) + Fockov(m,e)

      DO E = 1, Nvrt_b
      Do M = 1, Nocc_b
         T = 0.0D0
         Fme_b(m,e) = Fockov_b(m,e)
      Do N = 1, Nocc_b
      Do F = 1, Nvrt_b
         T = T + T1bb(f,n)*W_bb(m,n,e,f)
      Enddo
      Enddo
         Fme_b(m,e) = Fme_b(m,e) + T
      Enddo
      Enddo

      DO E = 1, Nvrt_b
      Do M = 1, Nocc_b
         T = 0.0D0
      Do N = 1, Nocc_a
      Do F = 1, Nvrt_a
         T =  T + T1aa(F,N)*W_ab(N,m,F,e)
      Enddo
      Enddo
         Fme_b(m,e) = Fme_b(m,e) + T
      Enddo
      Enddo

      call checksum("Htme_a    :",Fme_a,Nocc_a*Nvrt_a)
      call checksum("Htme_b    :",Fme_b,Nocc_b*Nvrt_b)
       
      Return 
      End
