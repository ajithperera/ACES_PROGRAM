










      Subroutine Form_fme(Fme_a,Fme_b,Fockov_a,Fockov_b,Nocc_a,
     +                    Nocc_b,Nvrt_a,Nvrt_b,Fme_on)

      Implicit Double Precision(A-H,O-Z)

      Dimension Fme_a(Nocc_a,Nvrt_a),Fme_b(Nocc_b,Nvrt_b)
      Dimension Fockov_a(Nocc_a,Nvrt_a),Fockov_b(Nocc_b,Nvrt_b)

      Integer M,E,N,F
      Logical Fme_on

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



      
      Suma = 0.0D0
      DO E = 1, Nvrt_a
      Do M = 1, Nocc_a
          Fme_a(M,E) = Fockov_a(M,E)
          Suma = Suma + Fme_a(M,E)
      Enddo
      Enddo

      Sumb = 0.0D0
      DO E = 1, Nvrt_b
      Do M = 1, Nocc_b
          Fme_b(m,e) = Fockov_b(m,e)
          Sumb = Sumb + Fme_b(m,e)
      Enddo
      Enddo

      Fme_on = Suma .Gt. 0.0D0 .OR. Sumb .Gt. 0.0D0 

      call checksum("Fme_a     :",Fme_a,Nocc_a*Nvrt_a)
      call checksum("Fme_b     :",Fme_b,Nocc_b*Nvrt_b)
       
      Return 
      End
