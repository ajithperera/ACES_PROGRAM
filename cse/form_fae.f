










      Subroutine Form_fae(T1aa,T1bb,Fae_a,Fae_b,Fockvv_od_a,
     +                    Fockvv_od_b,Fockov_a,Fockov_b,Nocc_a,
     +                    Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension Fae_a(Nvrt_a,Nvrt_a)
      Dimension Fae_b(Nvrt_b,Nvrt_b)

      Dimension Fockov_a(Nocc_a,Nvrt_a)
      Dimension Fockov_b(Nocc_b,Nvrt_b)
 
      Dimension Fockvv_od_a(Nvrt_a,Nvrt_a)
      Dimension Fockvv_od_b(Nvrt_b,Nvrt_b)
   
      Integer M,N,A,F,E

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




C Fae_a(A,E)  = -1/2f(M,E)T(A,M)

      sum = 0.0D0
      Do e = 1, Nvrt_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do m = 1, Nocc_a
         T = T + 0.50D0*T1aa(A,M)*Fockov_a(M,E)
      Enddo
         C = (1.0D0-Ocn_va(A))
         Fae_a(A,E) =  Fockvv_od_a(A,E)*C - T
      Enddo
      Enddo

C Fae_a(a,e)  = -1/2f(m,e)T(a,m)

      Do e = 1, Nvrt_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do m = 1, Nocc_b
         T = T + 0.50D0*T1bb(a,m)*Fockov_b(m,e)
      Enddo
         C = (1.0D0-Ocn_vb(a))
         Fae_b(a,e) = Fockvv_od_b(a,e)*C - T
      Enddo
      Enddo

      call checksum("Fae_a     :",Fae_a,Nvrt_a*Nvrt_a)
      call checksum("Fae_b     :",Fae_b,Nvrt_b*Nvrt_b)
      Write(6,"(a)") "From F(a,e)"
      do e=1,Nvrt_a
      do a=1,Nvrt_a
      Write(6,"(2I2,1x,F15.10)") a,e,Fae_a(a,e)
      enddo
      enddo

      Return 
      End 
