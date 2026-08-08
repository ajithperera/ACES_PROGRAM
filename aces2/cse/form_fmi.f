










      subroutine Form_fmi(T1aa,T1bb,Fmi_a,Fmi_b,Fockoo_od_a,
     +                    Fockoo_od_b,Fockov_a,Fockov_b,Nocc_a,
     +                    Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)
     
      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension Fockoo_od_a(Nocc_a,Nocc_a)
      Dimension Fockoo_od_b(Nocc_b,Nocc_b)

      Dimension Fmi_a(Nocc_a,Nocc_a)
      Dimension Fmi_b(Nocc_b,Nocc_b)

      Dimension Fockov_a(Nocc_a,Nvrt_a)
      Dimension Fockov_b(Nocc_b,Nvrt_b)

      Integer M,N,E,F,I

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




C  Fmi_a(M,I) =  +1/2f(M,E)*T1(E,I)

      Do i = 1, Nocc_a
      Do m = 1, Nocc_a
         T = 0.0D0
      Do e = 1, Nvrt_a
         T = T + 0.50D0*T1aa(E,I)*Fockov_a(M,E)
      Enddo
         C = Ocn_oa(I) 
         Fmi_a(m,i) = Fockoo_od_a(M,I)*C + T
      Enddo
      Enddo

C  Fmi_b(m,i) =  +1/2f(m,e)*T1(e,i)

      Do i = 1, Nocc_b
      Do m = 1, Nocc_b
         T = 0.0D0
      Do e = 1, Nvrt_b
         T = T + 0.50D0*T1bb(e,i)*Fockov_b(m,e)
      Enddo
         C = Ocn_ob(i)
         Fmi_b(m,i) = Fockoo_od_b(m,i)*C  + T
      Enddo
      Enddo

      call checksum("Fmi_a     :",Fmi_a,Nocc_a*Nocc_a)
      call checksum("Fmi_b     :",Fmi_b,Nocc_b*Nocc_b)
      do m=1,Nocc_a
      do i=1,Nocc_a
      Write(6,"(2I2,1x,F15.10)") m,i,Fmi_a(m,i)
      enddo
      enddo

      Return 
      End 
