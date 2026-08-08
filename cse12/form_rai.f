










       Subroutine Form_rai(T1aa,T1bb,R1ai_aa,R1ai_bb,Fockoo_a,Fockoo_b,
     +                     Fockvv_a,Fockvv_b,Nocc_a,Nocc_b,Nvrt_a,
     +                     Nvrt_b,Scale)

      Implicit Double Precision (A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension Fockoo_a(Nocc_a,Nocc_a)
      Dimension Fockoo_b(Nocc_b,Nocc_b)
      Dimension Fockvv_a(Nvrt_a,Nvrt_a)
      Dimension Fockvv_b(Nvrt_b,Nvrt_b)

      Dimension R1ai_aa(Nvrt_a,Nocc_a)
      Dimension R1ai_bb(Nvrt_b,Nocc_b)

      Integer A,B,I,K

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




C Rai_a(A,I) <- -T1(A,K)*Fockoo_a(K,I)

      DO I = 1, Nocc_a
      DO A = 1, Nvrt_a
         T = 0.0D0
         C = Ocn_oa(I)
      DO K = 1, Nocc_a
         If (K .EQ. I) Then
             T = T + T1aa(A,K)*Fockoo_a(K,I)*C
         Endif
      ENDDO
         R1ai_aa(A,I) = R1ai_aa(A,I) - T*Scale 
      ENDDO
      ENDDO

C Rai_a(A,I) <- +T1(B,I)*Fockvv_a(B,A)

      DO I = 1, Nocc_a
      DO A = 1, Nvrt_a
         T = 0.0D0
         C = (1.0D0-Ocn_va(A))
      DO B = 1, Nvrt_a
         If (A .EQ. B) Then
             T = T + T1aa(B,I)*Fockvv_a(B,A)*C
         Endif
      ENDDO
         R1ai_aa(A,I) = R1ai_aa(A,I) + T*Scale 
      ENDDO
      ENDDO

C Rai_b(a,i) <- -T1(a,k)*Fockoo_b(k,i)

      DO I = 1, Nocc_b
      DO A = 1, Nvrt_b
         T = 0.0D0
         C = Ocn_ob(i)
      DO K = 1, Nocc_b
         If (k .EQ. i) Then
             T = T + T1bb(a,k)*Fockoo_b(k,i)*C
         Endif
      ENDDO
         R1ai_bb(a,i) = R1ai_bb(a,i) - T*Scale 
      ENDDO
      ENDDO

C Rai_b(a,i) <- +T1(b,i)*Fockvv_a(b,a)

      DO I = 1, Nocc_b
      DO A = 1, Nvrt_b
         T = 0.0D0
         C = (1.0D0-Ocn_vb(a))
      DO B = 1, Nvrt_b
         If (a .EQ. b) Then
             T = T + T1aa(b,i)*Fockvv_b(b,a)*C
         Endif
      ENDDO
         R1ai_bb(a,i) = R1ai_bb(a,i) + T*Scale 
      ENDDO
      ENDDO

      Return 
      End

