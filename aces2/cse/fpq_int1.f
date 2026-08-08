










      Subroutine Fpq_int1(T1aa,T1bb,T1resid_aa,T1resid_bb,
     +                    Fockoo_a,Fockoo_b,Fockvv_a,Fockvv_b,
     +                    Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)
      
      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension T1resid_aa(Nvrt_a,Nocc_a)
      Dimension T1resid_bb(Nvrt_b,Nocc_b)

      Dimension Fockoo_a(Nocc_a,Nocc_a)
      Dimension Fockoo_b(Nocc_b,Nocc_b)
      Dimension Fockvv_a(Nvrt_a,Nvrt_a)
      Dimension Fockvv_b(Nvrt_b,Nvrt_b)

      Integer I,J,A,B

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




C AA block 

      DO I = 1, Nocc_a
      DO A = 1, Nvrt_a
         T = 0.0D0
         C = Ocn_oa(I)
      DO K = 1, Nocc_a
         If (K .EQ. I) Then
             T = T + T1aa(A,K)*Fockoo_a(K,I)*C
         Endif 
      ENDDO 
         T1resid_aa(A,I) = T1resid_aa(A,I) - T
      ENDDO
      ENDDO

      DO I = 1, Nocc_a
      DO A = 1, Nvrt_a
         T = 0.0D0
         C = (1.0D0-Ocn_va(A))
      DO B = 1, Nvrt_a
         If (A .EQ. B) Then
             T = T + T1aa(B,I)*Fockvv_a(B,A)*C 
         Endif 
      ENDDO 
         T1resid_aa(A,I) = T1resid_aa(A,I) + T
      ENDDO
      ENDDO

C BB block 

      DO I = 1, Nocc_b
      DO A = 1, Nvrt_b
         T = 0.0D0
         C = Ocn_ob(i)
      DO K = 1, Nocc_b
         If (k .EQ. i) Then
             T = T + T1bb(a,k)*Fockoo_b(k,i)*C
         Endif 
      ENDDO
         T1resid_bb(a,i) = T1resid_bb(a,i) - T
      ENDDO
      ENDDO

      DO I = 1, Nocc_b
      DO A = 1, Nvrt_b
         T = 0.0D0
         C = (1.0D0-Ocn_vb(a))
      DO B = 1, Nvrt_b
         If (a .EQ. b) Then
             T = T + T1aa(b,i)*Fockvv_b(b,a)*C
         Endif 
      ENDDO
         T1resid_bb(a,i) = T1resid_bb(a,i) + T
      ENDDO
      ENDDO

      L_aa = Nocc_a*Nvrt_a
      L_bb = Nocc_b*Nvrt_b
      call checksum("T1aa      :",T1resid_aa,L_aa)
      call checksum("T1bb      :",T1resid_bb,L_bb)
      do i=1,Nocc_a
      do a=1,Nvrt_a
      Write(6,"(2I2,1x,F15.10)") a,i,T1resid_aa(a,i)
      enddo
      enddo

      Return
      End
