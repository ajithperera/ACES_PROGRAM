










      Subroutine Active_sp_proj(T1aa,T1bb,T2aa,T2bb,T2ab,P,Length,
     +                          Nocc_a,Nvrt_a,Nocc_b,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)
      Dimension P(Length) 

      Integer I,J,A,B,Inill
      Double Precision Dnill 
      Data Ione,One,INill,DNill/1,1.0D0,0,0.0D0/

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




      Index  = Inill 
      Ioff   = Inill 

      Sum = 0.0D0
      Do i=1,Length
         Sum = Sum + P(Ioff+i)*P(Ioff+i)
      Enddo
      Write(6,"(a,1x,2(1x,ES8.2E2))") "@-act_sp_proj,P*P    : ",Sum
      L_aa = Nocc_a*Nvrt_a
      L_bb = Nocc_b*Nvrt_b
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("T1aa      :",T1aa,L_aa)
      call checksum("T1bb      :",T1bb,L_bb)
      call checksum("T2aa      :",T2aa,L_aaaa)
      call checksum("T2bb      :",T2bb,L_bbbb)
      call checksum("T2ab      :",T2ab,L_abab)

C AA block 

      DO I = 1, Nocc_a
      DO A = 1, Nvrt_a
         Index = Index + Ione
         T1aa(A,I) = T1aa(A,I)*P(ioff+Index) 
      ENDDO
      ENDDO

C BB block

      Index = Index - Ione 
      DO i = 1, Nocc_b
      DO a = 1, Nvrt_b
         Index = Index + Ione
         T1bb(a,i) = T1bb(a,i)*P(ioff+Index)
      ENDDO
      ENDDO

C AAAA block

      Index = Index - Ione 
      DO J = 1, Nocc_a
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_a
      DO A = 1, Nvrt_a
         Index = Index + Ione
         T2aa(A,B,I,J) = T2aa(A,B,I,J)*P(ioff+Index)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C BBBB block

      Index = Index - Ione 
      DO J = 1, Nocc_b
      DO I = 1, Nocc_b
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_b
         Index = Index + Ione
         T2bb(a,b,i,j) = T2bb(a,b,i,j)*P(ioff+Index)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C ABAB block

      Index = Index - Ione 
      DO J = 1, Nocc_b
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_a
         Index = Index + Ione
         T2ab(A,b,I,j) = T2ab(A,b,I,j)*P(ioff+Index) 
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      call checksum("T1aa      :",T1aa,L_aa)
      call checksum("T1bb      :",T1bb,L_bb)
      call checksum("T2aa      :",T2aa,L_aaaa)
      call checksum("T2bb      :",T2bb,L_bbbb)
      call checksum("T2ab      :",T2ab,L_abab)
      Return
      End
