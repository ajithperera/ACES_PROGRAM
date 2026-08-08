










      Subroutine Energy_frmt1(T1aa_old,T1bb_old,W0_aa,W0_bb,W0_ab,
     +                        Fockov_a,Fockov_b,Nocc_a,Nocc_b,Nvrt_a,
     +                        Nvrt_b,E1_aa,E1_bb,E1_ab) 

      Implicit Double Precision(A-H,O-Z)
      
      Dimension T1aa_old(Nvrt_a,Nocc_a)
      Dimension T1bb_old(Nvrt_b,Nocc_b)

      Dimension W0_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension w0_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension W0_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Fockov_a(Nocc_a,Nvrt_a)
      Dimension Fockov_b(Nocc_b,Nvrt_b)
 
      Integer A,B,I,J

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

      E1_aa = 0.0D0
      DO J = 1, Nocc_a
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_a 
      DO A = 1, Nvrt_a
         TC = T1aa_old(A,I)*T1aa_old(B,J)
         TE = T1aa_old(B,I)*T1aa_old(A,J) 
         E1_aa = E1_aa + (TC-TE)*W0_aa(A,B,I,J)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C BB block 

      E1_bb = 0.0D0
      DO J = 1, Nocc_b
      DO I = 1, Nocc_b
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_b
         TC = T1bb_old(a,i)*T1bb_old(b,j)
         TE = T1bb_old(b,i)*T1bb_old(a,j)
         E1_bb = E1_bb + (TC-TE)*W0_bb(a,b,i,j)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C AB block

      E1_ab = 0.0D0
      DO J = 1, Nocc_b
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_a
         E1_ab = E1_ab + T1aa_old(A,I)*T1bb_old(b,j)*W0_ab(A,b,I,j)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

  
      Enhf_aa = 0.0D0
      DO I = 1, Nocc_a
      DO A = 1, Nvrt_a
         Enhf_aa = Enhf_aa + T1aa_old(A,I)*Fockov_a(A,I)
      ENDDo
      ENDDO 
     
      Enhf_bb = 0.0D0
      DO I = 1, Nocc_b
      DO A = 1, Nvrt_b
         Enhf_bb = Enhf_bb + T1aa_old(A,I)*Fockov_a(A,I)
      ENDDo
      ENDDO 

      Write(6,"(a,2(1X,F15.10))") "Enhf_aa,Enhf_bb      :",Enhf_aa,
     +                                                     Enhf_bb
      E1_aa = Enhf_aa
      E1_bb = Enhf_bb
      E1_ab = 0.0D0 

      Return 
      End

