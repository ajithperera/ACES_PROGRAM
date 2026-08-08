










      Subroutine Energy(T2aa_old,T2bb_old,T2ab_old,W0_aa,W0_bb,
     +                  W0_ab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,E2_aa,
     +                  E2_bb,E2_ab) 

      Implicit Double Precision(A-H,O-Z)
      
      Dimension T2aa_old(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb_old(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab_old(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension W0_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension w0_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension W0_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)
 
      Integer A,B,I,J

c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end

      Integer cc_maxcyc
      Logical Ring_cc
      Double Precision ocn_oa,Ocn_ob,Ocn_va,Ocn_vb
      Double Precision Denom_tol  
      Dimension E_corr(0:500)

      Common /ccdlight_vars/Ring_cc,cc_conv,cc_maxcyc,
     +                      ocn_oa(Maxbasfn),ocn_ob(Maxbasfn),
     +                      ocn_va(Maxbasfn),ocn_vb(Maxbasfn),
     +                      E_corr,Denom_tol



      
C AAAA block 

      
      E2_aa = 0.0D0
      DO J = 1, Nocc_a
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_a 
      DO A = 1, Nvrt_a
         C = (1.0D0-Ocn_va(A))*(1.0D0-Ocn_va(B))*Ocn_oa(I)*Ocn_oa(J)
         E2_aa = E2_aa + T2aa_old(A,B,I,J)*W0_aa(A,B,I,J)*C
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C BBBB block 

      E2_bb = 0.0D0
      DO J = 1, Nocc_b
      DO I = 1, Nocc_b
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_b
         C = (1.0D0-Ocn_vb(A))*(1.0D0-Ocn_vb(B))*Ocn_ob(I)*Ocn_ob(J) 
         E2_bb = E2_bb+ T2bb_old(a,b,i,j)*W0_bb(a,b,i,j)*C
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C ABAB block 

      E2_ab = 0.0D0
      DO J = 1, Nocc_b
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_a
         C = (1.0D0-Ocn_va(A))*(1.0D0-Ocn_vb(B))*Ocn_oa(I)*Ocn_ob(J) 
         E2_ab = E2_ab + T2ab_old(A,b,I,j)*W0_ab(A,b,I,j)*C
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      Write(6,*)
      E2_a = E2_aa*0.25
      E2_b = E2_bb*0.25
      Write(6,"(a,3(1X,F15.10))") "E2aa,E2bb,E2ab       :",E2_a,E2_b,
     +                             E2_ab
     
      Return 
      End

