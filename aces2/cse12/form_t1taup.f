










      Subroutine Form_t1taup(T1taup_aa,T1taup_bb,T1taup_ab,T1aa,
     +                       T1bb,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)
      
      Dimension T1taup_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T1taup_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T1taup_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

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




C AAAA block 

      DO J = 1, Nocc_a
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_a
      DO A = 1, Nvrt_a
         T1taup_aa(A,B,I,J) = (T1aa(A,I)*T1aa(B,J)  -
     +                       T1aa(B,I)*T1aa(A,J))*0.50D0
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C BBBB block 

      DO J = 1, Nocc_b
      DO I = 1, Nocc_b
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_b
         T1taup_bb(a,b,i,j) = (T1bb(a,i)*T1bb(b,j) - 
     +                        T1bb(b,i)*T1bb(a,j))*0.50D0
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C ABAB block 

      DO J = 1, Nocc_b
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_a
         T1taup_ab(A,b,I,J) = T1aa(A,I)*T1bb(b,j)*0.50D0
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("T1taup_aa :",T1taup_aa,L_aaaa)
      call checksum("T1taup_bb :",T1taup_bb,L_bbbb)
      call checksum("T1taup_ab :",T1taup_ab,L_abab)
      Write(6,*) 

      Return
      End
