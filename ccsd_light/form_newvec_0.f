










      Subroutine Form_newvec_0(D1invt_aa,D1invt_bb,D2invt_aa,D2invt_bb,
     +                         D2invt_ab,Newt1_aa,Newt1_bb,Newt2_aa,
     +                         Newt2_bb,Newt2_ab,Nocc_a,Nocc_b,Nvrt_a,
     +                         Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Double Precision Newt1_aa,Newt1_bb,Newt2_aa,Newt2_bb,
     +                 Newt2_ab

c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end

      Integer cc_maxcyc
      Integer Act_min_a,Act_min_b,Act_max_a,Act_max_b
      Integer Lineq_mxcyc
      Logical Ring_cc,Brueck,Active_space,Regular
      Double Precision ocn_oa,Ocn_ob,Ocn_va,Ocn_vb
      Double Precision Denom_tol,Brueck_tol,Lineq_tol
      Double Precision Rfac
      Dimension E_corr(0:500)

      Common /ccsdlight_vars/Ring_cc,Brueck,cc_conv,cc_maxcyc,
     +                       ocn_oa(Maxbasfn),ocn_ob(Maxbasfn),
     +                       ocn_va(Maxbasfn),ocn_vb(Maxbasfn),
     +                       E_corr,Denom_tol,Brueck_tol,Lineq_tol,
     +                       Act_min_a,Act_min_b,Act_max_a,
     +                       Act_max_b,Active_space,Lineq_mxcyc,
     +                       Regular,Rfac
     +                       




      Dimension D1invt_aa(Nvrt_a,Nocc_a)
      Dimension D1invt_bb(Nvrt_b,Nocc_b)
      Dimension D2invt_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension D2invt_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension D2invt_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Newt1_aa(Nvrt_a,Nocc_a)
      Dimension Newt1_bb(Nvrt_b,Nocc_b)
      Dimension Newt2_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Newt2_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Newt2_ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

C AA block

      DO I = 1, Nocc_a
      DO A = 1, Nvrt_a
         Newt1_aa(A,I) = D1invt_aa(A,I)
      ENDDO
      ENDDO

C BB block

      DO I = 1, Nocc_b
      DO A = 1, Nvrt_b
         Newt1_bb(a,i) = D1invt_bb(a,i)
      ENDDO
      ENDDO

C AAAA block

      DO J = 1, Nocc_a
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_a
      DO A = 1, Nvrt_a
         Newt2_aa(A,B,I,J) = D2invt_aa(A,B,I,J)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C BBBB block

      DO J = 1, Nocc_b
      DO I = 1, Nocc_b
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_b
         Newt2_bb(a,b,i,j) = D2invt_bb(a,b,i,j)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C ABAB block

      DO J = 1, Nocc_b
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_a
         Newt2_ab(A,b,I,j) = D2invt_ab(A,b,I,j)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      L_aa = Nocc_a*Nvrt_a
      L_bb = Nocc_b*Nvrt_b
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("Newt1_aa  :",Newt1_aa,L_aa)
      call checksum("Newt1_ab  :",Newt1_bb,L_bb)
      call checksum("Newt2_aa  :",Newt2_aa,L_aaaa)
      call checksum("Newt2_bb  :",Newt2_bb,L_bbbb)
      call checksum("Newt2_ab  :",Newt2_ab,L_abab)
      Write(6,*) 

      Return
      End
