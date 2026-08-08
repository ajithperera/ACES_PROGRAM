










      Subroutine Symmetrize_t2(T2aa_old,T2bb_old,T2ab_old,Nocc_a,
     +                         Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)
      
      Dimension T2aa_old(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb_old(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab_old(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Integer I,J,A,B

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




C AAAA and BBBB block 

      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
       
      Call Dcopy(L_aaaa,T2ab_old,1,T2aa_old,1)
      Call Dcopy(L_bbbb,T2ab_old,1,T2bb_old,1)

      Write(6,*) 
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("T2aa      :",T2aa_old,L_aaaa)
      call checksum("T2bb      :",T2bb_old,L_bbbb)
      call checksum("T2ab      :",T2ab_old,L_abab)

      Return
      End
