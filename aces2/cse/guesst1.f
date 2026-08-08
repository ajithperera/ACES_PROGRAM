










      Subroutine Guesst1(T1aa_old,T1bb_old,Dens_aa,Dens_bb,Nocc_a,
     +                   Nocc_b,Nvrt_a,Nvrt_b,Nbasis)

      Implicit Double Precision(A-H,O-Z)
      
      Dimension T1aa_old(Nvrt_a,Nocc_a)
      Dimension T1bb_old(Nvrt_b,Nocc_b)
      Dimension Dens_aa(Nbasis,Nbasis)
      Dimension Dens_bb(Nbasis,Nbasis)

      Integer I,A

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





C#ifdef _NOSKIP
C AA block

      DO I = 1, Nocc_a
      DO A = 1, Nvrt_a
         T1aa_old(A,I) = 0.0D0
      ENDDO
      ENDDO

C BB block

      DO I = 1, Nocc_b
      DO A = 1, Nvrt_b
         T1bb_old(a,i) = 0.0D0
      ENDDO
      ENDDO
C#endif 


      Return
      End
