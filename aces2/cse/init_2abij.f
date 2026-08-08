










      Subroutine Init_2abij(T1resid_aa,T1resid_bb,Fockov_aa,
     +                      Fockov_bb,Nocc_a,Nocc_b,Nvrt_a,
     +                      Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension T1resid_aa(Nvrt_a,Nocc_a)
      Dimension T1resid_bb(Nvrt_b,Nocc_b)
      Dimension Fockov_aa(Nocc_a,Nvrt_a) 
      Dimension Fockov_bb(Nocc_b,Nvrt_b) 

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
         C = Ocn_oa(I)*(1.0D0-Ocn_va(A))
         T1resid_aa(A,I) = fockov_aa(I,A)*C
      ENDDO
      ENDDO

C BB block

      DO I = 1, Nocc_b
      DO A = 1, Nvrt_b
         C = Ocn_ob(i)*(1.0D0-Ocn_vb(a))
         T1resid_bb(a,i) = fockov_bb(i,a)*C 
      ENDDO
      ENDDO

      Write(6,*) 
      L_aa   = Nocc_a*Nvrt_a
      L_bb   = Nocc_b*Nvrt_b
      call checksum("T1Resid_aa:",T1Resid_aa,L_aa)
      call checksum("T1Resid_bb:",T1Resid_bb,L_bb)
      do i=1,Nocc_a
      do a=1,Nvrt_a
      Write(6,"(2I2,1x,F15.10)") a,i,T1resid_aa(a,i)
      enddo
      enddo

      Return
      End
