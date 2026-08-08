










      Subroutine Fae_int1(T1aa,T1bb,T1resid_aa,T1resid_bb,Fae_a,
     +                    Fae_b,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)
 
      Dimension T1resid_aa(Nvrt_a,Nocc_a)
      Dimension T1resid_bb(Nvrt_b,Nocc_b)

      Dimension Fae_a(Nvrt_a,Nvrt_a)
      Dimension Fae_b(Nvrt_b,Nvrt_b)

      Integer I,J,A,B,E

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




C AA block 

      Do I = 1, Nocc_a
      Do A = 1, Nvrt_a
         T = 0.0D0 
      Do E = 1, Nvrt_a 
CSSS         C = (1.0D0-Ocn_va(E))
         C = 1.0D0
         T = T + Fae_a(A,E)*T1aa(E,I)*C 
      Enddo 
         T1resid_aa(A,I) = T1resid_aa(A,I) + T
      Enddo 
      Enddo 

C BB block 

      Do I = 1, Nocc_b
      Do A = 1, Nvrt_b
         T = 0.0D0 
      Do E = 1, Nvrt_b
CSS         C = (1.0D0-Ocn_vb(e))
         C = 1.0D0
         T = T + Fae_b(a,e)*T1bb(e,i)*C
      Enddo
         T1resid_bb(a,i) = T1resid_bb(a,i) + T 
      Enddo
      Enddo

      L_aa = Nocc_a*Nvrt_a
      L_bb = Nocc_b*Nvrt_b
      call checksum("T1resid_aa:",T1resid_aa,L_aa)
      call checksum("T1resid_bb:",T1resid_bb,L_bb)
C      do i=1,Nocc_a
C      do a=1,Nvrt_a
C      Write(6,"(2I2,1x,F15.10)") a,i,T1resid_aa(a,i)
C      enddo
C      enddo

      Return
      end 
