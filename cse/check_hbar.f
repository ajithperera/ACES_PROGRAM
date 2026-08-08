










      Subroutine check_hbar(T1aa_old,T1bb_old,T1resid_aa,T1resid_bb,
     +                      H1bar_aa,H1bar_bb,Nocc_a,Nocc_b,Nvrt_a,
     +                      Nvrt_b,Nbasis)

      Implicit Double Precision(A-H,O-Z)
      
      Dimension T1aa_old(Nvrt_a,Nocc_a)
      Dimension T1bb_old(Nvrt_b,Nocc_b)

      Dimension T1resid_aa(Nvrt_a,Nocc_a)
      Dimension T1resid_bb(Nvrt_b,Nocc_b)

      Dimension H1bar_aa(2*Nbasis,2*Nbasis)
      Dimension H1bar_bb(2*Nbasis,2*Nbasis)

      Data Onem,Onep/-1.0D0,1.0D0/

      Integer A,B,E,I,J

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

      Do I = 1, Nocc_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do E = 1, Nvrt_a
         T = T + H1bar_aa(Nbasis+A,Nbasis+E)*T1aa_old(E,I)
      Enddo
         T1resid_aa(A,I) = T + H1bar_aa(Nocc_a+A,I)
      Enddo
      Enddo

      Do I = 1, Nocc_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
         T = T  + H1bar_aa(M,I)*T1aa_old(A,M)
      Enddo
         T1resid_aa(A,I) = T1resid_aa(A,I) - T
      Enddo
      Enddo

C BB block 

      Do I = 1, Nocc_b
      Do A = 1, Nvrt_b
         T = 0.0D0
      Do E = 1, Nvrt_b
         T = T + H1bar_bb(Nbasis+a,Nbasis+e)*T1bb_old(e,i)
      Enddo
         T1resid_bb(a,i) = T + H1bar_bb(Nocc_b+a,i)
      Enddo
      Enddo

      Do I = 1, Nocc_b
      Do A = 1, Nvrt_b
         T = 0.0D0
      Do M = 1, Nocc_b
         T = T + H1bar_bb(m,i)*T1bb_old(a,m)
      Enddo
         T1resid_bb(a,i) = T1resid_bb(a,i) - T
      Enddo
      Enddo

      Write(6,"(a)") "Printing from check_hbar"
      Ldm = 2*Nbasis 
      Call output(H1bar_aa,1,Ldm,1,Ldm,Ldm,Ldm,1)
      L_aa = Nocc_a*NvRT_a
      L_bb = Nocc_b*Nvrt_b
      call checksum("T1aa      :",T1resid_aa,L_aa)
      call checksum("T1bb      :",T1resid_bb,L_bb)

      Return
      End
