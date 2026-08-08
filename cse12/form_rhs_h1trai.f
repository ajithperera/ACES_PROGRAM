










      Subroutine Form_rhs_h1trai(Htem_aa,Htem_bb,Htme_aa,Htme_bb,
     +                           T1aa,T1bb,T2aa,T2bb,T2ab,W_aa,W_bb,
     +                           W_ab,W4_aa,W4_bb,W4_ab,W4_ba,W5_aa,
     +                           W5_bb,W5_ab,W5_ba,T1resid_aa,
     +                           T1resid_bb,Fockoo_a,Fockoo_b,Fockvv_a,
     +                           Fockvv_b,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b) 

      Implicit Double Precision(A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)
 
      Dimension T1resid_aa(Nvrt_a,Nocc_a)
      Dimension T1resid_bb(Nvrt_b,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension W4_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W4_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W4_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W4_ba(Nvrt_a,Nvrt_b,Nocc_a,Nvrt_b)

      Dimension W5_aa(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension W5_bb(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension W5_ab(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension W5_ba(Nocc_a,Nocc_b,Nvrt_a,Nocc_b)
     
      Dimension Fockoo_a(Nocc_a,Nocc_a)
      Dimension Fockoo_b(Nocc_b,Nocc_b)
      Dimension Fockvv_a(Nvrt_a,Nvrt_a)
      Dimension Fockvv_b(Nvrt_b,Nvrt_b)

      Dimension Htme_aa(Nocc_a,Nvrt_a)
      Dimension Htme_bb(Nocc_b,Nvrt_b)
      Dimension Htem_aa(Nvrt_a,Nocc_a)
      Dimension Htem_bb(Nvrt_b,Nocc_b)

      Integer A,B,I,J,M

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
C Rai_a(A,I) <- Htme_aa(M,E)*T2aa(A,E,I,M) +Htme_bb(m,e)*T(A,e,I,m)
C             + Htem_aa(A,I)

      Do I = 1, Nocc_a
      Do A = 1, Nvrt_a
         T = 0.0D0
         T1resid_aa(A,I)  = Htem_aa(A,I)
      Do E = 1, Nvrt_a 
      Do M = 1, Nocc_a 
         T = T  + Htme_aa(M,E)*T2aa(A,E,I,M)
      Enddo 
      Enddo 
         T1resid_aa(A,I) = T1resid_aa(A,I) + T 
      Enddo 
      Enddo 

      Do I = 1, Nocc_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do E = 1, Nvrt_b
      Do M = 1, Nocc_b
         T = T  + Htme_bb(m,e)*T2ab(A,e,I,m)
      Enddo
      Enddo
         T1resid_aa(A,I) = T1resid_aa(A,I) + T
      Enddo
      Enddo

C BB block  
C Rai_b(a,i) <- Htme_bb(m,e)*T2aa(a,e,i,m) +Htme_bb(M,E)*T2ab(E,a,M,i)
C             + Htem_bb(a,i)
      Do I = 1, Nocc_b
      Do A = 1, Nvrt_b
         T = 0.0D0
         T1resid_bb(a,i) = Htem_bb(a,i)
      Do E = 1, Nvrt_b
      Do M = 1, Nocc_b
         T = T  + Htme_bb(m,e)*T2bb(a,e,i,m)
      Enddo
      Enddo
         T1resid_bb(a,i) = T1resid_bb(a,i) + T
      Enddo
      Enddo

      Do I = 1, Nocc_b
      Do A = 1, Nvrt_b
         T = 0.0D0
      Do E = 1, Nvrt_a 
      Do M = 1, Nocc_a
         T = T  + Htme_aa(M,E)*T2ab(E,a,M,i)
      Enddo
      Enddo
         T1resid_bb(a,i) = T1resid_bb(a,i) + T
      Enddo
      Enddo

      L_aa = Nocc_a*Nvrt_a
      L_bb = Nocc_b*Nvrt_b
      call checksum("T1resid_aa:",T1resid_aa,L_aa)
      call checksum("T1resid_bb:",T1resid_bb,L_bb)

      Call Form_htabci(T1aa,T1bb,W4_aa,W4_bb,W4_ab,W4_ba,W_aa,
     +                 W_bb,W_ab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Call Form_htijka(T1aa,T1bb,W5_aa,W5_bb,W5_ab,W5_ba,W_aa,
     +                 W_bb,W_ab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C Two contributions -1/2T2(EF,IM)<MA||EF> - 1/2T2(AE,MN)<MN||EI>

      Call T2_inrai(T2aa,T2bb,T2ab,W4_aa,W4_bb,W4_ab,W4_ba,W5_aa,
     +              W5_bb,W5_ab,W5_ba,T1resid_aa,T1resid_bb,Nocc_a,
     +              Nocc_b,Nvrt_a,Nvrt_b)

C Undo the addition of T1 temrs to Ht(abci) and Ht(abef) since
C we need the bare integrals for subsequent use.

      Call Modf_htabci(T1aa,T1bb,W4_aa,W4_bb,W4_ab,W4_ba,W_aa,
     +                 W_bb,W_ab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Call Modf_htijka(T1aa,T1bb,W5_aa,W5_bb,W5_ab,W5_ba,W_aa,
     +                 W_bb,W_ab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)


      L_aa = Nocc_a*Nvrt_a
      L_bb = Nocc_b*Nvrt_b
      call checksum("T1resid_aa:",T1resid_aa,L_aa)
      call checksum("T1resid_bb:",T1resid_bb,L_bb)
      Write(6,*)
      Return
      End 
