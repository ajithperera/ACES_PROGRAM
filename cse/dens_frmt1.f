










      Subroutine Dens_frmt1(Densp_aa,Densp_bb,Densm_aa,Densm_bb,
     +                      T1aa,T1bb,Oneh_aa,Oneh_bb,Repuls,Work,
     +                      Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,Nbasis,
     +                      Memleft,Frac_occ,E1)

      Implicit Double Precision(A-H,O-Z)
      Double Precision Oneh_aa,Onehh_bb

      Dimension Work(Memleft)
      Dimension Nocca(8),Noccb(8)
 
      Dimension Densp_aa(Nbasis,Nbasis)
      Dimension Densp_bb(Nbasis,Nbasis)
      Dimension Densm_aa(Nbasis,Nbasis)
      Dimension Densm_bb(Nbasis,Nbasis)

      Dimension Oneh_aa(Nbasis,Nbasis)
      Dimension Oneh_bb(Nbasis,Nbasis)
 
      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)
 
      Integer P,Q,A,I,J,E,M
      Logical Frac_occ

      Data Ione /1/

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




      If (Frac_occ) Then

      Do P = 1, Nbasis
      Do Q = 1, Nbasis
         If (P .EQ. Q) Then
            Densp_aa(Q,P) =  T1aa(Q,P) + Ocn_oa(P)
            Densm_aa(Q,P) = -T1aa(Q,P) + (1.0D0-Ocn_oa(P))
         Else 
            Densp_aa(Q,P) =  T1aa(Q,P)
            Densm_aa(Q,P) = -T1aa(Q,P) 
         Endif 
      Enddo
      Enddo 
      
      Do P = 1, Nbasis
      Do Q = 1, Nbasis
         If (p .EQ. q) Then
            Densp_bb(q,p) =  T1bb(q,p) + Ocn_ob(p)
            Densm_bb(Q,P) = -T1bb(Q,P) + (1.0D0-Ocn_ob(p))
         Else
            Densp_bb(q,p) =  T1bb(q,p)
            Densm_bb(q,p) = -T1bb(q,p)
         Endif
      Enddo
      Enddo
 
      Else

      Call form_canonical_dens(Densp_aa,Densp_bb,Densm_aa,Densm_bb,
     +                         T1aa,T1bb,Nocc_a,Nocc_b,Nvrt_a,
     +                         Nvrt_b,Nbasis)
      Endif 
      E1_aa = Ddot(Nbasis*Nbasis,Oneh_aa,1,Densp_aa,1)
      E1_bb = Ddot(Nbasis*Nbasis,Oneh_bb,1,Densp_bb,1)
      E1    = E1_aa+E1_bb+Repuls 
      Write(6,"(a,(1X,F15.10))") "D(p,q)*h(p,q)        :",E1

      If (Frac_occ) Then
      Write(6,*) " The alpha and beat density matrices"
      Call Output(Densp_aa,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
      Call Output(Densp_bb,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
      Write(6,*) 
      Ipdh = 1
      Imdh = Ipdh + Nbasis*Nbasis 
      Iend = Imdh + Nbasis*Nbasis 
      Call Dgemm("N","N",Nbasis,Nbasis,Nbasis,1.0D0,Densm_aa,Nbasis,
     +            Oneh_aa,Nbasis,0.0D0,Work(ipdh),Nbasis)
      Call Dgemm("N","N",Nbasis,Nbasis,Nbasis,1.0D0,Work(ipdh),Nbasis,
     +            Densp_aa,Nbasis,0.0D0,Work(imdh),Nbasis)
      Write(6,"(a)") "The test (1-D)HD" 
      Call Output(Work(imdh),1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
C Lets look at the eigenvalues and the eigenvectors.
      Ievall = 1
      Ievalr = Ievall + Nbasis 
      Ievecr = Ievalr + Nbasis 
      Ievecl = Ievecr + Nbasis*Nbasis
      Iscr1  = Ievecl + Nbasis*Nbasis 
      Iend   = Iscr1  + 4*Nbasis 
     
      Call Dgeev("V","V",Nbasis,Densp_aa,Nbasis,Work(Ievalr),
     +            Work(Ievall),Work(Ievecl),Nbasis,Work(Ievecr),
     +            Nbasis,Work(Iscr1),4*Nbasis,Ierr)
      Write(6,"(a)") "The eigenvalues of D" 
      Write(6,"(6(1x,F15.10))") (Work(Ievalr-1+j),j=1,Nbasis)
      Write(6,"(a)") "The eigenvectors of D" 
      Call Output(Work(Ievecr),1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
      Endif 
 
      Return
      End
