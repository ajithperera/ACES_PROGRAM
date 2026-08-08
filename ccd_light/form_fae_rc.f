










      Subroutine Form_fae_rc(T2aa,T2bb,T2ab,Fae_a,Fae_b,W_aa,W_bb, 
     +                       W_ab,Fockvv_od_a,Fockvv_od_b,Nocc_a,
     +                       Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a) 
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b) 
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b) 

      Dimension Fae_a(Nvrt_a,Nvrt_a)
      Dimension Fae_b(Nvrt_b,Nvrt_b)
 
      Dimension Fockvv_od_a(Nvrt_a,Nvrt_a)
      Dimension Fockvv_od_b(Nvrt_b,Nvrt_b)
   
      Integer M,N,A,F,E

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




C Fae_aa(A,E)= -1/2 T2(AF,MN)*W(MN,EF) - T2(Af,Mn)*W(Mn,Ef)
      
      Fac = 1.0D0
      Do e = 1, Nvrt_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do m = 1, Nocc_a
      Do n = 1, Nocc_a
      Do f = 1, Nvrt_a
         T = T - T2aa(A,F,N,M)*W_aa(N,M,E,F)*Fac
      Enddo
      Enddo
      Enddo 
      Fae_a(A,E) = T + Fockvv_od_a(a,e)
      Enddo 
      Enddo 

      Do e = 1, Nvrt_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do m = 1, Nocc_b
      Do n = 1, Nocc_a
      Do f = 1, Nvrt_b
         T = T - T2ab(A,f,N,m)*W_ab(N,m,E,f)
      Enddo
      Enddo
      Enddo 
         Fae_a(A,E) = Fae_a(A,E) + T
      Enddo 
      Enddo 

C Fae_bb(a,e)= -1/2 T2(af,mn)*W2(mn,ef) - T2(mN,aF)*W2(mN,eF)

      Do e = 1, Nvrt_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do m = 1, Nocc_b
      Do n = 1, Nocc_b
      Do f = 1, Nvrt_b
         T = T - T2bb(a,f,n,m)*W_bb(n,m,e,f)*Fac 
      Enddo
      Enddo
      Enddo
         Fae_b(a,e)  = T + Fockvv_od_b(a,e)
      Enddo
      Enddo 

      Do e = 1, Nvrt_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do m = 1, Nocc_a
      Do n = 1, Nocc_b
      Do f = 1, Nvrt_a
         T = T - T2ab(F,a,M,n)*W_ab(M,n,F,e)
      Enddo
      Enddo
      Enddo
         Fae_b(a,e)  = Fae_b(a,e) + T
      Enddo
      Enddo 

      call checksum("Fae_a      :",Fae_a,Nvrt_a*Nvrt_a)
      call checksum("Fae_b      :",Fae_b,Nvrt_b*Nvrt_b)

      Return 
      End 
