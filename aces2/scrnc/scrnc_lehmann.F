      Subroutine Scrnc_lehmann(Work,Maxcor,Iuhf)
C
C Sum-over-states (Lehmann) correlation screened interaction and GW
C self-energy from the full EOM-EE-CCSD spectrum, replacing the
C imaginary-frequency linear-response W + 4-node quadrature + Pade
C continuation route (SCRNC_DEBUG_LOG.md 2026-10-15: that route's
C real-axis Im Sigma_c is a Pade artifact, so no satellite/spectral
C observable can be built on it).
C
C Input: eom_tdm.dat from xvee EE_ROOT_FIND=FULL_DIAG (every root's
C full MO left/right transition densities, gamma^L_s and gamma^R_s,
C and excitation energy Omega_s); lehmann_input (see below). Runs in
C place of the response pipeline when lehmann_input exists.
C
C With V^R_s(pq) = Sum_rs (pq|rs) gamma^R_s(rs) (spin-summed density,
C so triplets drop out) and V^L_s likewise, the EOM-CC density response
C chi_AB(w) = Sum_s [<0|(1+Lam)A|R_s><L_s|B|0>/(w-Omega_s)
C                  - <0|(1+Lam)B|R_s><L_s|A|0>/(w+Omega_s)]
C gives W_c = v chi v as
C   W_c(pq,rt;w) = Sum_s [V^R_s(pq)V^L_s(rt)/(w-Omega_s)
C                       - V^L_s(pq)V^R_s(rt)/(w+Omega_s)].
C The pipeline's Sigma_c(p,q;iw) = -(1/2pi) Int dw' Sum_m
C W_c(pm,qm;iw')/(iw+iw'-eps_m), closed in the upper half plane, is
C then exact and analytic (no quadrature):
C   Sigma_c(p,q;z) = Sum_s [ Sum_{m occ} V^R_s(pm)V^L_s(qm)/
C                                        (z-eps_m+Omega_s)
C                          + Sum_{m vrt} V^L_s(pm)V^R_s(qm)/
C                                        (z-eps_m-Omega_s) ]
C (occupied term from the w'=-w-i*eps_m pole, eps_m<0; virtual term
C from the w'=i*Omega_s pole; the other combinations cancel).
C
C Validation outputs, comparing against the pipeline's own files:
C   wcreduc_lehmann.dat  - Re W_c reduced (p,q,m) at each requested
C                          w' in exactly wcreduc.dat's format, one
C                          block per w', plus max|Im W_c| per block
C   sigma_lehmann_imag.dat - Sigma_c(p,q;iw) at each requested w:
C                          exact (analytic) and via the pipeline's
C                          own N-node quadrature with Re W_c (Lehmann)
C
C lehmann_input (list-directed):
C   Nwp, w'(1..Nwp)          w' values for W_c dumps
C   Nquad, Omega_c           pipeline quadrature (sigma_input)
C   Next, w(1..Next)         imaginary-axis points for Sigma_c
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Work(Maxcor)
#include "machsp.com"
#include "info.com"
C
      Parameter (Maxw=100)
      Dimension Wprime(Maxw),Wext(Maxw),Xnode(Maxw),Wnode(Maxw)
      Character*8 Labelc,Labele
      Data Labelc /"SCFEVECA"/
      Data Labele /"SCFEVALA"/
C
      Open(Unit=84,File="lehmann_input",Form="FORMATTED",
     +     Status="OLD")
      Read(84,*) Nwp,(Wprime(I),I=1,Nwp)
      Read(84,*) Nquad,Omega_c
      Read(84,*) Next,(Wext(I),I=1,Next)
      Close(Unit=84)
C
      Call Getrec(20,"JOBARC","NBASTOT ",1,Nao)
      Nmo       = Nocco(1) + Nvrto(1)
      Nocc      = Nocco(1)
      Nao_pairs = Nao*Nao
      Nmo_pairs = Nmo*Nmo
      Write(6,"(a,3I6)") " @-Scrnc_lehmann: Nao, Nmo, Nocc = ",
     +                   Nao,Nmo,Nocc
C
C First pass over eom_tdm.dat: count the roots.
C
      Open(Unit=83,File="eom_tdm.dat",Form="FORMATTED",Status="OLD")
      Read(83,*) Nmo_t,Iuhf_t
      If (Nmo_t .NE. Nmo) Then
         Write(6,"(a,2I5)") " @-Scrnc_lehmann: Nmo mismatch ",
     +                      Nmo_t,Nmo
         Call Errex
      Endif
      Nrec = 0
      Do While (.True.)
         Read(83,*,End=10) Iside,Ispin,Root,R0
         Read(83,*) (Dum,I=1,Nmo_pairs)
         Nrec = Nrec + 1
      Enddo
   10 Continue
      Nrecroot = 2*(1+Iuhf_t)
      Nroot    = Nrec/Nrecroot
      If (Nroot*Nrecroot .NE. Nrec) Then
         Write(6,"(a,I8)") " @-Scrnc_lehmann: incomplete root set ",
     +                     Nrec
         Call Errex
      Endif
      Write(6,"(a,I6,a,I3)") " @-Scrnc_lehmann: roots read ",Nroot,
     +                       "  IUHF of EOM run ",Iuhf_t
C
      I_c    = 1
      I_eps  = I_c    + Nao*Nmo
      I_vao  = I_eps  + Nmo
      I_g2   = I_vao  + Nao_pairs*Nao_pairs
      I_vmo  = I_g2   + Nao_pairs*Nao_pairs
      I_vr   = I_vmo  + Nmo_pairs*Nmo_pairs
      I_vl   = I_vr   + Nmo_pairs*Nroot
      I_om   = I_vl   + Nmo_pairs*Nroot
      I_gr   = I_om   + Nroot
      I_gl   = I_gr   + Nmo_pairs
      I_gt   = I_gl   + Nmo_pairs
      I_wre  = I_gt   + Nmo_pairs
      I_wim  = I_wre  + Nmo_pairs*Nmo_pairs
      I_vs   = I_wim  + Nmo_pairs*Nmo_pairs
      I_wred = I_vs   + Nmo_pairs*Nroot
      I_wq   = I_wred + Nmo*Nmo*Nmo
      I_sx   = I_wq   + Nmo*Nmo*Nmo*Nquad
      I_sq   = I_sx   + 2*Nmo*Nmo*Next
      I_scr  = I_sq   + 2*Nmo*Nmo*Next
      Iend   = I_scr  + Nmo_pairs*Nao_pairs + Nao*Nao + Nao*Nmo +
     +                  Nmo*Nmo
      If (Iend .GT. Maxcor) Call Insmem("Scrnc_lehmann",Iend,Maxcor)
C
      Call Getrec(20,"JOBARC",Labelc,Iintfp*Nao*Nmo,Work(I_c))
      Call Getrec(20,"JOBARC",Labele,Iintfp*Nmo,Work(I_eps))
C
C Bare Coulomb (pq|rs) in the MO-pair basis, through the pipeline's
C own (already verified) AO build and AO->MO transform.
C
      Call Scrnc_sigma_bare_ao(Work(I_vao),Work(I_g2),Nao,Nao_pairs)
      Call Scrnc_wc_ao2mo(Work(I_scr),Maxcor-I_scr+1,
     +                     Work(I_vao),Nao,Nao_pairs,
     +                     Work(I_c),Nmo,Work(I_vmo),Nmo_pairs)
C
C Second pass: spin-summed gamma^R_s, gamma^L_s -> V^R_s, V^L_s.
C
      Rewind(83)
      Read(83,*) Nmo_t,Iuhf_t
      Do Iroot = 1, Nroot
         Call Zero(Work(I_gr),Nmo_pairs)
         Call Zero(Work(I_gl),Nmo_pairs)
         Do Irec = 1, Nrecroot
            Read(83,*) Iside,Ispin,Root,R0
            Read(83,*) (Work(I_gt-1+I),I=1,Nmo_pairs)
            If (Iside .EQ. 1) Then
               Call Daxpy(Nmo_pairs,1.0D0,Work(I_gt),1,Work(I_gr),1)
            Else
               Call Daxpy(Nmo_pairs,1.0D0,Work(I_gt),1,Work(I_gl),1)
            Endif
         Enddo
C        Closed-shell RHF-based EOM stores one spin case only.
         If (Iuhf_t .EQ. 0) Then
            Call Dscal(Nmo_pairs,2.0D0,Work(I_gr),1)
            Call Dscal(Nmo_pairs,2.0D0,Work(I_gl),1)
         Endif
         Work(I_om-1+Iroot) = Root
         Ioff = (Iroot-1)*Nmo_pairs
         Call Dgemv("N",Nmo_pairs,Nmo_pairs,1.0D0,Work(I_vmo),
     +              Nmo_pairs,Work(I_gr),1,0.0D0,Work(I_vr+Ioff),1)
         Call Dgemv("N",Nmo_pairs,Nmo_pairs,1.0D0,Work(I_vmo),
     +              Nmo_pairs,Work(I_gl),1,0.0D0,Work(I_vl+Ioff),1)
      Enddo
      Close(Unit=83)
C
C W_c dumps at the requested w'.
C
      Open(Unit=87,File="wcreduc_lehmann.dat",Form="FORMATTED",
     +     Status="UNKNOWN")
      Do K = 1, Nwp
         Call Scrnc_lehmann_w(Nmo_pairs,Nroot,Work(I_vr),Work(I_vl),
     +                        Work(I_om),Wprime(K),Work(I_vs),
     +                        Work(I_wre),Work(I_wim),Wimmax)
         Call Scrnc_sigma_reduce(Work(I_wre),Nmo,Nmo_pairs,
     +                           Work(I_wred))
         Write(87,"(a,1PD22.14,a,1PD12.4)") " omega' = ",Wprime(K),
     +         "   max|Im W_c| = ",Wimmax
         Write(87,*) Nmo
         Do I = 1, Nmo*Nmo*Nmo
            Write(87,*) Work(I_wred-1+I)
         Enddo
         Do I = 1, Nmo
            Write(87,*) Work(I_eps-1+I)
         Enddo
      Enddo
      Close(Unit=87)
C
C Re W_c (Lehmann) at the pipeline's quadrature nodes.
C
      Call Scrnc_freq_grid(Nquad,Omega_c,Xnode,Wnode)
      Do K = 1, Nquad
         Call Scrnc_lehmann_w(Nmo_pairs,Nroot,Work(I_vr),Work(I_vl),
     +                        Work(I_om),Xnode(K),Work(I_vs),
     +                        Work(I_wre),Work(I_wim),Wimmax)
         Call Scrnc_sigma_reduce(Work(I_wre),Nmo,Nmo_pairs,
     +                           Work(I_wq+(K-1)*Nmo*Nmo*Nmo))
         Write(6,"(a,1PD12.4,a,1PD12.4)")
     +      " @-Scrnc_lehmann: quadrature node w' = ",Xnode(K),
     +      "  max|Im W_c| = ",Wimmax
      Enddo
C
      Call Scrnc_lehmann_sigma(Nmo,Nocc,Nmo_pairs,Nroot,Next,Nquad,
     +                         Work(I_vr),Work(I_vl),Work(I_om),
     +                         Work(I_eps),Wext,Xnode,Wnode,
     +                         Work(I_wq),Work(I_sx),Work(I_sq))
C
      Open(Unit=88,File="sigma_lehmann_imag.dat",Form="FORMATTED",
     +     Status="UNKNOWN")
      Write(88,"(a,I4,a,I4,a,I6)") "Nmo=",Nmo," Next=",Next,
     +                             " Nroot=",Nroot
      Write(88,"(a)") "  j   p   q      Re(Sigma_exact)"//
     +   "      Im(Sigma_exact)       Re(Sigma_quad)"//
     +   "       Im(Sigma_quad)"
      Do J = 1, Next
         Do Iq = 1, Nmo
            Do Ip = 1, Nmo
               Ix = I_sx + 2*((J-1)*Nmo*Nmo + (Iq-1)*Nmo + Ip-1)
               Iy = I_sq + 2*((J-1)*Nmo*Nmo + (Iq-1)*Nmo + Ip-1)
               Write(88,"(3I4,4(1x,1PD20.12))") J,Ip,Iq,
     +            Work(Ix),Work(Ix+1),Work(Iy),Work(Iy+1)
            Enddo
         Enddo
      Enddo
      Close(Unit=88)
C
      Write(6,"(a)") " @-Scrnc_lehmann: wrote wcreduc_lehmann.dat,"//
     +               " sigma_lehmann_imag.dat"
C
      Return
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_lehmann_w(Nmo_pairs,Nroot,Vr,Vl,Om,Wp,Vs,
     +                           Wre,Wim,Wimmax)
C
C W_c(pq,rt;iw') = Sum_s [Vr(pq)Vl(rt)/(iw'-Om) - Vl(pq)Vr(rt)/(iw'+Om)]
C   Re = -Sum_s [Vr Vl^T + Vl Vr^T] Om/(w'^2+Om^2)
C   Im =  Sum_s [Vl Vr^T - Vr Vl^T] w'/(w'^2+Om^2)
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Vr(Nmo_pairs,Nroot),Vl(Nmo_pairs,Nroot),Om(Nroot)
      Dimension Vs(Nmo_pairs,Nroot)
      Dimension Wre(Nmo_pairs,Nmo_pairs),Wim(Nmo_pairs,Nmo_pairs)
C
      Do Is = 1, Nroot
         Fac = -Om(Is)/(Wp*Wp+Om(Is)*Om(Is))
         Do I = 1, Nmo_pairs
            Vs(I,Is) = Fac*Vr(I,Is)
         Enddo
      Enddo
      Call Dgemm("N","T",Nmo_pairs,Nmo_pairs,Nroot,1.0D0,Vs,
     +           Nmo_pairs,Vl,Nmo_pairs,0.0D0,Wre,Nmo_pairs)
      Call Dgemm("N","T",Nmo_pairs,Nmo_pairs,Nroot,1.0D0,Vl,
     +           Nmo_pairs,Vs,Nmo_pairs,1.0D0,Wre,Nmo_pairs)
C
      Do Is = 1, Nroot
         Fac = Wp/(Wp*Wp+Om(Is)*Om(Is))
         Do I = 1, Nmo_pairs
            Vs(I,Is) = Fac*Vr(I,Is)
         Enddo
      Enddo
      Call Dgemm("N","T",Nmo_pairs,Nmo_pairs,Nroot,1.0D0,Vl,
     +           Nmo_pairs,Vs,Nmo_pairs,0.0D0,Wim,Nmo_pairs)
      Call Dgemm("N","T",Nmo_pairs,Nmo_pairs,Nroot,-1.0D0,Vs,
     +           Nmo_pairs,Vl,Nmo_pairs,1.0D0,Wim,Nmo_pairs)
C
      Wimmax = 0.0D0
      Do J = 1, Nmo_pairs
         Do I = 1, Nmo_pairs
            Wimmax = Max(Wimmax,Abs(Wim(I,J)))
         Enddo
      Enddo
C
      Return
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_lehmann_sigma(Nmo,Nocc,Nmo_pairs,Nroot,Next,
     +                               Nquad,Vr,Vl,Om,Eps,Wext,Xnode,
     +                               Wnode,Wq,Sx,Sq)
C
C Sx: exact analytic Sigma_c(p,q;iw_j) (see Scrnc_lehmann header).
C Sq: the pipeline's quadrature formula (scrnc_sigma_driver.F, the
C     two-denominator fold) fed with Re W_c (Lehmann) at its nodes --
C     reproduces sigma_c_output.dat iff the two W_c agree.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Vr(Nmo_pairs,Nroot),Vl(Nmo_pairs,Nroot),Om(Nroot)
      Dimension Eps(Nmo),Wext(Next),Xnode(Nquad),Wnode(Nquad)
      Dimension Wq(Nmo,Nmo,Nmo,Nquad)
      Double Complex Sx(Nmo,Nmo,Next),Sq(Nmo,Nmo,Next),Z
C
      Pi = Acos(-1.0D0)
C
      Do J = 1, Next
         Z = Dcmplx(0.0D0,Wext(J))
         Do Iq = 1, Nmo
            Do Ip = 1, Nmo
               Sx(Ip,Iq,J) = Dcmplx(0.0D0,0.0D0)
               Sq(Ip,Iq,J) = Dcmplx(0.0D0,0.0D0)
            Enddo
         Enddo
C
         Do Is = 1, Nroot
            Do Im = 1, Nmo
               If (Im .LE. Nocc) Then
                  Z = Dcmplx(-Eps(Im)+Om(Is),Wext(J))
               Else
                  Z = Dcmplx(-Eps(Im)-Om(Is),Wext(J))
               Endif
               Z = 1.0D0/Z
               Do Iq = 1, Nmo
                  Iqm = (Iq-1)*Nmo + Im
                  Do Ip = 1, Nmo
                     Ipm = (Ip-1)*Nmo + Im
                     If (Im .LE. Nocc) Then
                        Sx(Ip,Iq,J) = Sx(Ip,Iq,J) +
     +                                Vr(Ipm,Is)*Vl(Iqm,Is)*Z
                     Else
                        Sx(Ip,Iq,J) = Sx(Ip,Iq,J) +
     +                                Vl(Ipm,Is)*Vr(Iqm,Is)*Z
                     Endif
                  Enddo
               Enddo
            Enddo
         Enddo
C
         Do K = 1, Nquad
            Factor = -Wnode(K)/(2.0D0*Pi)
            Do Im = 1, Nmo
               Z = 1.0D0/Dcmplx(-Eps(Im),Wext(J)+Xnode(K)) +
     +             1.0D0/Dcmplx(-Eps(Im),Wext(J)-Xnode(K))
               Do Iq = 1, Nmo
                  Do Ip = 1, Nmo
                     Sq(Ip,Iq,J) = Sq(Ip,Iq,J) +
     +                             Factor*Z*Wq(Ip,Iq,Im,K)
                  Enddo
               Enddo
            Enddo
         Enddo
      Enddo
C
      Return
      End
