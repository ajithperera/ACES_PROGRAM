










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
C   Zprint                   (optional, default 1.0D-3) smallest pole
C                            strength listed in qp_lehmann.dat
C   Isrc                     (optional, default 0) spectrum source:
C                            0 = eom_tdm.dat (xvee FULL_DIAG)
C                            1 = closed-shell direct RPA solved here
C                                (Scrnc_lehmann_drpa), the dRPA/DRCCD
C                                screening; eom_tdm.dat is not read
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Work(Maxcor)


c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end



c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end
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
      Zprint = 1.0D-3
      Isrc   = 0
      Read(84,*,End=5) Zprint
      Read(84,*,End=5) Isrc
    5 Continue
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
C dRPA: one singlet root per (i,a) pair; triplets carry no spin-summed
C density and drop out of W_c exactly as in the EOM route.
C
      If (Isrc .EQ. 1) Then
         Nroot = Nocc*(Nmo-Nocc)
         Write(6,"(a,I6)") " @-Scrnc_lehmann: direct RPA spectrum,"//
     +                     " singlet roots ",Nroot
         Go To 20
      Endif
C
C Bare (orbital-difference) spin-orbital states, one per (spin,i,a):
C the vertex variants must then reproduce GF2 exactly.
C
      If (Isrc .EQ. 2) Then
         Nroot  = 2*Nocc*(Nmo-Nocc)
         Iuhf_t = 1
         Write(6,"(a,I6)") " @-Scrnc_lehmann: bare spin-orbital"//
     +                     " states (GF2 test), roots ",Nroot
         Go To 20
      Endif
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
   20 Continue
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
C Exchange vertices X, X' (right/left), each from the alpha density
C (slot 1) and from the triplet part t = (alpha-beta)/2 (slot 2), and
C the alpha, beta and t densities of one state.
      Nx     = 2*Nmo_pairs*Nroot
      I_xr   = Iend
      I_xl   = I_xr   + Nx
      I_xpr  = I_xl   + Nx
      I_xpl  = I_xpr  + Nx
      I_gra  = I_xpl  + Nx
      I_gla  = I_gra  + Nmo_pairs
      I_grb  = I_gla  + Nmo_pairs
      I_glb  = I_grb  + Nmo_pairs
      I_gtr  = I_glb  + Nmo_pairs
      I_gtl  = I_gtr  + Nmo_pairs
      I_vd   = I_gtl  + Nmo_pairs
C Complex-pair flag and 2x2 block per root (eom_evali.dat).
      I_mb   = I_vd   + Nmo_pairs
      Iend   = I_mb   + 5*Nroot
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
C dRPA source: V^R_s = V^L_s from the direct RPA X+Y, in the I_scr
C scratch (free again after the AO->MO transform).
C
      If (Isrc .EQ. 1) Then
         Nov    = Nroot
         Lwork  = 64*Nov
         I_zm   = I_scr
         I_gm   = I_zm + Nov*Nov
         I_dq   = I_gm + Nmo_pairs
         I_wk   = I_dq + Nov
         If (I_wk+Lwork .GT. Maxcor) Call Insmem("Scrnc_lehmann_drpa",
     +                                           I_wk+Lwork,Maxcor)
         Call Scrnc_lehmann_drpa(Nmo,Nocc,Nmo_pairs,Work(I_vmo),
     +                           Work(I_eps),Nov,Work(I_om),
     +                           Work(I_vr),Work(I_vl),Work(I_zm),
     +                           Work(I_gm),Work(I_dq),Work(I_wk),
     +                           Lwork)
         Go To 30
      Endif
C
      If (Isrc .EQ. 2) Then
         Nvrt = Nmo - Nocc
         Iroot = 0
         Do Isp = 1, 2
         Do Ii = 1, Nocc
         Do Ia = Nocc+1, Nmo
            Iroot = Iroot + 1
            Work(I_om-1+Iroot) = Work(I_eps-1+Ia) - Work(I_eps-1+Ii)
            Call Zero(Work(I_gr),Nmo_pairs)
            Call Zero(Work(I_gl),Nmo_pairs)
            Work(I_gr-1+(Ii-1)*Nmo+Ia) = 1.0D0
            Work(I_gl-1+(Ia-1)*Nmo+Ii) = 1.0D0
            Call Zero(Work(I_gra),Nmo_pairs)
            Call Zero(Work(I_gla),Nmo_pairs)
            Call Zero(Work(I_grb),Nmo_pairs)
            Call Zero(Work(I_glb),Nmo_pairs)
            If (Isp .EQ. 1) Then
               Call Dcopy(Nmo_pairs,Work(I_gr),1,Work(I_gra),1)
               Call Dcopy(Nmo_pairs,Work(I_gl),1,Work(I_gla),1)
            Else
               Call Dcopy(Nmo_pairs,Work(I_gr),1,Work(I_grb),1)
               Call Dcopy(Nmo_pairs,Work(I_gl),1,Work(I_glb),1)
            Endif
            Ioff = (Iroot-1)*Nmo_pairs
            Call Scrnc_lehmann_vx(Nmo,Nmo_pairs,Work(I_vmo),
     +           Work(I_gr),Work(I_gl),Work(I_gra),Work(I_gla),
     +           Work(I_vr+Ioff),Work(I_vl+Ioff),Work(I_xr+Ioff),
     +           Work(I_xl+Ioff),Work(I_xpr+Ioff),Work(I_xpl+Ioff))
            Call Scrnc_lehmann_vt(Nmo,Nmo_pairs,Nroot,Ioff,Work(I_vmo),
     +           Work(I_gra),Work(I_gla),Work(I_grb),Work(I_glb),
     +           Work(I_gtr),Work(I_gtl),Work(I_vd),Work(I_xr),
     +           Work(I_xl),Work(I_xpr),Work(I_xpl))
         Enddo
         Enddo
         Enddo
         Go To 30
      Endif
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
               If (Ispin .EQ. 1) Call Dcopy(Nmo_pairs,Work(I_gt),1,
     +                                      Work(I_gra),1)
               If (Ispin .EQ. 2) Call Dcopy(Nmo_pairs,Work(I_gt),1,
     +                                      Work(I_grb),1)
            Else
               Call Daxpy(Nmo_pairs,1.0D0,Work(I_gt),1,Work(I_gl),1)
               If (Ispin .EQ. 1) Call Dcopy(Nmo_pairs,Work(I_gt),1,
     +                                      Work(I_gla),1)
               If (Ispin .EQ. 2) Call Dcopy(Nmo_pairs,Work(I_gt),1,
     +                                      Work(I_glb),1)
            Endif
         Enddo
C        Closed-shell RHF-based EOM stores one spin case only (singlets;
C        the triplet part of the vertex variants is then missing).
         If (Iuhf_t .EQ. 0) Then
            Call Dscal(Nmo_pairs,2.0D0,Work(I_gr),1)
            Call Dscal(Nmo_pairs,2.0D0,Work(I_gl),1)
            Call Dcopy(Nmo_pairs,Work(I_gra),1,Work(I_grb),1)
            Call Dcopy(Nmo_pairs,Work(I_gla),1,Work(I_glb),1)
         Endif
         Work(I_om-1+Iroot) = Root
         Ioff = (Iroot-1)*Nmo_pairs
         Call Scrnc_lehmann_vx(Nmo,Nmo_pairs,Work(I_vmo),
     +        Work(I_gr),Work(I_gl),Work(I_gra),Work(I_gla),
     +        Work(I_vr+Ioff),Work(I_vl+Ioff),Work(I_xr+Ioff),
     +        Work(I_xl+Ioff),Work(I_xpr+Ioff),Work(I_xpl+Ioff))
         Call Scrnc_lehmann_vt(Nmo,Nmo_pairs,Nroot,Ioff,Work(I_vmo),
     +        Work(I_gra),Work(I_gla),Work(I_grb),Work(I_glb),
     +        Work(I_gtr),Work(I_gtl),Work(I_vd),Work(I_xr),Work(I_xl),
     +        Work(I_xpr),Work(I_xpl))
      Enddo
      Close(Unit=83)
   30 Continue
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
C
C Real-axis diagonal Dyson solve, all roots (QP + satellites), reusing
C the (now free) W_c/Sigma scratch from I_wre on.
C
      Npole = Nroot*Nmo
      I_pp  = I_wre
      I_pc  = I_pp + Npole
      I_ix  = I_pc + Npole
      Iend  = I_ix + Npole
      If (Iend .GT. Maxcor) Call Insmem("Scrnc_lehmann_qp",Iend,
     +                                  Maxcor)
      Call Scrnc_lehmann_qp(Nmo,Nocc,Nmo_pairs,Nroot,Npole,Zprint,
     +                      Work(I_vr),Work(I_vl),Work(I_om),
     +                      Work(I_eps),Work(I_pp),Work(I_pc),
     +                      Work(I_ix))
C
C
C Self-energy vertex variants (direct x exchange, antisymmetrized);
C dRPA supplies no transition densities, so they are skipped there.
C
      If (Isrc .NE. 1) Call Scrnc_lehmann_vtx(Nmo,Nocc,Nmo_pairs,
     +   Nroot,Isrc,Work(I_vr),Work(I_vl),Work(I_xr),Work(I_xl),
     +   Work(I_xpr),Work(I_xpl),Work(I_om),Work(I_eps),Work(I_vmo))
C
C Complex EOM pairs: same three variants with each pair entering
C through its exact 2x2 block, to qp_vertex_cpx.dat.
C
      If (Isrc .EQ. 0) Then
         Call Scrnc_lehmann_rdcpx(Nroot,Work(I_om),Work(I_mb),Npair)
         If (Npair .GT. 0) Call Scrnc_lehmann_vtxc(Nmo,Nocc,
     +      Nmo_pairs,Nroot,Work(I_vr),Work(I_vl),Work(I_xr),
     +      Work(I_xl),Work(I_xpr),Work(I_xpl),Work(I_om),Work(I_eps),
     +      Work(I_mb))
      Endif
C
      Write(6,"(a)") " @-Scrnc_lehmann: wrote wcreduc_lehmann.dat,"//
     +               " sigma_lehmann_imag.dat, qp_lehmann.dat"
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
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_lehmann_qp(Nmo,Nocc,Nmo_pairs,Nroot,Npole,
     +                            Zprint,Vr,Vl,Om,Eps,Pole,Res,Ix)
C
C Diagonal Dyson equation on the real axis with the exact Lehmann
C self-energy, Sigma_pp(w) = Sum_k c_k/(w-P_k):
C   P = eps_m - Omega_s, c = V^R_s(pm)V^L_s(pm)   (m occupied)
C   P = eps_m + Omega_s, c = V^L_s(pm)V^R_s(pm)   (m virtual)
C F(w) = w - eps_p - Sigma_pp(w) is rational, so every root is found
C by bracketing between neighbouring poles (plus one below the lowest
C and one above the highest); each gets the pole strength
C Z = 1/(1 + Sum_k c_k/(w-P_k)^2). For c_k >= 0 there is exactly one
C root per interval and Sum Z = 1 exactly -- printed as a check. Roots
C with Z >= Zprint go to qp_lehmann.dat; the largest-Z root per MO is
C its quasiparticle, the rest are satellites. No Pade, no quadrature.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Vr(Nmo_pairs,Nroot),Vl(Nmo_pairs,Nroot),Om(Nroot)
      Dimension Eps(Nmo),Pole(Npole),Res(Npole)
      Integer Ix(Npole)
C
      Parameter (Tiny=1.0D-14, Wout=50.0D0, Tol=1.0D-11)
C
      Open(Unit=89,File="qp_lehmann.dat",Form="FORMATTED",
     +     Status="UNKNOWN")
      Write(89,"(a)") "# p   eps_HF        root (Ha)      root (eV)"//
     +   "     Z         [QP = largest Z per p]"
      Write(6,"(a)") " @-Scrnc_lehmann_qp: p  eps_HF   eps_QP(Ha)"//
     +   "   Z_QP    Sum Z   Nroots  min c"
C
      Do Ip = 1, Nmo
C
C        Collect the non-negligible poles, sorted ascending.
C
         Np = 0
         Cmin = 0.0D0
         Do Is = 1, Nroot
            Do Im = 1, Nmo
               Ipm = (Ip-1)*Nmo + Im
               C = Vr(Ipm,Is)*Vl(Ipm,Is)
               If (Abs(C) .GT. Tiny) Then
                  Np = Np + 1
                  Res(Np) = C
                  If (Im .LE. Nocc) Then
                     Pole(Np) = Eps(Im) - Om(Is)
                  Else
                     Pole(Np) = Eps(Im) + Om(Is)
                  Endif
                  Cmin = Min(Cmin,C)
               Endif
            Enddo
         Enddo
         Call Scrnc_lehmann_sort(Np,Pole,Res,Ix)
C
         Zsum  = 0.0D0
         Zbest = -1.0D0
         Ebest = Eps(Ip)
         Nfound = 0
         Do Iint = 0, Np
            If (Iint .EQ. 0) Then
               A = Pole(1) - Wout
            Else
               A = Pole(Iint)
            Endif
            If (Iint .EQ. Np) Then
               B = Pole(Np) + Wout
            Else
               B = Pole(Iint+1)
            Endif
            Gap = B - A
            If (Gap .LT. 1.0D-12) Goto 100
            Xa = A + 1.0D-9*Gap
            Xb = B - 1.0D-9*Gap
            If (Iint .EQ. 0)  Xa = A
            If (Iint .EQ. Np) Xb = B
            Fa = Scrnc_lehmann_f(Xa,Eps(Ip),Np,Pole,Res,Dum)
            Fb = Scrnc_lehmann_f(Xb,Eps(Ip),Np,Pole,Res,Dum)
            If ((Fa.LT.0.0D0) .EQV. (Fb.LT.0.0D0)) Goto 100
            Do Iter = 1, 200
               Xm = 0.5D0*(Xa+Xb)
               Fm = Scrnc_lehmann_f(Xm,Eps(Ip),Np,Pole,Res,Dum)
               If ((Fm.LT.0.0D0) .EQV. (Fa.LT.0.0D0)) Then
                  Xa = Xm
                  Fa = Fm
               Else
                  Xb = Xm
               Endif
               If (Xb-Xa .LT. Tol) Goto 50
            Enddo
   50       Continue
            Xr = 0.5D0*(Xa+Xb)
            Fm = Scrnc_lehmann_f(Xr,Eps(Ip),Np,Pole,Res,Dsig)
            Z  = 1.0D0/(1.0D0+Dsig)
            Zsum = Zsum + Z
            Nfound = Nfound + 1
            If (Z .GT. Zbest) Then
               Zbest = Z
               Ebest = Xr
            Endif
            If (Z .GE. Zprint) Write(89,"(I4,1x,F12.6,1x,F14.8,1x,"//
     +         "F12.5,1x,F10.6)") Ip,Eps(Ip),Xr,Xr*27.2113957D0,Z
  100       Continue
         Enddo
         Write(6,"(a,I3,2F10.5,2F9.5,I7,1PD10.2)")
     +      " @-Scrnc_lehmann_qp:",Ip,Eps(Ip),Ebest,Zbest,Zsum,
     +      Nfound,Cmin
      Enddo
      Close(Unit=89)
C
      Return
      End
C
C-----------------------------------------------------------------
C
      Double Precision Function Scrnc_lehmann_f(W,Ep,Np,Pole,Res,
     +                                          Dsig)
C
C F(w) = w - eps_p - Sum_k c_k/(w-P_k); Dsig = Sum_k c_k/(w-P_k)^2
C (= -dSigma/dw), so the pole strength is 1/(1+Dsig).
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Pole(Np),Res(Np)
      Sig  = 0.0D0
      Dsig = 0.0D0
      Do K = 1, Np
         D    = 1.0D0/(W-Pole(K))
         Sig  = Sig  + Res(K)*D
         Dsig = Dsig + Res(K)*D*D
      Enddo
      Scrnc_lehmann_f = W - Ep - Sig
      Return
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_lehmann_sort(N,Pole,Res,Ix)
C
C Ascending sort of Pole (Res carried along) via an index heapsort,
C then exact-coincidence merge is unnecessary: equal poles give a
C zero-width interval, skipped by the caller.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Pole(N),Res(N)
      Integer Ix(N)
      Do I = 1, N
         Ix(I) = I
      Enddo
      If (N .LT. 2) Return
      L  = N/2 + 1
      Ir = N
   10 Continue
      If (L .GT. 1) Then
         L = L - 1
         Indxt = Ix(L)
         Q = Pole(Indxt)
      Else
         Indxt = Ix(Ir)
         Q = Pole(Indxt)
         Ix(Ir) = Ix(1)
         Ir = Ir - 1
         If (Ir .EQ. 1) Then
            Ix(1) = Indxt
            Goto 30
         Endif
      Endif
      I = L
      J = L + L
   20 If (J .LE. Ir) Then
         If (J .LT. Ir) Then
            If (Pole(Ix(J)) .LT. Pole(Ix(J+1))) J = J + 1
         Endif
         If (Q .LT. Pole(Ix(J))) Then
            Ix(I) = Ix(J)
            I = J
            J = J + J
         Else
            J = Ir + 1
         Endif
         Goto 20
      Endif
      Ix(I) = Indxt
      Goto 10
   30 Continue
C     Apply the permutation (Res/Pole copies via the tail of Ix order).
      Call Scrnc_lehmann_perm(N,Pole,Ix)
      Call Scrnc_lehmann_perm(N,Res,Ix)
      Return
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_lehmann_perm(N,A,Ix)
C
C A <- A(Ix(:)), in place, by cycle-following with a visited mark
C (Ix restored afterwards).
C
      Implicit Double Precision (A-H,O-Z)
      Dimension A(N)
      Integer Ix(N)
      Do I = 1, N
         If (Ix(I) .GT. 0) Then
            Tmp = A(I)
            J = I
   10       Continue
            K = Ix(J)
            Ix(J) = -K
            If (K .EQ. I) Then
               A(J) = Tmp
            Else
               A(J) = A(K)
               J = K
               Goto 10
            Endif
         Endif
      Enddo
      Do I = 1, N
         Ix(I) = -Ix(I)
      Enddo
      Return
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_lehmann_drpa(Nmo,Nocc,Nmo_pairs,Vmo,Eps,Nov,
     +                              Om,Vr,Vl,Zm,Gam,Dsq,Work,Lwork)
C
C Closed-shell direct RPA (no exchange) in the spatial MO basis, i.e.
C the screening of the dRPA/DRCCD pipeline. With D(ia) = eps_a-eps_i
C and K(ia,jb) = (ia|jb), singlet A = D + 2K, B = 2K, so A-B = D is
C diagonal and the problem reduces to the symmetric
C   D^1/2 (D + 4K) D^1/2 Z = Omega^2 Z,   (X+Y) = D^1/2 Z Omega^-1/2
C (X^T X - Y^T Y = 1). The spin-summed transition density of root s is
C sqrt(2)(X+Y)(ia); it is split evenly over the (i,a) and (a,i) pairs
C ((pq|rs) is symmetric in r,s), and V^R_s = V^L_s = (pq|gamma_s).
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Vmo(Nmo_pairs,Nmo_pairs),Eps(Nmo),Om(Nov)
      Dimension Vr(Nmo_pairs,Nov),Vl(Nmo_pairs,Nov)
      Dimension Zm(Nov,Nov),Gam(Nmo_pairs),Dsq(Nov),Work(Lwork)
C
      Nvrt = Nmo - Nocc
      Do Ia = 1, Nvrt
         Do I = 1, Nocc
            Iov      = (Ia-1)*Nocc + I
            Dsq(Iov) = Sqrt(Eps(Nocc+Ia) - Eps(I))
         Enddo
      Enddo
C
      Do Jb = 1, Nvrt
         Do J = 1, Nocc
            Jov  = (Jb-1)*Nocc + J
            Jpr  = (Nocc+Jb-1)*Nmo + J
            Do Ia = 1, Nvrt
               Do I = 1, Nocc
                  Iov = (Ia-1)*Nocc + I
                  Ipr = (Nocc+Ia-1)*Nmo + I
                  Zm(Iov,Jov) = 4.0D0*Dsq(Iov)*Vmo(Ipr,Jpr)*Dsq(Jov)
               Enddo
            Enddo
            Zm(Jov,Jov) = Zm(Jov,Jov) + Dsq(Jov)**4
         Enddo
      Enddo
C
C
C Tr A = Sum_ia [D(ia) + 2K(ia,ia)], for the plasmon-formula check
C E_c(dRPA) = (Sum_s Omega_s - Tr A)/2 = E(drCCD) below.
C
      Tra = 0.0D0
      Do Iov = 1, Nov
         D   = Dsq(Iov)**2
         Tra = Tra + D + 2.0D0*(Zm(Iov,Iov)-D*D)/(4.0D0*D)
      Enddo
C
      Call Dsyev("V","U",Nov,Zm,Nov,Om,Work,Lwork,Info)
      If (Info .NE. 0) Then
         Write(6,"(a,I6)") " @-Scrnc_lehmann_drpa: DSYEV info ",Info
         Call Errex
      Endif
      If (Om(1) .LE. 0.0D0) Then
         Write(6,"(a,1PD12.4)") " @-Scrnc_lehmann_drpa: unstable,"//
     +                          " lowest Omega^2 = ",Om(1)
         Call Errex
      Endif
C
      Write(6,"(a)") " @-Scrnc_lehmann_drpa: singlet dRPA excitation"//
     +               " energies (eV)"
      Do Is = 1, Nov
         Om(Is) = Sqrt(Om(Is))
      Enddo
      Write(6,"(5(1x,F14.8))") (Om(Is)*27.2113957D0,Is=1,Nov)
      Ec = -0.5D0*Tra
      Do Is = 1, Nov
         Ec = Ec + 0.5D0*Om(Is)
      Enddo
      Write(6,"(a,F20.12)") " @-Scrnc_lehmann_drpa: plasmon-formula"//
     +                      " E_c(dRPA) = ",Ec
C
      Do Is = 1, Nov
         Call Zero(Gam,Nmo_pairs)
         Fac = Sqrt(2.0D0/Om(Is))*0.5D0
         Do Ia = 1, Nvrt
            Do I = 1, Nocc
               Iov = (Ia-1)*Nocc + I
               X   = Fac*Dsq(Iov)*Zm(Iov,Is)
               Gam((Nocc+Ia-1)*Nmo + I) = X
               Gam((I-1)*Nmo + Nocc+Ia) = X
            Enddo
         Enddo
         Call Dgemv("N",Nmo_pairs,Nmo_pairs,1.0D0,Vmo,Nmo_pairs,
     +              Gam,1,0.0D0,Vr(1,Is),1)
         Call Dcopy(Nmo_pairs,Vr(1,Is),1,Vl(1,Is),1)
      Enddo
C
      Return
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_lehmann_vx(Nmo,Nmo_pairs,Vmo,Gr,Gl,Gra,Gla,
     +                            Vr,Vl,Xr,Xl,Xpr,Xpl)
C
C Direct and exchange vertices of one state, index (p-1)*Nmo+m:
C   V(pm)  = Sum_rs (pm|rs) g(rs)        (spin-summed g)
C   X(pm)  = Sum_rs (ps|rm) ga(rs)       (alpha g only)
C   X'(pm) = Sum_rs (pr|sm) ga(rs)
C with g(r,s) stored at (s-1)*Nmo+r, as in eom_tdm.dat.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Vmo(Nmo_pairs,Nmo_pairs),Gr(Nmo_pairs),Gl(Nmo_pairs)
      Dimension Gra(Nmo,Nmo),Gla(Nmo,Nmo),Vr(Nmo_pairs),Vl(Nmo_pairs)
      Dimension Xr(Nmo,Nmo),Xl(Nmo,Nmo),Xpr(Nmo,Nmo),Xpl(Nmo,Nmo)
C
      Call Dgemv("N",Nmo_pairs,Nmo_pairs,1.0D0,Vmo,Nmo_pairs,Gr,1,
     +           0.0D0,Vr,1)
      Call Dgemv("N",Nmo_pairs,Nmo_pairs,1.0D0,Vmo,Nmo_pairs,Gl,1,
     +           0.0D0,Vl,1)
      Do Ip = 1, Nmo
         Do Im = 1, Nmo
            Sr  = 0.0D0
            Sl  = 0.0D0
            Spr = 0.0D0
            Spl = 0.0D0
            Do Is = 1, Nmo
               Do Ir = 1, Nmo
                  Ex  = Vmo((Is-1)*Nmo+Ip,(Im-1)*Nmo+Ir)
                  Exq = Vmo((Ir-1)*Nmo+Ip,(Im-1)*Nmo+Is)
                  Sr  = Sr  + Ex *Gra(Ir,Is)
                  Sl  = Sl  + Ex *Gla(Ir,Is)
                  Spr = Spr + Exq*Gra(Ir,Is)
                  Spl = Spl + Exq*Gla(Ir,Is)
               Enddo
            Enddo
            Xr(Im,Ip)  = Sr
            Xl(Im,Ip)  = Sl
            Xpr(Im,Ip) = Spr
            Xpl(Im,Ip) = Spl
         Enddo
      Enddo
C
      Return
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_lehmann_res(Iv,Im,Nocc,Ipm,Is,Nmo_pairs,Nroot,
     +                             Vr,Vl,Xr,Xl,Xpr,Xpl,Om,Eps,C,P)
C
C Pole P and residue C of Sigma_pp for state Is and orbital Im:
C   Iv=1 (a) direct x direct          (G0W)
C   Iv=2 (b) direct x (direct - exchange)
C   Iv=3 (c) 1/2 (direct - exchange) x (direct - exchange)
C (b) and (c) are both exact GF2 for bare states; (c) also carries
C the triplet channel: X x X' from the Ms=0 states, plus the spin-flip
C (Ms=+-1) states, which for a closed-shell singlet reference give
C twice the triplet-part product At*Bt (spin isotropy).
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Vr(Nmo_pairs,Nroot),Vl(Nmo_pairs,Nroot)
      Dimension Xr(Nmo_pairs,Nroot,2),Xl(Nmo_pairs,Nroot,2)
      Dimension Xpr(Nmo_pairs,Nroot,2),Xpl(Nmo_pairs,Nroot,2)
      Dimension Om(Nroot),Eps(*)
C
      If (Im .LE. Nocc) Then
         P  = Eps(Im) - Om(Is)
         A  = Vr(Ipm,Is)
         B  = Vl(Ipm,Is)
         Ax = Xr(Ipm,Is,1)
         Bx = Xpl(Ipm,Is,1)
         At = Xr(Ipm,Is,2)
         Bt = Xpl(Ipm,Is,2)
      Else
         P  = Eps(Im) + Om(Is)
         A  = Vl(Ipm,Is)
         B  = Vr(Ipm,Is)
         Ax = Xl(Ipm,Is,1)
         Bx = Xpr(Ipm,Is,1)
         At = Xl(Ipm,Is,2)
         Bt = Xpr(Ipm,Is,2)
      Endif
      If (Iv .EQ. 1) Then
         C = A*B
      Else If (Iv .EQ. 2) Then
         C = A*(B-Bx)
      Else
         C = 0.5D0*((A-Ax)*(B-Bx) + 2.0D0*At*Bt)
      Endif
C
      Return
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_lehmann_vtx(Nmo,Nocc,Nmo_pairs,Nroot,Isrc,
     +                             Vr,Vl,Xr,Xl,Xpr,Xpl,Om,Eps,Vmo)
C
C Quasiparticle energies with the three vertex variants (diagonal
C Newton solve from eps_p, Z = 1/(1 - dSigma/dw)), to qp_vertex.dat.
C Isrc=2 (bare states): also checks (a) against the direct part of
C GF2 and (b), (c) against full GF2 at z = eps_p + 0.1i.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Vr(Nmo_pairs,Nroot),Vl(Nmo_pairs,Nroot)
      Dimension Xr(Nmo_pairs,Nroot,2),Xl(Nmo_pairs,Nroot,2)
      Dimension Xpr(Nmo_pairs,Nroot,2),Xpl(Nmo_pairs,Nroot,2)
      Dimension Om(Nroot),Eps(Nmo),Vmo(Nmo_pairs,Nmo_pairs)
      Dimension Qp(3),Zq(3)
      Double Complex Z,S(3),Gd,Gx
C
      Hev = 27.211386245988D0
      Open(Unit=89,File="qp_vertex.dat",Form="FORMATTED",
     +     Status="UNKNOWN")
      Write(89,"(a)") "# p   eps(eV)   QP_a(eV)   Z_a   QP_b(eV)"//
     +  "   Z_b   QP_c(eV)   Z_c   (a GW, b +exchange vertex,"//
     +  " c antisymmetrized)"
C
      Dmax = 0.0D0
      Do Ip = 1, Nmo
         Do Iv = 1, 3
            W = Eps(Ip)
            Do It = 1, 200
               Sg = 0.0D0
               Dg = 0.0D0
               Do Is = 1, Nroot
                  Do Im = 1, Nmo
                     Ipm = (Ip-1)*Nmo + Im
                     Call Scrnc_lehmann_res(Iv,Im,Nocc,Ipm,Is,
     +                    Nmo_pairs,Nroot,Vr,Vl,Xr,Xl,Xpr,Xpl,Om,
     +                    Eps,C,P)
                     Sg = Sg + C/(W-P)
                     Dg = Dg + C/((W-P)*(W-P))
                  Enddo
               Enddo
               F  = W - Eps(Ip) - Sg
               Fp = 1.0D0 + Dg
               Dw = -F/Fp
               If (Abs(Dw) .GT. 0.05D0) Dw = Sign(0.05D0,Dw)
               W = W + Dw
               If (Abs(Dw) .LT. 1.0D-10) Go To 40
            Enddo
            Write(6,"(a,2I4)") " @-Scrnc_lehmann_vtx: no convergence,"//
     +                         " p, variant ",Ip,Iv
   40       Continue
            Qp(Iv) = W
            Zq(Iv) = 1.0D0/Fp
         Enddo
         Write(89,"(I4,7F11.5)") Ip,Eps(Ip)*Hev,(Qp(K)*Hev,Zq(K),K=1,3)
C
         If (Isrc .EQ. 2) Then
            Z = Dcmplx(Eps(Ip),0.1D0)
            Do Iv = 1, 3
               S(Iv) = Dcmplx(0.0D0,0.0D0)
               Do Is = 1, Nroot
                  Do Im = 1, Nmo
                     Ipm = (Ip-1)*Nmo + Im
                     Call Scrnc_lehmann_res(Iv,Im,Nocc,Ipm,Is,
     +                    Nmo_pairs,Nroot,Vr,Vl,Xr,Xl,Xpr,Xpl,Om,
     +                    Eps,C,P)
                     S(Iv) = S(Iv) + C/(Z-P)
                  Enddo
               Enddo
            Enddo
            Call Scrnc_lehmann_gf2(Nmo,Nocc,Nmo_pairs,Ip,Z,Eps,Vmo,
     +                             Gd,Gx)
            D1 = Abs(S(1)-Gd)
            D2 = Abs(S(2)-Gd-Gx)
            D3 = Abs(S(3)-Gd-Gx)
            Write(6,"(a,I4,3(1PD12.3))") " @-Scrnc_lehmann_vtx: GF2"//
     +         " test p, |a-GF2dir| |b-GF2| |c-GF2| ",Ip,D1,D2,D3
            Dmax = Max(Dmax,D1,D2,D3)
         Endif
      Enddo
      Close(Unit=89)
      If (Isrc .EQ. 2) Write(6,"(a,1PD12.3)")
     +   " @-Scrnc_lehmann_vtx: GF2 test, max deviation ",Dmax
C
      Return
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_lehmann_gf2(Nmo,Nocc,Nmo_pairs,Ip,Z,Eps,Vmo,
     +                             Gd,Gx)
C
C Closed-shell second-order self-energy Sigma_pp(z), split into its
C direct (Gd) and exchange (Gx) parts:
C   Sum_ijb (pi|bj)[2(pi|bj)-(pj|bi)]/(z+e_b-e_i-e_j)
C + Sum_iab (pa|ib)[2(pa|ib)-(pb|ia)]/(z+e_i-e_a-e_b)
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Eps(Nmo),Vmo(Nmo_pairs,Nmo_pairs)
      Double Complex Z,Gd,Gx,Rd
C
      Gd = Dcmplx(0.0D0,0.0D0)
      Gx = Dcmplx(0.0D0,0.0D0)
      Do Ii = 1, Nocc
         Do Ij = 1, Nocc
            Do Ib = Nocc+1, Nmo
               V1 = Vmo((Ii-1)*Nmo+Ip,(Ij-1)*Nmo+Ib)
               V2 = Vmo((Ij-1)*Nmo+Ip,(Ii-1)*Nmo+Ib)
               Rd = 1.0D0/(Z+Eps(Ib)-Eps(Ii)-Eps(Ij))
               Gd = Gd + 2.0D0*V1*V1*Rd
               Gx = Gx - V1*V2*Rd
            Enddo
         Enddo
      Enddo
      Do Ii = 1, Nocc
         Do Ia = Nocc+1, Nmo
            Do Ib = Nocc+1, Nmo
               V1 = Vmo((Ia-1)*Nmo+Ip,(Ib-1)*Nmo+Ii)
               V2 = Vmo((Ib-1)*Nmo+Ip,(Ia-1)*Nmo+Ii)
               Rd = 1.0D0/(Z+Eps(Ii)-Eps(Ia)-Eps(Ib))
               Gd = Gd + 2.0D0*V1*V1*Rd
               Gx = Gx - V1*V2*Rd
            Enddo
         Enddo
      Enddo
C
      Return
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_lehmann_vt(Nmo,Nmo_pairs,Nroot,Ioff,Vmo,Gra,
     +                            Gla,Grb,Glb,Gtr,Gtl,Vd,Xr,Xl,Xpr,Xpl)
C
C Exchange vertices of the triplet part t = (alpha-beta)/2 of one
C state's densities, into slot 2 of Xr, Xl, Xpr, Xpl (column Ioff).
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Vmo(Nmo_pairs,Nmo_pairs),Gra(Nmo_pairs),Gla(Nmo_pairs)
      Dimension Grb(Nmo_pairs),Glb(Nmo_pairs),Gtr(Nmo_pairs)
      Dimension Gtl(Nmo_pairs),Vd(Nmo_pairs)
      Dimension Xr(*),Xl(*),Xpr(*),Xpl(*)
C
      Do I = 1, Nmo_pairs
         Gtr(I) = 0.5D0*(Gra(I)-Grb(I))
         Gtl(I) = 0.5D0*(Gla(I)-Glb(I))
      Enddo
      I2 = Nmo_pairs*Nroot + Ioff + 1
      Call Scrnc_lehmann_vx(Nmo,Nmo_pairs,Vmo,Gtr,Gtl,Gtr,Gtl,Vd,Vd,
     +                      Xr(I2),Xl(I2),Xpr(I2),Xpl(I2))
C
      Return
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_lehmann_rdcpx(Nroot,Om,Fm,Npair)
C
C Reads eom_evali.dat (xvee DIAGALL, one line per root, in eom_tdm.dat
C order). Fm(1,Is) = 1 marks Is as the 2nd member w of a complex pair
C (u = Is-1, w = Is); Fm(2:5,Is) = Muu, Muw, Mwu, Mww, the block
C <L_k|Hbar|R_l> of the stored vectors. No file: no pairs.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Om(Nroot),Fm(5,Nroot)
      Logical Yes
C
      Npair = 0
      Call Zero(Fm,5*Nroot)
      Inquire(File="eom_evali.dat",Exist=Yes)
      If (.Not. Yes) Return
      Open(Unit=86,File="eom_evali.dat",Form="FORMATTED",
     +     Status="OLD")
      Dmax = 0.0D0
      Do Is = 1, Nroot
         Read(86,*) Ir,I,Ev,Evi,(Fm(K,Is),K=2,5)
         If (Abs(Ev-Om(Is)) .GT. 1.0D-8) Then
            Write(6,"(a,I7,2F16.10)") " @-Scrnc_lehmann_rdcpx: root"//
     +         " order mismatch ",Is,Ev,Om(Is)
            Call Errex
         Endif
         If (Evi .LT. 0.0D0 .AND. Is .GT. 1) Then
            Npair = Npair + 1
            Fm(1,Is) = 1.0D0
            Write(6,"(a,I7,a,F12.6,a,1PE10.3,a,2(1PE11.3),a,"//
     +         "1PE9.2)") " @-Scrnc_lehmann_rdcpx: pair at root ",Is,
     +         " Re w (eV) ",Ev*27.211386245988D0," Im w ",-Evi,
     +         " Muw Mwu ",Fm(3,Is),Fm(4,Is)," |Muu-Re w| ",
     +         Abs(Fm(2,Is)-Ev)
            Dmax = Max(Dmax,Abs(Fm(2,Is)-Ev),Abs(Fm(5,Is)-Ev))
         Endif
      Enddo
      Close(Unit=86)
      Write(6,"(a,I6,a,1PE10.3)") " @-Scrnc_lehmann_rdcpx: complex"//
     +   " pairs ",Npair,"  max |Mkk - Re w| ",Dmax
C
      Return
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_lehmann_fac(Iv,Im,Nocc,Ipm,Is,Nmo_pairs,Nroot,
     +                             Vr,Vl,Xr,Xl,Xpr,Xpl,Nt,Cf,Fr,Fl)
C
C The residue of Scrnc_lehmann_res as Sum_t Cf(t) Fr(t) Fl(t), with Fr
C built from the right and Fl from the left transition density of
C state Is, so that a 2x2 block can couple Fr of one pair member with
C Fl of the other.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Vr(Nmo_pairs,Nroot),Vl(Nmo_pairs,Nroot)
      Dimension Xr(Nmo_pairs,Nroot,2),Xl(Nmo_pairs,Nroot,2)
      Dimension Xpr(Nmo_pairs,Nroot,2),Xpl(Nmo_pairs,Nroot,2)
      Dimension Cf(2),Fr(2),Fl(2)
C
      Nt    = 1
      Cf(1) = 1.0D0
      If (Im .LE. Nocc) Then
         Fr(1) = Vr(Ipm,Is)
         Fl(1) = Vl(Ipm,Is)
         If (Iv .GE. 2) Fl(1) = Vl(Ipm,Is) - Xpl(Ipm,Is,1)
         If (Iv .EQ. 3) Then
            Fr(1) = Vr(Ipm,Is) - Xr(Ipm,Is,1)
            Fr(2) = Xr(Ipm,Is,2)
            Fl(2) = Xpl(Ipm,Is,2)
         Endif
      Else
         Fl(1) = Vl(Ipm,Is)
         Fr(1) = Vr(Ipm,Is)
         If (Iv .GE. 2) Fr(1) = Vr(Ipm,Is) - Xpr(Ipm,Is,1)
         If (Iv .EQ. 3) Then
            Fl(1) = Vl(Ipm,Is) - Xl(Ipm,Is,1)
            Fl(2) = Xl(Ipm,Is,2)
            Fr(2) = Xpr(Ipm,Is,2)
         Endif
      Endif
      If (Iv .EQ. 3) Then
         Nt    = 2
         Cf(1) = 0.5D0
         Cf(2) = 1.0D0
      Endif
C
      Return
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_lehmann_vtxc(Nmo,Nocc,Nmo_pairs,Nroot,
     +                              Vr,Vl,Xr,Xl,Xpr,Xpl,Om,Eps,Fm)
C
C As Scrnc_lehmann_vtx, but each complex pair (u,w) enters exactly:
C Sum_kl Fr_k G_kl Fl_l with G = (x + M)**-1, x = w - eps_m (m occ),
C or G = (y - M)**-1, y = w - eps_m (m virtual), M the 2x2 block of
C Hbar on the pair. Real states are treated as before. The pair sum
C is real for real w (it is the sum of the two complex-conjugate
C poles). To qp_vertex_cpx.dat, same layout as qp_vertex.dat.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Vr(Nmo_pairs,Nroot),Vl(Nmo_pairs,Nroot)
      Dimension Xr(Nmo_pairs,Nroot,2),Xl(Nmo_pairs,Nroot,2)
      Dimension Xpr(Nmo_pairs,Nroot,2),Xpl(Nmo_pairs,Nroot,2)
      Dimension Om(Nroot),Eps(Nmo),Fm(5,Nroot)
      Dimension Qp(3),Zq(3),G(2,2),G2(2,2)
      Dimension Cfu(2),Fru(2),Flu(2),Cfw(2),Frw(2),Flw(2)
C
      Hev = 27.211386245988D0
      Open(Unit=89,File="qp_vertex_cpx.dat",Form="FORMATTED",
     +     Status="UNKNOWN")
      Write(89,"(a)") "# p   eps(eV)   QP_a(eV)   Z_a   QP_b(eV)"//
     +  "   Z_b   QP_c(eV)   Z_c   (complex EOM pairs as exact 2x2"//
     +  " blocks)"
C
      Do Ip = 1, Nmo
         Do Iv = 1, 3
            W = Eps(Ip)
            Do It = 1, 200
               Sg = 0.0D0
               Dg = 0.0D0
               Do Is = 1, Nroot
                  Ipair = 0
                  If (Fm(1,Is) .GT. 0.5D0) Ipair = 1
                  If (Is .LT. Nroot) Then
                     If (Fm(1,Is+1) .GT. 0.5D0) Ipair = -1
                  Endif
                  If (Ipair .EQ. -1) Go To 30
                  Do Im = 1, Nmo
                     Ipm = (Ip-1)*Nmo + Im
                     If (Ipair .EQ. 0) Then
                        Call Scrnc_lehmann_res(Iv,Im,Nocc,Ipm,Is,
     +                       Nmo_pairs,Nroot,Vr,Vl,Xr,Xl,Xpr,Xpl,Om,
     +                       Eps,C,P)
                        Sg = Sg + C/(W-P)
                        Dg = Dg + C/((W-P)*(W-P))
                     Else
                        X = W - Eps(Im)
                        If (Im .LE. Nocc) Then
                           A = X + Fm(2,Is)
                           B = Fm(3,Is)
                           Cc = Fm(4,Is)
                           D = X + Fm(5,Is)
                        Else
                           A = X - Fm(2,Is)
                           B = -Fm(3,Is)
                           Cc = -Fm(4,Is)
                           D = X - Fm(5,Is)
                        Endif
                        Det = A*D - B*Cc
                        G(1,1) =  D/Det
                        G(1,2) = -B/Det
                        G(2,1) = -Cc/Det
                        G(2,2) =  A/Det
                        Do K = 1, 2
                           Do L = 1, 2
                              G2(K,L) = G(K,1)*G(1,L) + G(K,2)*G(2,L)
                           Enddo
                        Enddo
                        Call Scrnc_lehmann_fac(Iv,Im,Nocc,Ipm,Is-1,
     +                       Nmo_pairs,Nroot,Vr,Vl,Xr,Xl,Xpr,Xpl,
     +                       Nt,Cfu,Fru,Flu)
                        Call Scrnc_lehmann_fac(Iv,Im,Nocc,Ipm,Is,
     +                       Nmo_pairs,Nroot,Vr,Vl,Xr,Xl,Xpr,Xpl,
     +                       Nt,Cfw,Frw,Flw)
                        Do It2 = 1, Nt
                           Sg = Sg + Cfu(It2)*(
     +                          Fru(It2)*G(1,1)*Flu(It2)
     +                        + Fru(It2)*G(1,2)*Flw(It2)
     +                        + Frw(It2)*G(2,1)*Flu(It2)
     +                        + Frw(It2)*G(2,2)*Flw(It2))
                           Dg = Dg + Cfu(It2)*(
     +                          Fru(It2)*G2(1,1)*Flu(It2)
     +                        + Fru(It2)*G2(1,2)*Flw(It2)
     +                        + Frw(It2)*G2(2,1)*Flu(It2)
     +                        + Frw(It2)*G2(2,2)*Flw(It2))
                        Enddo
                     Endif
                  Enddo
   30             Continue
               Enddo
               F  = W - Eps(Ip) - Sg
               Fp = 1.0D0 + Dg
               Dw = -F/Fp
               If (Abs(Dw) .GT. 0.05D0) Dw = Sign(0.05D0,Dw)
               W = W + Dw
               If (Abs(Dw) .LT. 1.0D-10) Go To 40
            Enddo
            Write(6,"(a,2I4)") " @-Scrnc_lehmann_vtxc: no"//
     +                         " convergence, p, variant ",Ip,Iv
   40       Continue
            Qp(Iv) = W
            Zq(Iv) = 1.0D0/Fp
         Enddo
         Write(89,"(I4,7F11.5)") Ip,Eps(Ip)*Hev,(Qp(K)*Hev,Zq(K),K=1,3)
      Enddo
      Close(Unit=89)
C
      Return
      End
