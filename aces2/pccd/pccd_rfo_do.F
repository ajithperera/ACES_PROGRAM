      Subroutine Pccd_rfo_do(G,C,Cprv,F,U,Scr,Gv,Gp,Sp,Step,H,Ah,Ev,
     +                       Wk,Lwork,Vb,Sb,Am,Wm,Cs,Gs,Res,Hf,Irp,
     +                       Ih0,Mxhv,Nbas,Nocc,Npar,Grdcnv,Tol,Etol,
     +                       Gmax,Grms,Ncycle,Wr,Lrsp,Lah,Lhes)
C Lah: length of Ah; Lhes = 0: H and Ah are not Npar**2 arrays (large
C Npar with the exact Hessian, see Pccd_rfo_step), Ah is only the
C Davidson workspace.
      Implicit Double Precision(A-H,O-Z)
      Logical Grdcnv,Reject,Escape,Jump,Accept,Hok
      Dimension Vb(Npar,Mxhv),Sb(Npar,Mxhv),Am(Mxhv,Mxhv),Wm(Mxhv)
      Dimension Cs(Nbas*Nbas),Gs(Npar),Res(Npar)
      Dimension G(Nbas,Nbas),C(Nbas,Nbas),Cprv(Nbas,Nbas)
      Dimension F(Nbas,Nbas),U(Nbas,Nbas),Scr(Nbas,Nbas)
      Dimension Gv(Npar),Gp(Npar),Sp(Npar),Step(Npar)
      Dimension H(Npar,Npar),Ah(Npar+1,Npar+1),Ev(Npar+1),Wk(Lwork)
      Dimension Irp(Nbas),Hf(Npar,Npar),Wr(*),Gsh(2),Stp(2)
#include "machsp.com"
      Parameter (Hoovv=0.05D0, Trust0=0.1D0, Trstmx=0.5D0)
      Parameter (Ngvec=3)
      Parameter (Mxstb=5)
C Branch jump: an energy drop far beyond the model prediction means the
C amplitudes left the solution branch being followed (seen for stretched
C N2: a branch with |T2| ~ 30 and 0.23 Eh lower energy).
      Parameter (Rhojmp=10.0D0, Ejmp=1.0D-4)
C No step can be accepted once the trust radius has collapsed (N2
C cc-pVDZ 3.0 A from RHF orbitals: 600 macro iterations at trust 0,
C every point rejected by Pccd_rfo_ampchk); stop instead. No run of
C the 10/2026 N2/H2O surveys that went below 1e-6 recovered (all were
C the unphysical OO-pCCD branch of stretched N2, 2.52-3.0 A).
      Parameter (Trmin=1.0D-6)
      Data Zero,Half,One,Two,Four /0.0D0,0.5D0,1.0D0,2.0D0,4.0D0/

      Nsq = Nbas*Nbas
      Np1 = Npar + 1
C Idav = 1: Davidson instead of the dense diagonalizations (large Npar).
      Call Pccd_rfo_eigmode(Npar,Idav)
      Call Getrec(20,"JOBARC","OBRTGRDA",Nsq*Iintfp,G)
      Call Getrec(20,"JOBARC","TOTENERG",Iintfp,E)
      Call Getrec(20,"JOBARC","SCFEVCA0",Nsq*Iintfp,C)

      K = 0
      Do Ip = 1, Nbas
         Do Iq = Ip+1, Nbas
            If (Irp(Ip) .Eq. Irp(Iq)) Then
               K = K + 1
               Gv(K) = -Two*(G(Ip,Iq) - G(Iq,Ip))
            Endif
         Enddo
      Enddo

      Call Getrec(0,"JOBARC","RFO_NITR",Length,Junk)
      Nit = 0
      If (Length .Gt. 0) Call Getrec(20,"JOBARC","RFO_NITR",1,Nit)
      Call Getrec(0,"JOBARC","RFO_NTOT",Length,Junk)
      Ntot = 0
      If (Length .Gt. 0) Call Getrec(20,"JOBARC","RFO_NTOT",1,Ntot)
      Ntot = Ntot + 1
      Call Putrec(20,"JOBARC","RFO_NTOT",1,Ntot)
      Ncycle = Ntot
      Call Getrec(0,"JOBARC","RFO_PHAS",Length,Junk)
      Iphase = 0
      If (Length .Gt. 0) Call Getrec(20,"JOBARC","RFO_PHAS",1,Iphase)
      Reject = .False.
      Escape = .False.
      Jump   = .False.
C Iexp = 1 (2): the previous step came from the exact (frozen) Hessian
C (Ih0 = 3; no model Hessian is kept); Iescp = 1: it was a step along a
C negative mode found by the stability check.
      Iexp  = 0
      Iescp = 0
      If (Nit .Gt. 0) Then
         Call Getrec(0,"JOBARC","RFO_HEXA",Length,Junk)
         If (Length .Gt. 0) Call Getrec(20,"JOBARC","RFO_HEXA",1,Iexp)
         Call Getrec(0,"JOBARC","RFO_ESCP",Length,Junk)
         If (Length .Gt. 0) Call Getrec(20,"JOBARC","RFO_ESCP",1,Iescp)
      Endif
C RFO_FAIL = 1: the pCCD amplitude or lambda equations did not converge
C at these orbitals or reached an unphysical solution (Pccd_rfo_ampchk);
C set in Pccd. E and G are not available.
      Call Getrec(0,"JOBARC","RFO_FAIL",Length,Junk)
      Ifail = 0
      If (Length .Gt. 0) Call Getrec(20,"JOBARC","RFO_FAIL",1,Ifail)
      If (Ifail .Eq. 1) Then
         Call Putrec(20,"JOBARC","RFO_FAIL",1,0)
         If ((Iphase .Eq. 1 .Or. Nit .Eq. 0) .And. Iphase .Ne. 2) Then
            Write(6,"(2a)") " @-Pccd_rfo_do: pCCD amplitudes not",
     +           " converged or unphysical at an accepted/stability"//
     +           " point; stopping."
            Call Errex
         Endif
      Endif
      Iacpt  = 0
      Call Putrec(20,"JOBARC","RFO_ACPT",1,Iacpt)
      Rho    = Zero
C Analytic frozen-amplitude Hessian of this point (Pccd_ohess), when
C the model uses it (Ih0 = 1: first point only; 2: every accepted
C point) and it was formed (not at failed or Davidson points).
      Ifrz = 0
      If (Ih0 .Ge. 1 .And. Ifail .Eq. 0 .And. Iphase .Ne. 1) Then
         Call Getrec(0,"JOBARC","OHESOK  ",Length,Junk)
         If (Length .Gt. 0) Call Getrec(20,"JOBARC","OHESOK  ",1,Ifrz)
         If (Ifrz .Eq. 1) Call Getrec(20,"JOBARC","OHESFRZ ",
     +                                Npar*Npar*Iintfp,Hf)
      Endif

      If (Iphase .Eq. 2) Then
C Soft-mode probe (Pccd_rfo_prb): this point is a probe.
         Call Pccd_rfo_prb("CONT ",Ipr,C,E,Gv,Sb,Wm,0,Step,F,U,Scr,G,
     +                     Irp,Nbas,Npar,Ifail)
         If (Ipr .Eq. 0) Return
         If (Ipr .Eq. 1) Then
C Back at the minimum (the last xpccd run was at a probe).
            Grdcnv = .True.
            Gmax = Zero
            Grms = Zero
            Do K = 1, Npar
               Gmax = Max(Gmax,Dabs(Gv(K)))
               Grms = Grms + Gv(K)*Gv(K)
            Enddo
            Grms = Dsqrt(Grms/Dble(Npar))
            Call Putrec(20,"JOBARC","TOTENERG",Iintfp,E)
            Call Putrec(20,"JOBARC","SCFEVCA0",Nsq*Iintfp,C)
            Call Putrec(20,"JOBARC","SCFEVCB0",Nsq*Iintfp,C)
            Call Putrec(20,"JOBARC","GRADTEST",1,1)
            Write(6,"(2a,i5,a)") " @-RFO orbital optimization",
     +           " converged to a minimum in",Ntot," gradients."
            Return
         Endif
C A lower probe: a new start from it.
         Nit    = 0
         Iphase = 0
      Endif

      If (Iphase .Eq. 1) Then
C Stability phase: this gradient belongs to a displaced point of the
C Davidson iterations.
         Call Getrec(20,"JOBARC","RFO_HESS",Npar*Npar*Iintfp,H)
         Call Getrec(20,"JOBARC","RFO_TRST",Iintfp,Trust)
         Call Pccd_rfo_stab("CONT ",Iact,Alam,C,Gv,E,H,Vb,Sb,Am,Wm,Cs,
     +                      Gs,Step,Res,F,U,Scr,G,Wk,Lwork,Irp,Nbas,
     +                      Npar,Mxhv)
         Grdcnv = .False.
         If (Iact .Eq. 0) Then
            Call Putrec(20,"JOBARC","GRADTEST",1,0)
            Return
         Elseif (Iact .Eq. 1 .Or. Iact .Eq. 3) Then
C The last xpccd run was at a displaced point: put the energy and the
C gradient statistics of the stationary point back for Pccd_finalize.
C Iact = 3: the Davidson budget ran out without a negative eigenvalue;
C the stationary point is not verified to be a minimum.
            Grdcnv = .True.
            Gmax = Zero
            Grms = Zero
            Do K = 1, Npar
               Gmax = Max(Gmax,Dabs(Gv(K)))
               Grms = Grms + Gv(K)*Gv(K)
            Enddo
            Grms = Dsqrt(Grms/Dble(Npar))
            Call Putrec(20,"JOBARC","TOTENERG",Iintfp,E)
            Call Putrec(20,"JOBARC","SCFEVCA0",Nsq*Iintfp,C)
            Call Putrec(20,"JOBARC","SCFEVCB0",Nsq*Iintfp,C)
            Call Putrec(20,"JOBARC","GRADTEST",1,1)
            If (Iact .Eq. 1) Then
            Write(6,"(2a,i5,a)") " @-RFO orbital optimization",
     +             " converged to a minimum in",Ntot," gradients."
            Else
            Write(6,"(2a,i5,a)") " @-RFO orbital optimization",
     +             " converged in",Ntot," gradients; stability check"//
     +             " INCONCLUSIVE (possibly a saddle point)."
            Endif
            Return
         Endif
C Negative curvature: put the measured eigenvalue into the model Hessian
C along the mode and take a trust-region step from the stationary point.
         Call Dgemv("N",Npar,Npar,One,H,Npar,Step,1,Zero,Gp,1)
         Vhv = Ddot(Npar,Step,1,Gp,1)
         Call Dger(Npar,Npar,Alam-Vhv,Step,1,Step,1,H,Npar)
         Trust  = Trust0
         Escape = .True.
         Write(6,"(a,1p,E12.4,a)") " @-RFO negative Hessian eigenvalue",
     +        Alam," found: following the mode."
      Elseif (Nit .Eq. 0) Then
C Initial model Hessian: 4|f(q,q)-f(p,p)| for ov pairs (close to the
C exact pCCD values), a small constant for oo and vv pairs whose
C curvature is purely correlation (orders of magnitude below 4*de).
C (Tried 10/2026: the occupation-difference diagonal
C 2|D(p,p)-D(q,q)||f(q,q)-f(p,p)| fits the H2O STO-3G FD diagonal
C better but makes the stretched-N2 runs (1.5-3.0 A) fail.)
         If (Lhes .Eq. 1) Then
         Call Getrec(20,"JOBARC","FOCKA   ",Nsq*Iintfp,F)
         Call Dgemm("N","N",Nbas,Nbas,Nbas,One,F,Nbas,C,Nbas,Zero,
     +              Scr,Nbas)
         Call Dzero(H,Npar*Npar)
         K = 0
         Do Ip = 1, Nbas
            Do Iq = Ip+1, Nbas
               If (Irp(Ip) .Eq. Irp(Iq)) Then
               K = K + 1
               If (Ip .Le. Nocc .And. Iq .Gt. Nocc) Then
                  Fpp = Ddot(Nbas,C(1,Ip),1,Scr(1,Ip),1)
                  Fqq = Ddot(Nbas,C(1,Iq),1,Scr(1,Iq),1)
                  H(K,K) = Four*Dabs(Fqq-Fpp)
               Else
                  H(K,K) = Hoovv
               Endif
               Endif
            Enddo
         Enddo
         Endif
         If (Ifrz .Eq. 1) Then
            If (Lhes .Eq. 1) Call Dcopy(Npar*Npar,Hf,1,H,1)
            If (Ih0 .Eq. 2) Call Putrec(20,"JOBARC","RFO_HFRZ",
     +                                  Npar*Npar*Iintfp,Hf)
            Write(6,"(2a)") " @-RFO starting Hessian: analytic",
     +                      " frozen-amplitude Hessian."
         Elseif (Ih0 .Ge. 1) Then
            Write(6,"(2a)") " @-RFO no analytic frozen Hessian (pCCD",
     +           " only, memory); constant starting Hessian."
         Endif
         Trust = Trust0
         Iacpt = 1
      Else
         If (Iexp .Eq. 0) Call Getrec(20,"JOBARC","RFO_HESS",
     +                                Npar*Npar*Iintfp,H)
         Call Getrec(20,"JOBARC","RFO_GPRV",Npar*Iintfp,Gp)
         Call Getrec(20,"JOBARC","RFO_SPRV",Npar*Iintfp,Sp)
         Call Getrec(20,"JOBARC","RFO_EPRV",Iintfp,Eprv)
         Call Getrec(20,"JOBARC","RFO_PRED",Iintfp,Pred)
         Call Getrec(20,"JOBARC","RFO_TRST",Iintfp,Trust)
         Call Getrec(20,"JOBARC","RFO_CPRV",Nsq*Iintfp,Cprv)
         Snrm  = Dsqrt(Ddot(Npar,Sp,1,Sp,1))
         Dele  = Zero
         If (Ifail .Eq. 0) Then
            Dele = E - Eprv
            If (Pred .Ne. Zero) Rho = Dele/Pred
            Jump = (Dele .Lt. -Ejmp) .And. (Rho .Gt. Rhojmp)
         Endif
C Ih0 = 2: the model at an accepted point is its frozen Hessian plus
C the correction carried so far, H - Hf(old) + Hf(new) (Ah: scratch).
         Accept = .Not. (Dele .Gt. Etol .Or. Ifail .Eq. 1 .Or. Jump)
         If (Ih0 .Eq. 2 .And. Accept .And. Ifrz .Eq. 1 .And.
     +       Iexp .Eq. 0) Then
            Call Getrec(20,"JOBARC","RFO_HFRZ",Npar*Npar*Iintfp,Ah)
            Call Daxpy(Npar*Npar,-One,Ah,1,H,1)
            Call Daxpy(Npar*Npar,One,Hf,1,H,1)
            Call Putrec(20,"JOBARC","RFO_HFRZ",Npar*Npar*Iintfp,Hf)
         Endif
C SR1 update with s = previous step, y = g(new) - g(old); Step is
C used as scratch for r = y - H*s. Not with the gradient of a failed
C or branch-jumped point, nor of the exact Hessian (Iexp = 1).
         If (Ifail .Eq. 0 .And. .Not. Jump .And. Iexp .Eq. 0) Then
            Do K = 1, Npar
               Step(K) = Gv(K) - Gp(K)
            Enddo
            Call Dgemv("N",Npar,Npar,-One,H,Npar,Sp,1,One,Step,1)
            Den   = Ddot(Npar,Step,1,Sp,1)
            Rnrm  = Dsqrt(Ddot(Npar,Step,1,Step,1))
            If (Dabs(Den) .Gt. 1.0D-8*Rnrm*Snrm) Call Dger(Npar,Npar,
     +          One/Den,Step,1,Step,1,H,Npar)
         Endif
C Etol sits above the energy noise of the amplitude solution.
         If (Dele .Gt. Etol .Or. Ifail .Eq. 1 .Or. Jump) Then
            Reject = .True.
            Trust  = Half*Min(Trust,Snrm)
            Call Dcopy(Nsq,Cprv,1,C,1)
            Call Dcopy(Npar,Gp,1,Gv,1)
            E = Eprv
         Else
            Iacpt = 1
            If (Rho .Gt. 0.75D0 .And. Snrm .Gt. 0.8D0*Trust) Then
               Trust = Min(Two*Trust,Trstmx)
            Elseif (Rho .Lt. 0.25D0) Then
               Trust = Half*Trust
            Endif
         Endif
      Endif

C Ih0 = 3: the exact Hessian (frozen part + amplitude response,
C Pccd_rsp) at every accepted point: dense (Idav = 0, Npar products)
C or applied within the Davidson iterations of the step (Idav = 1).
C Iexa = 1 when in use; where its checks fail (Pccd_rsp_init) the
C frozen Hessian of the point alone, Iexa = 2. Either way a rejected
C step is shortened along its direction (the RFO direction does not
C depend on the trust radius) and no model Hessian is kept.
C Without the frozen Hessian (memory) the SR1 model, Iexa = 0.
      Iexa = Iexp
      Irsp = 0
      Hok  = Lhes .Eq. 1 .And. Iexa .Eq. 0
      If (Ih0 .Eq. 3 .And. Iacpt .Eq. 1 .And. Iphase .Ne. 1 .And.
     +    Ifrz .Eq. 1) Then
         Iexa = 2
         If (Lrsp .Gt. 0)
     +      Call Pccd_rsp_init(Wr,Lrsp,Nbas,Nocc,Npar,Irp,Gv,Irsp,0)
         If (Irsp .Eq. 1) Then
            Iexa = 1
            If (Idav .Eq. 0) Then
               Call Dcopy(Npar*Npar,Hf,1,H,1)
               Do K = 1, Npar
                  Call Dzero(Step,Npar)
                  Step(K) = One
                  Call Pccd_rsp_rlx(Wr,Npar,Step,H(1,K))
               Enddo
               Call Pccd_rsp_sym(H,Npar)
            Endif
         Else
            Write(6,"(2a)") " @-RFO exact Hessian not available at",
     +           " this point; frozen Hessian step."
            If (Idav .Eq. 0) Call Dcopy(Npar*Npar,Hf,1,H,1)
         Endif
         Hok = Idav .Eq. 0
      Endif

      Call Putrec(20,"JOBARC","RFO_ACPT",1,Iacpt)
      Gmax = Zero
      Grms = Zero
      Do K = 1, Npar
         Gmax = Max(Gmax,Dabs(Gv(K)))
         Grms = Grms + Gv(K)*Gv(K)
      Enddo
      Grms = Dsqrt(Grms/Dble(Npar))

C Number of negative model-Hessian eigenvalues (diagnostic; -1 when
C not computed, Idav = 1).
      Nneg = -1
      If (Idav .Eq. 0 .And. Hok) Then
         Call Dcopy(Npar*Npar,H,1,Ah,1)
         Call Dsyev("N","U",Npar,Ah,Npar,Ev,Wk,Lwork,Info)
         Nneg = 0
         Do K = 1, Npar
            If (Ev(K) .Lt. Zero) Nneg = Nneg + 1
         Enddo
      Endif

      If (Ifail .Eq. 1) Then
         Write(6,"(a,i4,2a)") " @-RFO it",Nit," step rejected:",
     +      " pCCD amplitudes unconverged or unphysical."
      Elseif (Jump) Then
         Write(6,"(a,i4,a,1p,E10.3,a,E10.3,0p,a,F9.2)")
     +      " @-RFO it",Nit," step rejected (amplitude branch jump):"//
     +      " dE =",Dele,"  pred =",Pred,"  rho =",Rho
      Elseif (Reject) Then
         Write(6,"(a,i4,a,F18.12,a,1p,E10.3,a,E10.3)")
     +      " @-RFO it",Nit," step rejected: E =",E+Dele,
     +      "  dE =",Dele,"  pred =",Pred
      Endif
      Write(6,"(a,i4,a,F18.12,1p,a,E10.3,a,E10.3,a,E9.2,0p,a,F6.2,
     +          a,i4)")
     +      " @-RFO it",Nit,"  E =",E,"  max|g| =",Gmax,"  rms =",Grms,
     +      "  trust =",Trust,"  rho =",Rho,"  #neg(H) =",Nneg

      If (Reject .And. Trust .Lt. Trmin) Then
         Write(6,"(2a,F18.12)") " @-RFO trust radius collapsed: every",
     +        " step is rejected; stopping at E =",E
         Call Getrec(0,"JOBARC","RFO_T2MX",Length,Junk)
         T2mx = Zero
         If (Length .Gt. 0) Call Getrec(20,"JOBARC","RFO_T2MX",Iintfp,
     +                                  T2mx)
         If (T2mx .Gt. One) Write(6,"(a,F8.3,a,1p,E10.3,2a)")
     +      " @-RFO max|T2| =",T2mx," and max|g| =",Gmax,": probably",
     +      " the unphysical OO-pCCD branch (seen for stretched N2)."
         Call Getrec(0,"JOBARC","RFO_STPI",Length,Junk)
         If (Length .Gt. 0) Then
            Call Getrec(20,"JOBARC","RFO_STPI",2*Iintfp,Stp)
            Write(6,"(a,F18.12,a,1p,E12.4)") " @-RFO last stationary"//
     +         " point checked: E =",Stp(1),", lowest Hessian"//
     +         " eigenvalue",Stp(2)
         Endif
         Write(6,"(2a)") " @-RFO the orbital optimization did not",
     +        " converge (see the rejection reasons above)."
         Call Errex
      Endif
C Not after a rejected step along a negative mode (back at the
C stationary point that was already checked).
      Grdcnv = (Gmax .Lt. Tol) .And. .Not. Escape .And.
     +         .Not. (Reject .And. Iescp .Eq. 1)
      If (Grdcnv) Then
         Call Getrec(0,"JOBARC","RFO_NSTB",Length,Junk)
         Nstb = 0
         If (Length .Gt. 0) Call Getrec(20,"JOBARC","RFO_NSTB",1,Nstb)
         If (Nstb .Lt. Mxstb) Then
C Stationary point: start the stability (lowest eigenvalue) check.
            Write(6,"(a,i5,a)") " @-RFO gradient converged after",Ntot,
     +            " gradient evaluations; checking the Hessian."
            Call Putrec(20,"JOBARC","RFO_NSTB",1,Nstb+1)
C Exact Hessian (frozen Hessian of this point + amplitude response,
C Pccd_rsp): no displaced points. Iact = -1: finite differences below.
            Iact = -1
            If (Ifrz .Eq. 1 .And. Lrsp .Gt. 0)
     +         Call Pccd_rsp_stab(Hf,Gv,Wr,Lrsp,Sb,Wm,Irp,Nbas,Nocc,
     +                            Npar,Iact,Alam,Step,1-Irsp)
C No finite-difference check without the model Hessian arrays.
            If (Iact .Eq. -1 .And. Lhes .Eq. 0) Then
               Call Putrec(20,"JOBARC","SCFEVCA0",Nsq*Iintfp,C)
               Call Putrec(20,"JOBARC","SCFEVCB0",Nsq*Iintfp,C)
               Call Putrec(20,"JOBARC","GRADTEST",1,1)
               Write(6,"(2a,i5,a)") " @-RFO orbital optimization",
     +             " converged in",Ntot," gradients; no Hessian check"//
     +             " (exact Hessian not available)."
               Return
            Endif
            If (Iact .Ge. 1) Then
               Stp(1) = E
               Stp(2) = Alam
               Call Putrec(20,"JOBARC","RFO_STPI",2*Iintfp,Stp)
            Endif
            If (Iact .Eq. 1) Then
C Minimum: probe its softest modes first (Sb, Wm: the eigenpairs).
               Call Pccd_rfo_prb("START",Ipr,C,E,Gv,Sb,Wm,Min(8,Npar),
     +                           Step,F,U,Scr,G,Irp,Nbas,Npar,0)
               If (Ipr .Eq. 0) Return
               Call Putrec(20,"JOBARC","SCFEVCA0",Nsq*Iintfp,C)
               Call Putrec(20,"JOBARC","SCFEVCB0",Nsq*Iintfp,C)
               Call Putrec(20,"JOBARC","GRADTEST",1,1)
               Write(6,"(2a,i5,a)") " @-RFO orbital optimization",
     +             " converged to a minimum in",Ntot," gradients."
               Return
            Elseif (Iact .Eq. 2) Then
C Negative exact eigenvalue: into the model along the mode (SR1 model
C only), then a step along the mode (below).
               If (Iexa .Eq. 0) Then
               Call Dgemv("N",Npar,Npar,One,H,Npar,Step,1,Zero,Gp,1)
               Vhv = Ddot(Npar,Step,1,Gp,1)
               Call Dger(Npar,Npar,Alam-Vhv,Step,1,Step,1,H,Npar)
               Endif
               Trust  = Trust0
               Escape = .True.
               Grdcnv = .False.
               Write(6,"(a,1p,E12.4,a)") " @-RFO negative Hessian"//
     +              " eigenvalue",Alam," found: following the mode."
               Go to 400
            Endif
C Davidson start vectors (after a random one): the lowest eigenvectors
C of the SR1 model Hessian, which carries the curvature seen along the
C path (unit vectors of its lowest diagonal elements can miss a
C negative mode entirely: N2 cc-pVDZ 1.1016 A, symmetric saddle
C accepted as a minimum).
            If (.Not. Hok) Call Dcopy(Npar*Npar,Hf,1,H,1)
            Ng = Min(Ngvec,Npar)
            Iconv = 0
            If (Idav .Eq. 1) Then
C Sb is free until Pccd_rfo_stab START; Ah is the workspace.
C Tight: approximate vectors mix the exact invariance zero modes with
C the soft physical modes and the stability check then cannot resolve
C them (N2 cc-pVDZ sym off: 1e-6 -> third zero mode unresolved after
C 200 Hv, INCONCLUSIVE; 1e-10: 60 Hv here).
               Call Pccd_rfo_dav(H,Gv,Wr,Npar,0,Ng,Sb,Ev,1.0D-10,Ah,
     +                           Lah,Iconv,Nmv)
               If (Iconv .Eq. 1) Then
                  Call Putrec(20,"JOBARC","RFO_GVEC",Npar*Ng*Iintfp,Sb)
               Else
                  Write(6,"(2a,i6)") " @-RFO Davidson (model Hessian",
     +                 " vectors) not converged; dense, Hv =",Nmv
               Endif
            Endif
            If (Iconv .Ne. 1) Then
            Call Dcopy(Npar*Npar,H,1,Ah,1)
            Call Dsyev("V","U",Npar,Ah,Npar,Ev,Wk,Lwork,Info)
            Call Putrec(20,"JOBARC","RFO_GVEC",Npar*Ng*Iintfp,Ah)
            Endif
            Call Pccd_rfo_pdia(Gp,C,F,Scr,Wk,Lwork,Irp,Nbas,Nocc,Npar)
            Call Putrec(20,"JOBARC","RFO_PDIA",Npar*Iintfp,Gp)
            Call Putrec(20,"JOBARC","RFO_HESS",Npar*Npar*Iintfp,H)
            Call Putrec(20,"JOBARC","RFO_TRST",Iintfp,Trust)
            Call Putrec(20,"JOBARC","RFO_NITR",1,Nit+1)
            Call Pccd_rfo_stab("START",Iact,Alam,C,Gv,E,H,Vb,Sb,Am,Wm,
     +                         Cs,Gs,Step,Res,F,U,Scr,G,Wk,Lwork,Irp,
     +                         Nbas,Npar,Mxhv)
            Grdcnv = .False.
            Call Putrec(20,"JOBARC","GRADTEST",1,0)
            Return
         Endif
         Call Putrec(20,"JOBARC","SCFEVCA0",Nsq*Iintfp,C)
         Call Putrec(20,"JOBARC","SCFEVCB0",Nsq*Iintfp,C)
         Call Putrec(20,"JOBARC","GRADTEST",1,1)
         Write(6,"(2a,i5,a)") " @-RFO orbital optimization converged",
     +         " (stability check limit reached) in",Ntot," gradients."
         Return
      Endif
 400  Continue
      Call Putrec(20,"JOBARC","GRADTEST",1,0)

      If (Escape) Then
C Along the negative mode found by the stability check (in Step), length
C Trust, downhill to first order. Not the RFO eigenvector: at a
C stationary point it is undefined when g is orthogonal to the mode
C (a symmetry-breaking mode at a symmetric saddle point).
         Vn  = Dsqrt(Ddot(Npar,Step,1,Step,1))
         Gsv = Ddot(Npar,Gv,1,Step,1)/Vn
         Fac = Trust/Vn
         If (Gsv .Gt. Zero) Fac = -Fac
         Call Dscal(Npar,Fac,Step,1)
         Snrm   = Trust
         Gsh(1) = -Dabs(Gsv)*Trust
         Gsh(2) = Alam*Trust*Trust
         Go to 500
      Endif
      If (Reject .And. (Iexp .Ge. 1 .Or. Iescp .Eq. 1)) Then
C The rejected step, shortened to the new trust radius (Snrm = |Sp|).
         Call Getrec(20,"JOBARC","RFO_GSHS",2*Iintfp,Gsh)
         Fac = Trust/Snrm
         Do K = 1, Npar
            Step(K) = Fac*Sp(K)
         Enddo
         Snrm   = Trust
         Gsh(1) = Fac*Gsh(1)
         Gsh(2) = Fac*Fac*Gsh(2)
         Go to 500
      Endif

C RFO step: lowest eigenvector of [[H,g],[g^T,0]], scaled to the trust
C radius when longer. Idav = 1: Davidson (eigenvector in Ev, Ah the
C workspace), dense if it does not converge; with the exact Hessian
C (Iexa = 1) its products come from Pccd_rsp_hv.
      Iconv = 0
      If (Idav .Eq. 1 .And. Iexa .Eq. 1) Then
C Residual tolerance relative to |g| (an absolute 1e-9 needed 600+
C products near convergence, N2 cc-pVTZ SYMMETRY=OFF, |g| 1e-5). If
C not converged, the last Ritz vector (the trust region judges it).
         Rtol = Max(1.0D-12,1.0D-4*Dsqrt(Ddot(Npar,Gv,1,Gv,1)))
         Call Pccd_rfo_dav(Hf,Gv,Wr,Npar,3,1,Ev,Wm,Rtol,Ah,Lah,
     +                     Iconv,Nmv)
         If (Iconv .Eq. 0 .And. Nmv .Gt. 0) Iconv = 2
         If (Iconv .Ge. 1 .And. Ev(Np1) .Ne. Zero) Then
            Do K = 1, Npar
               Step(K) = Ev(K)/Ev(Np1)
            Enddo
            If (Ddot(Npar,Gv,1,Step,1) .Gt. Zero)
     +         Call Dscal(Npar,-One,Step,1)
            If (Iconv .Eq. 1) Then
            Write(6,"(a,i5,a)") " @-RFO exact-Hessian step from"//
     +                          " Davidson (",Nmv," Hv)."
            Else
            Write(6,"(a,i5,a)") " @-RFO exact-Hessian step from"//
     +           " Davidson (",Nmv," Hv, not converged: Ritz vector)."
            Endif
            Iconv = 1
         Else
            Iconv = 0
            Iexa  = 2
            Write(6,"(2a,i6)") " @-RFO exact-Hessian Davidson (step)",
     +           " not converged; frozen Hessian, Hv =",Nmv
         Endif
      Endif
      If (Idav .Eq. 1 .And. Iconv .Ne. 1) Then
C Model (Iexa = 0) or frozen (Iexa = 2) Hessian matrix.
         If (Iexa .Eq. 0) Then
            Call Pccd_rfo_dav(H,Gv,Wr,Npar,1,1,Ev,Wm,1.0D-10,Ah,Lah,
     +                        Iconv,Nmv)
         Else
            Call Pccd_rfo_dav(Hf,Gv,Wr,Npar,1,1,Ev,Wm,1.0D-10,Ah,Lah,
     +                        Iconv,Nmv)
         Endif
         If (Iconv .Eq. 1 .And. Ev(Np1) .Ne. Zero) Then
            Do K = 1, Npar
               Step(K) = Ev(K)/Ev(Np1)
            Enddo
            Write(6,"(a,i5,a)") " @-RFO step from Davidson (",Nmv,
     +                          " Hv)."
         Else
            Iconv = 0
            Write(6,"(2a,i6)") " @-RFO Davidson (step) not converged;",
     +                         " dense, Hv =",Nmv
         Endif
      Endif
      If (Iconv .Ne. 1) Then
      If (Lah .Lt. Np1*Np1) Then
         Write(6,"(2a)") " @-Pccd_rfo_step: Davidson failed and no",
     +        " memory for the dense RFO step."
         Call Errex
      Endif
      If (Iexa .Eq. 2) Call Dcopy(Npar*Npar,Hf,1,H,1)
      Call Dzero(Ah,Np1*Np1)
      Do J = 1, Npar
         Do I = 1, Npar
            Ah(I,J) = H(I,J)
         Enddo
         Ah(Np1,J) = Gv(J)
         Ah(J,Np1) = Gv(J)
      Enddo
      Call Dsyev("V","U",Np1,Ah,Np1,Ev,Wk,Lwork,Info)
      If (Info .Ne. 0 .Or. Ah(Np1,1) .Eq. Zero) Then
         Write(6,"(a,i6)") " @-Pccd_rfo_step: DSYEV failed, info =",Info
         Call Errex
      Endif
      Do K = 1, Npar
         Step(K) = Ah(K,1)/Ah(Np1,1)
      Enddo
      Endif
      Snrm = Dsqrt(Ddot(Npar,Step,1,Step,1))
      If (Snrm .Gt. Trust) Then
         Call Dscal(Npar,Trust/Snrm,Step,1)
         Snrm = Trust
      Endif
      If (Iexa .Eq. 1 .And. Idav .Eq. 1) Then
         Call Pccd_rsp_hv(Hf,Wr,Npar,Step,Gp)
      Elseif (Iexa .Eq. 2 .And. Idav .Eq. 1) Then
         Call Dgemv("N",Npar,Npar,One,Hf,Npar,Step,1,Zero,Gp,1)
      Else
         Call Dgemv("N",Npar,Npar,One,H,Npar,Step,1,Zero,Gp,1)
      Endif
      Gsh(1) = Ddot(Npar,Gv,1,Step,1)
      Gsh(2) = Ddot(Npar,Step,1,Gp,1)
 500  Continue
      Call Dcopy(Npar,Step,1,Sp,1)
      Pred = Gsh(1) + Half*Gsh(2)

C Save the state at the current (accepted) point.
      Call Putrec(20,"JOBARC","RFO_CPRV",Nsq*Iintfp,C)
      Call Putrec(20,"JOBARC","RFO_GPRV",Npar*Iintfp,Gv)
      Call Putrec(20,"JOBARC","RFO_SPRV",Npar*Iintfp,Sp)
      If (Iexa .Eq. 0) Call Putrec(20,"JOBARC","RFO_HESS",
     +                             Npar*Npar*Iintfp,H)
      Call Putrec(20,"JOBARC","RFO_HEXA",1,Iexa)
      Iesc = 0
      If (Escape .Or. (Reject .And. Iescp .Eq. 1)) Iesc = 1
      Call Putrec(20,"JOBARC","RFO_ESCP",1,Iesc)
      Call Putrec(20,"JOBARC","RFO_GSHS",2*Iintfp,Gsh)
      Call Putrec(20,"JOBARC","RFO_EPRV",Iintfp,E)
      Call Putrec(20,"JOBARC","RFO_PRED",Iintfp,Pred)
      Call Putrec(20,"JOBARC","RFO_TRST",Iintfp,Trust)
      Call Putrec(20,"JOBARC","RFO_NITR",1,Nit+1)

C C(new) = C*exp(K(step)).
      Call Pccd_rfo_disp(C,Step,One,F,U,Scr,G,Irp,Nbas,Npar)
      Write(6,"(a,1p,E10.3,a,E11.3)") " @-RFO step |s| =",Snrm,
     +                                "  predicted dE =",Pred
      Return
      End


