      Subroutine Pccd_rfo_step(Work,Maxcor,Nbas,Nocc,Iuhf,Grdcnv,Tol,
     +                         Etol,Gmax,Grms,Ncycle)
C
C OO_METHOD=RFO: trust-region rational-function (augmented Hessian)
C orbital step with an SR1-updated model Hessian over all orbital pairs
C p<q (oo, ov and vv rotations, fully coupled). One call per macro
C iteration (xvscf/xvtran/xintprc/xpccd); the state between calls is
C kept in JOBARC (RFO_* records). Uphill steps are rejected: the
C previous orbitals are restored and the trust radius is reduced.
C At every convergence the lowest eigenvalue of the exact orbital
C Hessian is found by Davidson iterations (Pccd_rfo_stab) using finite-
C difference Hessian-vector products of the analytic gradient, one per
C macro iteration. A negative eigenvalue is injected into the model
C Hessian and the optimization continues along that mode (a gradient-
C only method can never leave e.g. a symmetric saddle point).
C
C Parametrization: C(new) = C * exp(K), K(p,q) = x(pq) = -K(q,p), p<q,
C dE/dx(pq) = -2*(G(p,q)-G(q,p)) with G = OBRTGRDA (verified against
C finite-difference energies, 10/2026).
C
C RHF, SYMMETRY=OFF, no dropped MOs.
C
      Implicit Double Precision(A-H,O-Z)
      Logical Grdcnv,Symmetry
      Dimension Work(Maxcor)
      Common/Symm/Symmetry
#include "machsp.com"

      If (Iuhf .Ne. 0 .Or. Symmetry) Then
         Write(6,"(2a)") " @-Pccd_rfo_step: OO_METHOD=RFO requires",
     +                   " REF=RHF and SYMMETRY=OFF."
         Call Errex
      Endif
      Call Getrec(20,"JOBARC","NBASTOT ",1,Nao)
      If (Nao .Ne. Nbas) Then
         Write(6,"(2a)") " @-Pccd_rfo_step: number of MOs and AOs",
     +                   " differ (dropped/linear dependent MOs)."
         Call Errex
      Endif

      Npar = Nbas*(Nbas-1)/2
      Nsq  = Nbas*Nbas
      Np1  = Npar + 1
      Mxhv = Min(Npar,60)
      Lwork = Max(64*Np1,3*Np1)
      I000 = 1
      I010 = I000 + Nsq
      I020 = I010 + Nsq
      I030 = I020 + Nsq
      I040 = I030 + Nsq
      I050 = I040 + Nsq
      I060 = I050 + Nsq
      I070 = I060 + Npar
      I080 = I070 + Npar
      I090 = I080 + Npar
      I100 = I090 + Npar
      I110 = I100 + Npar*Npar
      I120 = I110 + Np1*Np1
      I130 = I120 + Np1
      I140 = I130 + Lwork
      I150 = I140 + Npar*Mxhv
      I160 = I150 + Npar*Mxhv
      I170 = I160 + Mxhv*Mxhv
      I180 = I170 + Mxhv
      I190 = I180 + Nsq
      I200 = I190 + Npar
      Iend = I200 + Npar
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_rfo_step",Iend,Maxcor)

      Call Pccd_rfo_do(Work(I000),Work(I010),Work(I020),Work(I030),
     +                 Work(I040),Work(I050),Work(I060),Work(I070),
     +                 Work(I080),Work(I090),Work(I100),Work(I110),
     +                 Work(I120),Work(I130),Lwork,Work(I140),
     +                 Work(I150),Work(I160),Work(I170),Work(I180),
     +                 Work(I190),Work(I200),Mxhv,Nbas,Nocc,Npar,
     +                 Grdcnv,Tol,Etol,Gmax,Grms,Ncycle)
C The T2 guess of the next macro iterations follows accepted points only.
      Call Getrec(20,"JOBARC","RFO_ACPT",1,Iacpt)
      If (Iacpt .Eq. 1) Call Pccd_rfo_tguess(Work,Maxcor,"PROMO")
      Return
      End

