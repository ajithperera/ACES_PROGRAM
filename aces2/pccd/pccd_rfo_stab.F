      Subroutine Pccd_rfo_stab(Mode,Iact,Alam,C,Gv,E,H,Vb,Sb,Am,Wm,Cs,
     +                         Gs,Vec,Res,F,U,Scr,G,Wk,Lwork,Nbas,Npar,
     +                         Mxhv)
C
C Lowest eigenpair of the exact orbital Hessian at a stationary point
C (Cs, Gs, Es) by Davidson iterations with forward-difference
C Hessian-vector products of the analytic gradient,
C     H*b = (g(Cs*exp(Delta*b)) - g(Cs))/Delta ,
C one per xpccd run (macro iteration). Mode:
C   "START": store the stationary point, displace along the first
C            (pseudo-random, all-symmetry) trial vector.        Iact=0
C   "CONT ": Gv is the gradient at the last displaced point. Either
C            displace along a new trial vector (Iact=0), or finish:
C            Iact=1 stable (lowest eigenvalue > Xlamth), Iact=2 a
C            negative eigenvalue Alam with eigenvector Vec was found.
C            For Iact=1,2 C, Gv and E are reset to the stationary point.
C Rtol=5e-3 sits above the forward-difference noise floor of the Hv
C products (~1e-3 for Delta=1e-3). The lowest three Ritz pairs are
C tracked so that a nearly converged
C positive root cannot mask a negative one that enters late.
C
      Implicit Double Precision(A-H,O-Z)
      Character*5 Mode
      Logical Ok
      Dimension C(Nbas*Nbas),Gv(Npar),H(Npar,Npar)
      Dimension Vb(Npar,Mxhv),Sb(Npar,Mxhv),Am(Mxhv,Mxhv),Wm(Mxhv)
      Dimension Cs(Nbas*Nbas),Gs(Npar),Vec(Npar),Res(Npar)
      Dimension F(Nbas,Nbas),U(Nbas,Nbas),Scr(Nbas,Nbas),G(Nbas,Nbas)
      Dimension Wk(Lwork),Rn(3)
#include "machsp.com"
      Parameter (Delta=1.0D-3, Rtol=5.0D-3, Xnegst=-1.0D-3,
     +           Xlamth=-1.0D-4, Nguess=3, Mxroot=3)
      Data Zero,Half,One /0.0D0,0.5D0,1.0D0/

      Nsq  = Nbas*Nbas
      Iact = 0
      Alam = Zero

      If (Mode .Eq. "START") Then
         Call Dcopy(Nsq,C,1,Cs,1)
         Call Dcopy(Npar,Gv,1,Gs,1)
         Es  = E
         Nhv = 0
         Do K = 1, Npar
            X = Dsin(Dble(K)*12.9898D0+0.3D0)*43758.5453D0
            Vec(K) = X - Aint(X) - Half*Dsign(One,X)
         Enddo
         Call Pccd_rfo_addvec(Vec,Vb,Npar,Mxhv,Nhv,Ok)
         Call Pccd_rfo_disp(Cs,Vb(1,Nhv),Delta,F,U,Scr,G,Nbas,Npar)
         Call Putrec(20,"JOBARC","RFO_CSTB",Nsq*Iintfp,Cs)
         Call Putrec(20,"JOBARC","RFO_GSTB",Npar*Iintfp,Gs)
         Call Putrec(20,"JOBARC","RFO_ESTB",Iintfp,Es)
         Call Putrec(20,"JOBARC","RFO_VB  ",Npar*Mxhv*Iintfp,Vb)
         Call Putrec(20,"JOBARC","RFO_NHV ",1,Nhv)
         Call Putrec(20,"JOBARC","RFO_PHAS",1,1)
         Return
      Endif

      Call Getrec(20,"JOBARC","RFO_CSTB",Nsq*Iintfp,Cs)
      Call Getrec(20,"JOBARC","RFO_GSTB",Npar*Iintfp,Gs)
      Call Getrec(20,"JOBARC","RFO_ESTB",Iintfp,Es)
      Call Getrec(20,"JOBARC","RFO_VB  ",Npar*Mxhv*Iintfp,Vb)
      Call Getrec(20,"JOBARC","RFO_NHV ",1,Nhv)
      If (Nhv .Gt. 1) Call Getrec(20,"JOBARC","RFO_SB  ",
     +                            Npar*(Nhv-1)*Iintfp,Sb)
      Do K = 1, Npar
         Sb(K,Nhv) = (Gv(K) - Gs(K))/Delta
      Enddo
      Call Putrec(20,"JOBARC","RFO_SB  ",Npar*Nhv*Iintfp,Sb)

      If (Nhv .Lt. Nguess) Then
C Further start vectors: the unit vectors of the lowest model-Hessian
C diagonal elements (the Nhv-th smallest, ties broken by index).
         Kprev = 0
         Do Ig = 1, Nhv
            Dmin  = 1.0D30
            Kbest = 0
            Do K = 1, Npar
               If (K .Ne. Kprev .And. H(K,K) .Lt. Dmin) Then
                  If (Ig .Eq. 1 .Or. H(K,K) .Gt. Dprev .Or.
     +               (H(K,K) .Eq. Dprev .And. K .Gt. Kprev)) Then
                     Dmin  = H(K,K)
                     Kbest = K
                  Endif
               Endif
            Enddo
            Kprev = Kbest
            Dprev = Dmin
         Enddo
         Call Dzero(Vec,Npar)
         Vec(Kbest) = One
         Call Pccd_rfo_addvec(Vec,Vb,Npar,Mxhv,Nhv,Ok)
         If (Ok) Then
            Call Pccd_rfo_disp(Cs,Vb(1,Nhv),Delta,F,U,Scr,G,Nbas,Npar)
            Call Putrec(20,"JOBARC","RFO_VB  ",Npar*Mxhv*Iintfp,Vb)
            Call Putrec(20,"JOBARC","RFO_NHV ",1,Nhv)
            Return
         Endif
      Endif

C Rayleigh-Ritz in the current subspace.
      Do J = 1, Nhv
         Do I = 1, Nhv
            Am(I,J) = Half*(Ddot(Npar,Vb(1,I),1,Sb(1,J),1) +
     +                      Ddot(Npar,Vb(1,J),1,Sb(1,I),1))
         Enddo
      Enddo
      Call Dsyev("V","U",Nhv,Am,Mxhv,Wm,Wk,Lwork,Info)
      Nroot = Min(Mxroot,Nhv)
      Iunc  = 0
      Do Ir = Nroot, 1, -1
         Call Dgemv("N",Npar,Nhv,One,Sb,Npar,Am(1,Ir),1,Zero,Res,1)
         Call Dgemv("N",Npar,Nhv,-Wm(Ir),Vb,Npar,Am(1,Ir),1,One,Res,1)
         Rn(Ir) = Dsqrt(Ddot(Npar,Res,1,Res,1))
         If (Rn(Ir) .Gt. Rtol) Iunc = Ir
      Enddo
      Alam = Wm(1)
      Call Dgemv("N",Npar,Nhv,One,Vb,Npar,Am(1,1),1,Zero,Vec,1)
      Write(6,"(a,i4,a,1p,3E12.4,a,3E9.2)") " @-RFO Davidson",Nhv,
     +      " Hv: lowest roots",(Wm(Ir),Ir=1,Nroot),"  |r|",
     +      (Rn(Ir),Ir=1,Nroot)

      If (Alam .Lt. Xnegst .And. Rn(1) .Lt. 0.3D0*Dabs(Alam)) Then
         Iact = 2
      Elseif (Iunc .Eq. 0 .Or. Nhv .Ge. Mxhv) Then
         Iact = 1
         If (Alam .Lt. Xlamth) Iact = 2
      Else
C Expand with the preconditioned residual of the lowest unconverged root.
         Call Dgemv("N",Npar,Nhv,One,Sb,Npar,Am(1,Iunc),1,Zero,Res,1)
         Call Dgemv("N",Npar,Nhv,-Wm(Iunc),Vb,Npar,Am(1,Iunc),1,One,
     +              Res,1)
         Do K = 1, Npar
            Den = H(K,K) - Wm(Iunc)
            If (Dabs(Den) .Lt. 1.0D-2) Den = Dsign(1.0D-2,Den)
            Wk(K) = -Res(K)/Den
         Enddo
         Call Pccd_rfo_addvec(Wk,Vb,Npar,Mxhv,Nhv,Ok)
         If (.Not. Ok) Call Pccd_rfo_addvec(Res,Vb,Npar,Mxhv,Nhv,Ok)
         If (Ok) Then
            Call Pccd_rfo_disp(Cs,Vb(1,Nhv),Delta,F,U,Scr,G,Nbas,Npar)
            Call Putrec(20,"JOBARC","RFO_VB  ",Npar*Mxhv*Iintfp,Vb)
            Call Putrec(20,"JOBARC","RFO_NHV ",1,Nhv)
            Return
         Endif
         Iact = 1
         If (Alam .Lt. Xlamth) Iact = 2
      Endif

      If (Iact .Eq. 1) Write(6,"(a,1p,E12.4,a)")
     +   " @-RFO lowest Hessian eigenvalue",Alam,": minimum."
      Call Dcopy(Nsq,Cs,1,C,1)
      Call Dcopy(Npar,Gs,1,Gv,1)
      E = Es
      Call Putrec(20,"JOBARC","RFO_PHAS",1,0)
      Return
      End

      Subroutine Pccd_rfo_addvec(B,Vb,Npar,Mxhv,Nhv,Ok)
C Orthonormalize B against the Nhv columns of Vb (twice) and append it.
      Implicit Double Precision(A-H,O-Z)
      Logical Ok
      Dimension B(Npar),Vb(Npar,Mxhv)
      Ok = .False.
      If (Nhv .Ge. Mxhv) Return
      Bn0 = Dsqrt(Ddot(Npar,B,1,B,1))
      If (Bn0 .Eq. 0.0D0) Return
      Do Ipass = 1, 2
         Do J = 1, Nhv
            Ov = Ddot(Npar,Vb(1,J),1,B,1)
            Call Daxpy(Npar,-Ov,Vb(1,J),1,B,1)
         Enddo
      Enddo
      Bn = Dsqrt(Ddot(Npar,B,1,B,1))
      If (Bn .Lt. 1.0D-6*Bn0) Return
      Nhv = Nhv + 1
      Do K = 1, Npar
         Vb(K,Nhv) = B(K)/Bn
      Enddo
      Ok = .True.
      Return
      End
