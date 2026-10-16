      Subroutine Pccd_rfo_do(G,C,Cprv,F,U,Scr,Gv,Gp,Sp,Step,H,Ah,Ev,
     +                       Wk,Lwork,Nbas,Nocc,Npar,Grdcnv,Tol,Gmax,
     +                       Grms,Ncycle)
      Implicit Double Precision(A-H,O-Z)
      Logical Grdcnv,Reject
      Dimension G(Nbas,Nbas),C(Nbas,Nbas),Cprv(Nbas,Nbas)
      Dimension F(Nbas,Nbas),U(Nbas,Nbas),Scr(Nbas,Nbas)
      Dimension Gv(Npar),Gp(Npar),Sp(Npar),Step(Npar)
      Dimension H(Npar,Npar),Ah(Npar+1,Npar+1),Ev(Npar+1),Wk(Lwork)
#include "machsp.com"
      Parameter (Hoovv=0.05D0, Trust0=0.1D0, Trstmx=0.5D0)
      Data Zero,Half,One,Two,Four /0.0D0,0.5D0,1.0D0,2.0D0,4.0D0/

      Nsq = Nbas*Nbas
      Np1 = Npar + 1
      Call Getrec(20,"JOBARC","OBRTGRDA",Nsq*Iintfp,G)
      Call Getrec(20,"JOBARC","TOTENERG",Iintfp,E)
      Call Getrec(20,"JOBARC","SCFEVCA0",Nsq*Iintfp,C)

      K = 0
      Do Ip = 1, Nbas
         Do Iq = Ip+1, Nbas
            K = K + 1
            Gv(K) = -Two*(G(Ip,Iq) - G(Iq,Ip))
         Enddo
      Enddo

      Call Getrec(0,"JOBARC","RFO_NITR",Length,Junk)
      Nit = 0
      If (Length .Gt. 0) Call Getrec(20,"JOBARC","RFO_NITR",1,Nit)
      Reject = .False.
      Rho    = Zero

      If (Nit .Eq. 0) Then
C Initial model Hessian: 4|f(q,q)-f(p,p)| for ov pairs (close to the
C exact pCCD values), a small constant for oo and vv pairs whose
C curvature is purely correlation (orders of magnitude below 4*de).
         Call Getrec(20,"JOBARC","FOCKA   ",Nsq*Iintfp,F)
         Call Dgemm("N","N",Nbas,Nbas,Nbas,One,F,Nbas,C,Nbas,Zero,
     +              Scr,Nbas)
         Call Dzero(H,Npar*Npar)
         K = 0
         Do Ip = 1, Nbas
            Do Iq = Ip+1, Nbas
               K = K + 1
               If (Ip .Le. Nocc .And. Iq .Gt. Nocc) Then
                  Fpp = Ddot(Nbas,C(1,Ip),1,Scr(1,Ip),1)
                  Fqq = Ddot(Nbas,C(1,Iq),1,Scr(1,Iq),1)
                  H(K,K) = Four*Dabs(Fqq-Fpp)
               Else
                  H(K,K) = Hoovv
               Endif
            Enddo
         Enddo
         Trust = Trust0
      Else
         Call Getrec(20,"JOBARC","RFO_HESS",Npar*Npar*Iintfp,H)
         Call Getrec(20,"JOBARC","RFO_GPRV",Npar*Iintfp,Gp)
         Call Getrec(20,"JOBARC","RFO_SPRV",Npar*Iintfp,Sp)
         Call Getrec(20,"JOBARC","RFO_EPRV",Iintfp,Eprv)
         Call Getrec(20,"JOBARC","RFO_PRED",Iintfp,Pred)
         Call Getrec(20,"JOBARC","RFO_TRST",Iintfp,Trust)
         Call Getrec(20,"JOBARC","RFO_CPRV",Nsq*Iintfp,Cprv)
C SR1 update with s = previous step, y = g(new) - g(old); Step is
C used as scratch for r = y - H*s.
         Do K = 1, Npar
            Step(K) = Gv(K) - Gp(K)
         Enddo
         Call Dgemv("N",Npar,Npar,-One,H,Npar,Sp,1,One,Step,1)
         Den   = Ddot(Npar,Step,1,Sp,1)
         Rnrm  = Dsqrt(Ddot(Npar,Step,1,Step,1))
         Snrm  = Dsqrt(Ddot(Npar,Sp,1,Sp,1))
         If (Dabs(Den) .Gt. 1.0D-8*Rnrm*Snrm) Call Dger(Npar,Npar,
     +       One/Den,Step,1,Step,1,H,Npar)
         Dele = E - Eprv
         If (Pred .Ne. Zero) Rho = Dele/Pred
         If (Dele .Gt. 1.0D-11) Then
            Reject = .True.
            Trust  = Half*Min(Trust,Snrm)
            Call Dcopy(Nsq,Cprv,1,C,1)
            Call Dcopy(Npar,Gp,1,Gv,1)
            E = Eprv
         Else
            If (Rho .Gt. 0.75D0 .And. Snrm .Gt. 0.8D0*Trust) Then
               Trust = Min(Two*Trust,Trstmx)
            Elseif (Rho .Lt. 0.25D0) Then
               Trust = Half*Trust
            Endif
         Endif
      Endif

      Gmax = Zero
      Grms = Zero
      Do K = 1, Npar
         Gmax = Max(Gmax,Dabs(Gv(K)))
         Grms = Grms + Gv(K)*Gv(K)
      Enddo
      Grms = Dsqrt(Grms/Dble(Npar))

C Number of negative model-Hessian eigenvalues (diagnostic).
      Call Dcopy(Npar*Npar,H,1,Ah,1)
      Call Dsyev("N","U",Npar,Ah,Npar,Ev,Wk,Lwork,Info)
      Nneg = 0
      Do K = 1, Npar
         If (Ev(K) .Lt. Zero) Nneg = Nneg + 1
      Enddo

      If (Reject) Then
         Write(6,"(a,i4,a,F18.12,a,1p,E10.3,a,E10.3)")
     +      " @-RFO it",Nit," step rejected: E =",E+Dele,
     +      "  dE =",Dele,"  pred =",Pred
      Endif
      Write(6,"(a,i4,a,F18.12,1p,a,E10.3,a,E10.3,a,E9.2,0p,a,F6.2,
     +          a,i4)")
     +      " @-RFO it",Nit,"  E =",E,"  max|g| =",Gmax,"  rms =",Grms,
     +      "  trust =",Trust,"  rho =",Rho,"  #neg(H) =",Nneg

      Grdcnv = (Gmax .Lt. Tol)
      Ncycle = Nit + 1
      If (Grdcnv) Then
         Call Putrec(20,"JOBARC","SCFEVCA0",Nsq*Iintfp,C)
         Call Putrec(20,"JOBARC","SCFEVCB0",Nsq*Iintfp,C)
         Call Putrec(20,"JOBARC","GRADTEST",1,1)
         Write(6,"(a,i4,a)") " @-RFO orbital optimization converged in",
     +                       Nit+1," gradient evaluations."
         Return
      Endif
      Call Putrec(20,"JOBARC","GRADTEST",1,0)

C RFO step: lowest eigenvector of [[H,g],[g^T,0]], scaled to the trust
C radius when longer.
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
      Snrm = Dsqrt(Ddot(Npar,Step,1,Step,1))
      If (Snrm .Gt. Trust) Then
         Call Dscal(Npar,Trust/Snrm,Step,1)
         Snrm = Trust
      Endif
      Call Dcopy(Npar,Step,1,Sp,1)
      Call Dgemv("N",Npar,Npar,One,H,Npar,Step,1,Zero,Gp,1)
      Pred = Ddot(Npar,Gv,1,Step,1) + Half*Ddot(Npar,Step,1,Gp,1)

C Save the state at the current (accepted) point.
      Call Putrec(20,"JOBARC","RFO_CPRV",Nsq*Iintfp,C)
      Call Putrec(20,"JOBARC","RFO_GPRV",Npar*Iintfp,Gv)
      Call Putrec(20,"JOBARC","RFO_SPRV",Npar*Iintfp,Sp)
      Call Putrec(20,"JOBARC","RFO_HESS",Npar*Npar*Iintfp,H)
      Call Putrec(20,"JOBARC","RFO_EPRV",Iintfp,E)
      Call Putrec(20,"JOBARC","RFO_PRED",Iintfp,Pred)
      Call Putrec(20,"JOBARC","RFO_TRST",Iintfp,Trust)
      Call Putrec(20,"JOBARC","RFO_NITR",1,Nit+1)

C U = exp(K) by scaling and squaring of the Taylor series, then
C C(new) = C*U, re-orthonormalized (modified Gram-Schmidt in the
C MO metric, i.e. on the columns of U).
      Call Dzero(F,Nsq)
      K = 0
      Do Ip = 1, Nbas
         Do Iq = Ip+1, Nbas
            K = K + 1
            F(Ip,Iq) =  Step(K)
            F(Iq,Ip) = -Step(K)
         Enddo
      Enddo
      Call Pccd_rfo_expm(F,U,Scr,G,Nbas)
      Call Dgemm("N","N",Nbas,Nbas,Nbas,One,C,Nbas,U,Nbas,Zero,
     +           Scr,Nbas)
      Call Putrec(20,"JOBARC","SCFEVCA0",Nsq*Iintfp,Scr)
      Call Putrec(20,"JOBARC","SCFEVCB0",Nsq*Iintfp,Scr)
      Write(6,"(a,1p,E10.3,a,E11.3)") " @-RFO step |s| =",Snrm,
     +                                "  predicted dE =",Pred
      Return
      End

