










      Subroutine Pccd_rfo_do(G,C,Cprv,F,U,Scr,Gv,Gp,Sp,Step,H,Ah,Ev,
     +                       Wk,Lwork,Vb,Sb,Am,Wm,Cs,Gs,Res,Mxhv,
     +                       Nbas,Nocc,Npar,Grdcnv,Tol,Gmax,Grms,Ncycle)
      Implicit Double Precision(A-H,O-Z)
      Logical Grdcnv,Reject,Escape
      Dimension Vb(Npar,Mxhv),Sb(Npar,Mxhv),Am(Mxhv,Mxhv),Wm(Mxhv)
      Dimension Cs(Nbas*Nbas),Gs(Npar),Res(Npar)
      Dimension G(Nbas,Nbas),C(Nbas,Nbas),Cprv(Nbas,Nbas)
      Dimension F(Nbas,Nbas),U(Nbas,Nbas),Scr(Nbas,Nbas)
      Dimension Gv(Npar),Gp(Npar),Sp(Npar),Step(Npar)
      Dimension H(Npar,Npar),Ah(Npar+1,Npar+1),Ev(Npar+1),Wk(Lwork)


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



      Parameter (Hoovv=0.05D0, Trust0=0.1D0, Trstmx=0.5D0)
      Parameter (Mxstb=5)
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
      Rho    = Zero

      If (Iphase .Eq. 1) Then
C Stability phase: this gradient belongs to a displaced point of the
C Davidson iterations.
         Call Getrec(20,"JOBARC","RFO_HESS",Npar*Npar*Iintfp,H)
         Call Getrec(20,"JOBARC","RFO_TRST",Iintfp,Trust)
         Call Pccd_rfo_stab("CONT ",Iact,Alam,C,Gv,E,H,Vb,Sb,Am,Wm,Cs,
     +                      Gs,Step,Res,F,U,Scr,G,Wk,Lwork,Nbas,Npar,
     +                      Mxhv)
         Grdcnv = .False.
         If (Iact .Eq. 0) Then
            Call Putrec(20,"JOBARC","GRADTEST",1,0)
            Return
         Elseif (Iact .Eq. 1) Then
            Grdcnv = .True.
            Call Putrec(20,"JOBARC","SCFEVCA0",Nsq*Iintfp,C)
            Call Putrec(20,"JOBARC","SCFEVCB0",Nsq*Iintfp,C)
            Call Putrec(20,"JOBARC","GRADTEST",1,1)
            Write(6,"(2a,i5,a)") " @-RFO orbital optimization",
     +             " converged to a minimum in",Ntot," gradients."
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

      Grdcnv = (Gmax .Lt. Tol) .And. .Not. Escape
      If (Grdcnv) Then
         Call Getrec(0,"JOBARC","RFO_NSTB",Length,Junk)
         Nstb = 0
         If (Length .Gt. 0) Call Getrec(20,"JOBARC","RFO_NSTB",1,Nstb)
         If (Nstb .Lt. Mxstb) Then
C Stationary point: start the stability (lowest eigenvalue) check.
            Write(6,"(a,i5,a)") " @-RFO gradient converged after",Ntot,
     +            " gradient evaluations; checking the Hessian."
            Call Putrec(20,"JOBARC","RFO_NSTB",1,Nstb+1)
            Call Putrec(20,"JOBARC","RFO_HESS",Npar*Npar*Iintfp,H)
            Call Putrec(20,"JOBARC","RFO_TRST",Iintfp,Trust)
            Call Putrec(20,"JOBARC","RFO_NITR",1,Nit+1)
            Call Pccd_rfo_stab("START",Iact,Alam,C,Gv,E,H,Vb,Sb,Am,Wm,
     +                         Cs,Gs,Step,Res,F,U,Scr,G,Wk,Lwork,Nbas,
     +                         Npar,Mxhv)
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

C C(new) = C*exp(K(step)).
      Call Pccd_rfo_disp(C,Step,One,F,U,Scr,G,Nbas,Npar)
      Write(6,"(a,1p,E10.3,a,E11.3)") " @-RFO step |s| =",Snrm,
     +                                "  predicted dE =",Pred
      Return
      End

