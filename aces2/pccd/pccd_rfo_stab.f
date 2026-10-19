










      Subroutine Pccd_rfo_stab(Mode,Iact,Alam,C,Gv,E,H,Vb,Sb,Am,Wm,Cs,
     +                         Gs,Vec,Res,F,U,Scr,G,Wk,Lwork,Irp,Nbas,
     +                         Npar,Mxhv)
C
C Lowest eigenpair of the exact orbital Hessian at a stationary point
C (Cs, Gs, Es) by Davidson iterations with central-difference
C Hessian-vector products of the analytic gradient,
C     H*b = (g(Cs*exp(Delta*b)) - g(Cs*exp(-Delta*b)))/(2*Delta) ,
C one displaced point per xpccd run (macro iteration), two per Hv.
C RFO_HVPH = 1 (2): the next gradient belongs to the +Delta (-Delta)
C point of the newest trial vector. Mode:
C   "START": store the stationary point, displace along the first
C            (pseudo-random, all-symmetry) trial vector.        Iact=0
C   "CONT ": Gv is the gradient at the last displaced point. Either
C            displace along a new trial vector (Iact=0), or finish:
C            Iact=1 stable (lowest eigenvalue > Xlamth), Iact=2 a
C            negative eigenvalue Alam with eigenvector Vec was found.
C            For Iact=1,2 C, Gv and E are reset to the stationary point.
C Forward differences (noise ~1e-3 for Delta=1e-3, Rtol 5e-3) could
C not resolve curvatures below ~1e-3: N2 cc-pVDZ 1.1016 A, symmetric
C saddle with lowest eigenvalue -3.6e-4 (full FD Hessian) estimated as
C +1.4e-4. Central differences (Hessian asymmetry ~1e-5 there) allow
C Rtol=5e-4. The lowest three Ritz pairs are tracked so that a nearly
C converged positive root cannot mask a negative one that enters late.
C When the subspace is full (Mxhv) it is collapsed onto the lowest
C Mkeep Ritz vectors (thick restart), up to Mxtot Hv products in all;
C at the 60-product cap the N2 1.1016 A symmetric saddle was accepted
C as a minimum while the second root was still falling toward the
C exact -3.6e-4. Iact=3: budget exhausted with no eigenvalue below
C Xlamth and unconverged roots (inconclusive). The preconditioner is
C RFO_PDIA (Pccd_rfo_pdia), not the model Hessian diagonal.
C
      Implicit Double Precision(A-H,O-Z)
      Character*5 Mode
      Logical Ok
      Dimension C(Nbas*Nbas),Gv(Npar),H(Npar,Npar)
      Dimension Vb(Npar,Mxhv),Sb(Npar,Mxhv),Am(Mxhv,Mxhv),Wm(Mxhv)
      Dimension Cs(Nbas*Nbas),Gs(Npar),Vec(Npar),Res(Npar)
      Dimension F(Nbas,Nbas),U(Nbas,Nbas),Scr(Nbas,Nbas),G(Nbas,Nbas)
      Dimension Wk(Lwork),Rn(3),Irp(Nbas)


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



      Parameter (Delta=1.0D-3, Rtol=5.0D-4, Xnegst=-1.0D-3,
     +           Xlamth=-1.0D-4, Nguess=4, Mxroot=3)
      Parameter (Mxtot=200, Mkeep=Mxroot+2)
      Data Zero,Half,One,Two /0.0D0,0.5D0,1.0D0,2.0D0/

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
         Call Pccd_rfo_disp(Cs,Vb(1,Nhv),Delta,F,U,Scr,G,Irp,Nbas,
     +                        Npar)
         Call Putrec(20,"JOBARC","RFO_CSTB",Nsq*Iintfp,Cs)
         Call Putrec(20,"JOBARC","RFO_GSTB",Npar*Iintfp,Gs)
         Call Putrec(20,"JOBARC","RFO_ESTB",Iintfp,Es)
         Call Putrec(20,"JOBARC","RFO_VB  ",Npar*Mxhv*Iintfp,Vb)
         Call Putrec(20,"JOBARC","RFO_NHV ",1,Nhv)
         Call Putrec(20,"JOBARC","RFO_NHVT",1,0)
         Call Putrec(20,"JOBARC","RFO_HVPH",1,1)
         Call Putrec(20,"JOBARC","RFO_PHAS",1,1)
         Return
      Endif

      Call Getrec(20,"JOBARC","RFO_CSTB",Nsq*Iintfp,Cs)
      Call Getrec(20,"JOBARC","RFO_GSTB",Npar*Iintfp,Gs)
      Call Getrec(20,"JOBARC","RFO_ESTB",Iintfp,Es)
      Call Getrec(20,"JOBARC","RFO_VB  ",Npar*Mxhv*Iintfp,Vb)
      Call Getrec(20,"JOBARC","RFO_NHV ",1,Nhv)
      Call Getrec(20,"JOBARC","RFO_HVPH",1,Ihvph)
      If (Ihvph .Eq. 1) Then
C Gradient at +Delta: keep it and go to the -Delta point.
         Call Putrec(20,"JOBARC","RFO_GPLS",Npar*Iintfp,Gv)
         Call Pccd_rfo_disp(Cs,Vb(1,Nhv),-Delta,F,U,Scr,G,Irp,Nbas,
     +                      Npar)
         Call Putrec(20,"JOBARC","RFO_HVPH",1,2)
         Return
      Endif
      Call Putrec(20,"JOBARC","RFO_HVPH",1,1)
      Call Getrec(20,"JOBARC","RFO_GPLS",Npar*Iintfp,Res)
      If (Nhv .Gt. 1) Call Getrec(20,"JOBARC","RFO_SB  ",
     +                            Npar*(Nhv-1)*Iintfp,Sb)
      Do K = 1, Npar
         Sb(K,Nhv) = (Res(K) - Gv(K))/(Two*Delta)
      Enddo
      Call Putrec(20,"JOBARC","RFO_SB  ",Npar*Nhv*Iintfp,Sb)
      Call Getrec(20,"JOBARC","RFO_NHVT",1,Nhvt)
      Nhvt = Nhvt + 1
      Call Putrec(20,"JOBARC","RFO_NHVT",1,Nhvt)

      Ng = Min(Nguess-1,Npar)
      If (Nhv .Le. Ng) Then
C Further start vectors: the lowest eigenvectors of the model Hessian
C (RFO_GVEC, from Pccd_rfo_do), the Nhv-th one now.
         Call Getrec(20,"JOBARC","RFO_GVEC",Npar*Ng*Iintfp,Wk)
         Call Dcopy(Npar,Wk((Nhv-1)*Npar+1),1,Vec,1)
         Call Pccd_rfo_addvec(Vec,Vb,Npar,Mxhv,Nhv,Ok)
         If (Ok) Then
            Call Pccd_rfo_disp(Cs,Vb(1,Nhv),Delta,F,U,Scr,G,Irp,Nbas,
     +                        Npar)
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
      Write(6,"(a,i4,a,1p,3E12.4,a,3E9.2)") " @-RFO Davidson",Nhvt,
     +      " Hv: lowest roots",(Wm(Ir),Ir=1,Nroot),"  |r|",
     +      (Rn(Ir),Ir=1,Nroot)

      If (Alam .Lt. Xnegst .And. Rn(1) .Lt. 0.3D0*Dabs(Alam)) Then
         Iact = 2
      Elseif (Iunc .Eq. 0) Then
         Iact = 1
         If (Alam .Lt. Xlamth) Iact = 2
      Elseif (Nhvt .Ge. Mxtot) Then
         Iact = 3
         If (Alam .Lt. Xlamth) Iact = 2
      Else
         If (Nhv .Ge. Mxhv) Then
C Thick restart: keep the lowest Mkeep Ritz vectors and their H*b.
            M = Min(Mkeep,Nhv)
            Call Dgemm("N","N",Npar,M,Nhv,One,Vb,Npar,Am,Mxhv,Zero,
     +                 Wk,Npar)
            Call Dcopy(Npar*M,Wk,1,Vb,1)
            Call Dgemm("N","N",Npar,M,Nhv,One,Sb,Npar,Am,Mxhv,Zero,
     +                 Wk,Npar)
            Call Dcopy(Npar*M,Wk,1,Sb,1)
            Nhv = M
            Do J = 1, M
               Do I = 1, M
                  Am(I,J) = Zero
               Enddo
               Am(J,J) = One
            Enddo
            Call Putrec(20,"JOBARC","RFO_SB  ",Npar*Nhv*Iintfp,Sb)
            Write(6,"(a,i4,a)") " @-RFO Davidson restart after",Nhvt,
     +            " Hv products."
         Endif
C Expand with the preconditioned residual of the lowest unconverged root.
         Call Dgemv("N",Npar,Nhv,One,Sb,Npar,Am(1,Iunc),1,Zero,Res,1)
         Call Dgemv("N",Npar,Nhv,-Wm(Iunc),Vb,Npar,Am(1,Iunc),1,One,
     +              Res,1)
         Call Getrec(20,"JOBARC","RFO_PDIA",Npar*Iintfp,Wk(Npar+1))
         Do K = 1, Npar
            Den = Wk(Npar+K) - Wm(Iunc)
            If (Dabs(Den) .Lt. 1.0D-2) Den = Dsign(1.0D-2,Den)
            Wk(K) = -Res(K)/Den
         Enddo
         Call Pccd_rfo_addvec(Wk,Vb,Npar,Mxhv,Nhv,Ok)
         If (.Not. Ok) Call Pccd_rfo_addvec(Res,Vb,Npar,Mxhv,Nhv,Ok)
         If (Ok) Then
            Call Pccd_rfo_disp(Cs,Vb(1,Nhv),Delta,F,U,Scr,G,Irp,Nbas,
     +                        Npar)
            Call Putrec(20,"JOBARC","RFO_VB  ",Npar*Mxhv*Iintfp,Vb)
            Call Putrec(20,"JOBARC","RFO_NHV ",1,Nhv)
            Return
         Endif
         Iact = 3
         If (Alam .Lt. Xlamth) Iact = 2
      Endif

      If (Iact .Eq. 1) Write(6,"(a,1p,E12.4,a)")
     +   " @-RFO lowest Hessian eigenvalue",Alam,": minimum."
      If (Iact .Eq. 3) Write(6,"(a,i4,a,1p,E12.4,a)")
     +   " @-RFO stability check not converged after",Nhvt,
     +   " Hv products (lowest Ritz value",Alam,"): inconclusive."
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
