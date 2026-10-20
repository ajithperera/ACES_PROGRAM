










      Subroutine Pccd_rfo_dav(H,Gv,Wr,Npar,Iaug,Nroot,X,Ev,Rtol,Wk,
     +                        Lwk,Iconv,Niter)
C
C Lowest Nroot eigenpairs of the RFO model by Davidson iterations, for
C large Npar where the dense diagonalizations (9 Npar**3 flops with
C vectors) dominate each macro iteration (C2H4 cc-pVTZ SYMMETRY=OFF,
C Npar = 6670: about 80% of an xpccd run in DSYEV). One matrix-vector
C product costs Npar**2.
C   Iaug = 1: the augmented Hessian [[H,g],[g^T,0]] (order Npar+1),
C             start vector e(Npar+1) (its first correction is the
C             diagonal Newton step);
C   Iaug = 0: H itself, start vectors the unit vectors of the lowest
C             diagonal elements plus one pseudo-random vector.
C   Iaug = 2: as 0, but the operator is the exact orbital Hessian
C             (Pccd_rsp_hv: H = frozen Hessian, Wr = its workspace);
C   Iaug = 3: as 1 with the exact orbital Hessian.
C X(N,Nroot) (N = Npar+1 for Iaug = 1, 3, else Npar): normalized
C eigenvectors, Ev: eigenvalues,
C converged when every residual norm is below Rtol. Wk(Lwk): workspace
C (the dense scratch Ah of Pccd_rfo_do). Iconv = 1 converged, 0 not
C converged within Maxit, -1 workspace too small; Niter = number of
C matrix-vector products.
C
      Implicit Double Precision(A-H,O-Z)
      Dimension H(Npar,Npar),Gv(*),Wr(*),X(*),Ev(Nroot),Wk(Lwk)
      Parameter (Mx=40)

      N    = Npar
      If (Iaug .Eq. 1 .Or. Iaug .Eq. 3) N = Npar + 1
      Lsw  = 64*Mx
      Iv   = 1
      Iav  = Iv  + N*Mx
      Idg  = Iav + N*Mx
      Ir   = Idg + N
      Iam  = Ir  + N
      Iam2 = Iam + Mx*Mx
      Iwm  = Iam2 + Mx*Mx
      Isw  = Iwm + Mx
      Iend = Isw + Lsw
      Iconv = -1
      Niter = 0
      If (Iend .Gt. Lwk .Or. 2*Nroot+1 .Gt. Mx .Or. Nroot .Gt. N) Return
      Call Pccd_rfo_dav1(H,Gv,Wr,Npar,Iaug,N,Nroot,X,Ev,Rtol,Wk(Iv),
     +                   Wk(Iav),Wk(Idg),Wk(Ir),Wk(Iam),Wk(Iam2),
     +                   Wk(Iwm),Wk(Isw),Lsw,Mx,Iconv,Niter)
      Return
      End

      Subroutine Pccd_rfo_dav1(H,Gv,Wr,Npar,Iaug,N,Nroot,X,Ev,Rtol,V,
     +                         Av,Dg,R,Am,Am2,Wm,Sw,Lsw,Mx,Iconv,Niter)
      Implicit Double Precision(A-H,O-Z)
      Dimension H(Npar,Npar),Gv(*),Wr(*),X(N,Nroot),Ev(Nroot),V(N,Mx),
     +          Av(N,Mx),Dg(N),R(N),Am(Mx,Mx),Am2(Mx,Mx),Wm(Mx),
     +          Sw(Lsw)
      Parameter (Maxit=600, Dmin=1.0D-4)
      Data Zero,One /0.0D0,1.0D0/

      Do K = 1, Npar
         Dg(K) = H(K,K)
      Enddo
      If (N .Gt. Npar) Dg(N) = Zero

C Start vectors.
      Call Dzero(V,N*Mx)
      If (N .Gt. Npar) Then
         V(N,1) = One
         Nv = 1
      Else
         Do Kr = 1, Nroot
            Kmin = 0
            Do K = 1, N
               If (V(K,Mx) .Eq. Zero) Then
                  If (Kmin .Eq. 0) Then
                     Kmin = K
                  Elseif (Dg(K) .Lt. Dg(Kmin)) Then
                     Kmin = K
                  Endif
               Endif
            Enddo
            V(Kmin,Kr) = One
C V(.,Mx) marks the used unit vectors until the subspace gets there.
            V(Kmin,Mx) = One
         Enddo
         Do K = 1, N
            V(K,Mx) = Zero
            V(K,Nroot+1) = Dsin(Dble(7*K+3))
         Enddo
         Nv = Nroot + 1
         Call Pccd_rfo_dorth(V,N,Nroot,V(1,Nv),Bn)
         Call Dscal(N,One/Bn,V(1,Nv),1)
      Endif
      Nb = 0

      Do It = 1, Maxit
C Products for the new subspace vectors and the new rows/columns of
C the projected matrix.
         Call Pccd_rfo_dmv(H,Gv,Wr,Npar,Iaug,N,Nv-Nb,V(1,Nb+1),
     +                     Av(1,Nb+1))
         Niter = Niter + Nv - Nb
         Do J = Nb+1, Nv
            Do I = 1, J
               Am(I,J) = Ddot(N,V(1,I),1,Av(1,J),1)
               Am(J,I) = Am(I,J)
            Enddo
            Do I = J+1, Nb
               Am(I,J) = Ddot(N,V(1,I),1,Av(1,J),1)
               Am(J,I) = Am(I,J)
            Enddo
         Enddo
         Nb = Nv
         Do J = 1, Nv
            Do I = 1, Nv
               Am2(I,J) = Am(I,J)
            Enddo
         Enddo
         Call Dsyev("V","U",Nv,Am2,Mx,Wm,Sw,Lsw,Info)
         If (Info .Ne. 0) Then
            Iconv = 0
            Return
         Endif

C Ritz vectors, residuals and the corrections of the unconverged roots
C (appended after the current Nv vectors).
         Nconv = 0
         Nnew  = Nv
         Do Kr = 1, Nroot
            Ev(Kr) = Wm(Kr)
            Call Dgemv("N",N,Nv,One,V,N,Am2(1,Kr),1,Zero,X(1,Kr),1)
            Call Dgemv("N",N,Nv,One,Av,N,Am2(1,Kr),1,Zero,R,1)
            Call Daxpy(N,-Wm(Kr),X(1,Kr),1,R,1)
            Rn = Dsqrt(Ddot(N,R,1,R,1))
            If (Rn .Lt. Rtol) Then
               Nconv = Nconv + 1
            Elseif (Nnew .Lt. Mx) Then
               Do K = 1, N
                  Den = Dg(K) - Wm(Kr)
                  If (Dabs(Den) .Lt. Dmin) Den = Dsign(Dmin,Den)
                  V(K,Nnew+1) = -R(K)/Den
               Enddo
C The correction scales with the residual: normalize it first, so that
C the linear-dependence test below is relative (an absolute 1e-8 threw
C away every correction once the residual fell below ~1e-8).
               Tn = Dsqrt(Ddot(N,V(1,Nnew+1),1,V(1,Nnew+1),1))
               Call Dscal(N,One/Tn,V(1,Nnew+1),1)
               Call Pccd_rfo_dorth(V,N,Nnew,V(1,Nnew+1),Bn)
               If (Bn .Gt. 1.0D-8) Then
                  Nnew = Nnew + 1
                  Call Dscal(N,One/Bn,V(1,Nnew),1)
               Endif
            Endif
         Enddo
         If (Nconv .Eq. Nroot) Then
            Iconv = 1
            Return
         Endif
         If (Nnew .Eq. Nv .Or. Nnew+Nroot .Gt. Mx) Then
C No room (or no new direction): collapse to the Ritz vectors (and the
C new corrections, if any fit) and recompute their products.
            Ncol = Nnew - Nv
            Do J = 1, Ncol
               Call Dcopy(N,V(1,Nv+J),1,Av(1,J),1)
            Enddo
            Do Kr = 1, Nroot
               Call Dcopy(N,X(1,Kr),1,V(1,Kr),1)
            Enddo
            Nv = Nroot
            Do J = 1, Ncol
               Call Pccd_rfo_dorth(V,N,Nv,Av(1,J),Bn)
               If (Bn .Gt. 1.0D-8 .And. Nv .Lt. Mx) Then
                  Nv = Nv + 1
                  Call Dcopy(N,Av(1,J),1,V(1,Nv),1)
                  Call Dscal(N,One/Bn,V(1,Nv),1)
               Endif
            Enddo
            Nb = 0
         Else
            Nv = Nnew
         Endif
      Enddo
      Iconv = 0
      Return
      End

      Subroutine Pccd_rfo_dmv(H,Gv,Wr,Npar,Iaug,N,Nc,Xv,Yv)
C Yv = A Xv for Nc vectors, A = H or [[H,g],[g^T,0]]; Iaug = 2, 3: H
C is the frozen Hessian and the amplitude relaxation is added (exact
C Hessian, Pccd_rsp_rlx). One DGEMM for the block: H is read once.
      Implicit Double Precision(A-H,O-Z)
      Dimension H(Npar,Npar),Gv(*),Wr(*),Xv(N,Nc),Yv(N,Nc)
      Data Zero,One /0.0D0,1.0D0/
      If (Nc .Le. 0) Return
      Call Dgemm("N","N",Npar,Nc,Npar,One,H,Npar,Xv,N,Zero,Yv,N)
      Do J = 1, Nc
         If (Iaug .Ge. 2) Call Pccd_rsp_rlx(Wr,Npar,Xv(1,J),Yv(1,J))
         If (N .Gt. Npar) Then
            Call Daxpy(Npar,Xv(N,J),Gv,1,Yv(1,J),1)
            Yv(N,J) = Ddot(Npar,Gv,1,Xv(1,J),1)
         Endif
      Enddo
      Return
      End

      Subroutine Pccd_rfo_dorth(V,N,Nv,B,Bn)
C Orthogonalize B against the orthonormal V(.,1:Nv) (twice); Bn = |B|.
      Implicit Double Precision(A-H,O-Z)
      Dimension V(N,Nv),B(N)
      Do Ipass = 1, 2
         Do J = 1, Nv
            Ov = Ddot(N,V(1,J),1,B,1)
            Call Daxpy(N,-Ov,V(1,J),1,B,1)
         Enddo
      Enddo
      Bn = Dsqrt(Ddot(N,B,1,B,1))
      Return
      End
