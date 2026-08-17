










      Subroutine Scrnc_freq_grid(N,Omega_c,Xnode,Wnode)
C
C Builds an N-point imaginary-frequency quadrature grid on [0,Infinity)
C for the contour-deformation self-energy integral (route decided
C 2026-08-17, see SCRNC_DEBUG_LOG.md "self-energy architecture" entry;
C grid choice -- Gauss-Legendre with a rational map, Ajith 2026-08-17 --
C decided over a minimax grid since minimax needs external tabulated
C points/weights this codebase doesn't have):
C
C   Integral(0,Infinity) f(omega') d(omega')
C     = Sum(k=1,N) Wnode(k) * f(Xnode(k))
C
C Standard rational map from Gauss-Legendre's native [-1,1] to
C [0,Infinity): x -> omega' = Omega_c*(1+x)/(1-x), with Jacobian
C d(omega')/dx = Omega_c*2/(1-x)**2, so Wnode(k) = Wleg(k) *
C Omega_c*2/(1-Xleg(k))**2. Omega_c sets the frequency scale (map's
C "half-width" point, omega'=Omega_c at x=0) -- pick it comparable to
C the orbital-energy differences appearing in the Sigma_c denominator
C (1/(i*omega+i*omega'-eps_m)) for a well-conditioned grid; not fixed
C by this routine, caller's responsibility.
C
C Uses a local copy of Scrnc_gaussl (Numerical-Recipes-derived
C Gauss-Legendre root/weight finder, originally aces2/vdint/gaussl.f,
C JG 4/93) rather than linking vdint -- same multiply-duplicated-
C per-module convention already established throughout this codebase
C (see GNUmakefile.src's SYMEXP/DODIIS comments).
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Xnode(N),Wnode(N)
      Dimension Xleg(N),Wleg(N)
C
      Call Scrnc_gaussl(-1.0D0,1.0D0,N,Xleg,Wleg)
C
      Do K = 1, N
         Xnode(K) = Omega_c * (1.0D0+Xleg(K)) / (1.0D0-Xleg(K))
         Wnode(K) = Wleg(K) * Omega_c * 2.0D0 /
     +              ((1.0D0-Xleg(K)) * (1.0D0-Xleg(K)))
      Enddo
C
      Return
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_gaussl(X1,X2,N,X,W)
C
C Local copy of aces2/vdint/gaussl.f (Numerical-Recipes-derived
C Gauss-Legendre root/weight finder, JG 4/93) -- see Scrnc_freq_grid's
C header comment above for why this is duplicated here rather than
C linked from vdint.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension X(N),W(N)
      Parameter (Eps=3.D-14)
      M = (N+1)/2
      Xm = 0.5D0*(X2+X1)
      Xl = 0.5D0*(X2-X1)
      Do I = 1, M
         Z = Cos(3.141592654D0*(I-0.25D0)/(N+0.5D0))
    1    Continue
            P1 = 1.D0
            P2 = 0.D0
            Do J = 1, N
               P3 = P2
               P2 = P1
               P1 = ((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
            Enddo
            Pp = N*(Z*P1-P2)/(Z*Z-1.D0)
            Z1 = Z
            Z = Z1-P1/Pp
         If (Abs(Z-Z1) .GT. Eps) Goto 1
         X(I)     = Xm-Xl*Z
         X(N+1-I) = Xm+Xl*Z
         W(I)     = 2.D0*Xl/((1.D0-Z*Z)*Pp*Pp)
         W(N+1-I) = W(I)
      Enddo
      Return
      End
