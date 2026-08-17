










      Subroutine Scrnc_freq_grid_test
C
C Standing diagnostic (guarded by 1, permanent, not
C exploratory cruft -- same convention as Scrnc_bcktran_test.F),
C validating Scrnc_freq_grid against two closed-form integrals on
C [0,Infinity):
C
C   Integral 1/(1+omega**2) d(omega) = Pi/2   (slowly-decaying, tests
C                                              the map's large-omega tail)
C   Integral Exp(-omega)    d(omega) = 1      (fast-decaying, tests the
C                                              map's small-x/near-origin
C                                              region)
C
C Run at three grid sizes (N=8,16,24) and two Omega_c scales (0.5,2.0
C a.u.) to check both N-convergence and that the result isn't
C pathologically sensitive to the frequency-scale choice a caller picks.
C This is pure numerics (no molecular/Work-array data needed), unlike
C every other _test routine in this module -- kept as a plain
C Implicit-None-free Subroutine anyway, matching this module's overall
C convention, called from Scrnc_driver.F.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Xnode(24),Wnode(24)
      Dimension Ntest(3)
      Data Ntest /8,16,24/
C
      Pi = 3.14159265358979324D0
      Exact1 = Pi/2.0D0
      Exact2 = 1.0D0
C
      Write(6,*)
      Write(6,"(a)") " ---- Scrnc_freq_grid validation ----"
C
      Do Iscale = 1, 2
         If (Iscale .EQ. 1) Omega_c = 0.5D0
         If (Iscale .EQ. 2) Omega_c = 2.0D0
C
         Do Itest = 1, 3
            N = Ntest(Itest)
            Call Scrnc_freq_grid(N,Omega_c,Xnode,Wnode)
C
            Sum1 = 0.0D0
            Sum2 = 0.0D0
            Do K = 1, N
               Sum1 = Sum1 + Wnode(K)/(1.0D0+Xnode(K)*Xnode(K))
               Sum2 = Sum2 + Wnode(K)*Exp(-Xnode(K))
            Enddo
C
            Write(6,"(a,F5.2,a,I3,a,D12.5,a,D12.5)")
     +         " Omega_c=",Omega_c," N=",N,
     +         "  |sum1-Pi/2|=",Abs(Sum1-Exact1),
     +         "  |sum2-1|=",Abs(Sum2-Exact2)
         Enddo
      Enddo
C
      Write(6,"(a)") " ---- end Scrnc_freq_grid validation ----"
      Write(6,*)
C
      Return
      End
