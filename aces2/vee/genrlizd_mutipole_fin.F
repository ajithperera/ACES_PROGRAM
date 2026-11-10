      Subroutine Genrlizd_multipole_fin(Q_strength,O_strength,
     +                                  A_strength,D_strength,
     +                                  M_strength,Root,D_iso,A_iso,
     +                                  O_iso,Q_iso,M_iso)
C
C Isotropically averaged oscillator-strength contributions from the
C products of genrlizd_multipole_init (omega = Root > 0, a.u.,
C alpha = 1/c):
C
C  f(mu)   =  (2/3) w Sum_a Dl Dr
C  f(m)    = -(1/6) a^2 w Sum_a Al Ar           (magnetic dipole,
C                                                orbital part)
C  f(Q)    =  (1/45) a^2 w^3 Sum_ab Theta_ab^2     (traceless Theta;
C             = (1/60) a^2 w^3 [3 Sum_ab Q_ab^2 - (Tr Q)^2], raw Q)
C  f(mu,O) = -(1/45) a^2 w^3 Sum_c (mu_c <r^2 r_c>)_sym
C  f(mu,M) = -(1/9)  a^2 w^2 Sum e_bad (mu_d M_ab)_antisym
C
C (from the expansion of exp(ik.r) e.p to second order in k, with
C [H,r]=-ip relating the e.p, (k.r)(e.p) and (k.r)^2(e.p) matrix
C elements to r, rr, rrr and {r,l}). The sum is origin independent
C for exact states. The sign of f(mu,M) is the one that makes the
C sum origin independent with the vee AO transition densities (the
C only sym x antisym contraction, so it depends on the D(pq) index
C order); checked for HeH+ aug-cc-pVQZ EOM-CCSD (= FCI) with the
C property origin moved by 1.5 bohr: total f(2) varies by < 0.5%
C while f(Q), f(mu,O), f(mu,M) individually change by up to 100%.
C
C 11/2026 (Ajith Perera): f(Q) used the raw-Q formula with the
C traceless Theta integrals (2.25 times too large); f(m) had +2/3 w (no a^2/4, wrong sign for
C the antisymmetric l), f(mu,O) had the wrong sign and f(mu,M) the
C 27-term sum (see init); single-precision fine-structure constant.
C
      Implicit Double Precision(A-H,O-Z)
      Double Precision M_strength,M_iso

      Dimension A_strength(3),Q_strength(9),O_strength(3)
      Dimension D_strength(3),M_strength(27)

      Data Dzero,Done,Two,Three,Six/0.0D0,1.0D0,2.0D0,3.0D0,6.0D0/
      Data Recp_fsc/137.035999084D0/

      Alpha2 = (Done/Recp_fsc)**2

      D_iso = Dzero
      A_iso = Dzero
      Do Ic = 1, 3
         D_iso = D_iso + D_strength(Ic)
         A_iso = A_iso + A_strength(Ic)
      Enddo
      D_iso = (Two/Three)*Root*D_iso
      A_iso = -(Done/Six)*Alpha2*Root*A_iso

C Sum_ab Theta_ab^2 = Sum_a Theta_aa^2 + 2 Sum_a<b Theta_ab^2
      Q_iso = (Q_strength(1)+Q_strength(2)+Q_strength(3))
     +      + Two*(Q_strength(4)+Q_strength(5)+Q_strength(6))
      Q_iso = (Alpha2/45.0D0)*Root**3*Q_iso

      O_iso = Dzero
      Do Ic = 1, 3
         O_iso = O_iso + O_strength(Ic)
      Enddo
      O_iso = -(Alpha2/45.0D0)*Root**3*O_iso

      M_iso = Dzero
      Do Ic = 1, 27
         M_iso = M_iso + M_strength(Ic)
      Enddo
      M_iso = -(Alpha2/9.0D0)*Root**2*M_iso

      Write(6,*)
      Write(6,100) "Electric dipole contribution                   :",
     +              D_iso
      Write(6,100) "Magnetic dipole contribution                   :",
     +              A_iso
      Write(6,100) "Electric quadrupole contribution               :",
     +              Q_iso
      Write(6,100) "Electric dipole/electric octupole contribution :",
     +             O_iso
      Write(6,100) "Electric dipole/magnetic quadrupole contrib.   :",
     +             M_iso
      Write(6,100) "Total (second order in k)                      :",
     +             D_iso+A_iso+Q_iso+O_iso+M_iso
 100  Format(1x,a,1x,es12.4e2)

      Return
      End
