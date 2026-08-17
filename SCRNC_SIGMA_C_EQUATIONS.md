# Contour-deformation Sigma_c(p,q;iomega) -- equations reference

Plain-text durable companion to the rendered version (published as a
Claude artifact, 2026-08-17, for actual review/discussion -- this file
exists so the content survives independent of that). Covers exactly
what is implemented in `aces2/scrnc/scrnc_sigma_prep.F` and
`aces2/scrnc/scrnc_sigma_driver.F` as of `SCRNC_DEBUG_LOG.md`
"continued 10" through "continued 13". **Status: NOT manuscript-ready**
-- working notes for discussion, see Section 11 for the open item that
needs a real decision before this is trusted further.

---

## 1. Where this starts: W(omega) from CC response

scrnc builds the dynamically screened Coulomb interaction directly
from CC/ring-CC linear response:

    W(i*omega') = v + v*chi(i*omega')*v                              (1)

held as the full unpacked AO-pair matrix `Jj_lrespn`, Nao_pairs =
Nao**2, `Iao_pair = (Mu-1)*Nao + Nu` (Mu outer, Nu inner -- pinned
down against `scrnc_prep_2ints.F`, not assumed). v is the bare Coulomb
kernel; chi is the already-validated CC/ring-CC response (Symmet2 fix,
shared Hbar^-1, imaginary-axis squaring, all from earlier sessions).
This document takes W(omega) as a given, already-verified input.

## 2. AO->MO transform of the interaction

Because `Jj_lrespn` is the full unpacked Kronecker square (not
symmetry-packed), the 4-AO-index -> 4-MO-index transform collapses to
two ordinary matrix multiplies against T = C(x)C:

    T((mu,nu),(p,q)) = C(mu,p)*C(nu,q)                                (2)
    W_mo = T^T * W * T                                                (3)

`Scrnc_wc_ao2mo.F` implements Eq. 3 as two block Xgemm passes (T is
never materialized). Verified against an independent brute-force
quadruple loop to machine precision on every test case (Section 10).

Fixed this session: `Jj_lrespn`'s argument was silently read as
Integer*8, not Double Precision (name starts with 'J', in Fortran's
implicit I-N integer range, not covered by `Implicit Double Precision
(A-H,O-Z)`). Root-caused via gdb, fixed by explicit declaration. See
SCRNC_DEBUG_LOG.md "continued 8".

## 3. Subtracting the bare Coulomb piece

    W_c(i*omega') = W(i*omega') - v                                   (4)

Done at the AO-pair level (valid since Eq. 3 is linear). Bare v built
from `Scrnc_get2ints`'s raw AO integrals, re-indexed into scrnc's own
convention -- safe because bare 2-electron integrals are already fully
permutation-symmetric by construction (confirmed by inspection, not
assumed).

Self-check, every run: bare-v matrix must be exactly symmetric under
(mu,nu)<->(lam,sig) pair swap. Measured max|v(a,b)-v(b,a)| = 0.0
exactly, every frequency node, all three systems tested (Section 10).

## 4. Reduction to the self-energy contraction

The self-energy sum only ever needs the diagonal-in-second/fourth-
index slice:

    Wc_reduced(p,q,m) = W_mo(p,m,q,m)                                 (5)

Nmo**3 instead of Nmo**4. This plus the alpha-spin SCFEVALA orbital
energies is everything written to `wcreduc.dat` for the external
driver.

## 5. The self-energy integral

Standard contour-deformation form, G0 diagonal in the MO basis
(G0_mn(i*nu) = delta_mn/(i*nu-eps_m)):

    Sigma_c(p,q; i*omega) = -(1/2*pi) * Integral(0,Infinity) d(omega')
         Sum_m [1/(i*omega+i*omega'-eps_m)] * W_c(p,m,q,m; i*omega')  (6)

Decided over the sum-over-poles/spectral alternative (Hbar
diagonalization) to stay closest to "W computed via CC linear-response
machinery" per the manuscript's own novelty framing, and to avoid
drifting toward the diagrammatic/upfolding GW-CC style this project
means to stay distinct from. See SCRNC_DEBUG_LOG.md "continued 5".

## 6. Imaginary-axis quadrature

Gauss-Legendre on [-1,1], rationally mapped to [0,Infinity):

    omega'_k = Omega_c * (1+x_k)/(1-x_k)                              (7a)
    w_k      = w_leg_k * 2*Omega_c/(1-x_k)**2                         (7b)

(x_k, w_leg_k) standard Gauss-Legendre nodes/weights (reused from
`aces2/vdint/gaussl.f`). Omega_c sets the map's scale (omega'=Omega_c
at x=0), should be comparable to the orbital-energy gaps in Eq. 6's
denominator for a well-conditioned grid -- fixed at 1.0 a.u. for every
run in Section 10; see Section 11 for why this matters more than it
looks.

Validated against two closed-form integrals (Integral 1/(1+omega'^2) =
pi/2, Integral exp(-omega') = 1) at N=8,16,24, Omega_c=0.5,2.0 BEFORE
ever touching real data -- both converge cleanly, down to 1.7D-13 and
4.1D-9 respectively at N=24.

## 7. Why N pipeline runs, not N x Nomega

Eq. 6 has two frequencies playing different roles: the integration
variable omega' (needs the full CC pipeline rerun, drives the actual
dynamical response) and the external evaluation point omega (only
enters algebraically, in the denominator). So the expensive part -- N
pipeline subprocess runs, one per quadrature node -- never depends on
how many external omega points are wanted: each node's
Wc_reduced(p,q,m) is read once and folded into EVERY requested
Sigma_c(p,q;i*omega_j) in the same pass, then discarded. Confirmed
directly: H2O/STO-3G with Nomega=6 still made exactly 4 pipeline
subprocess calls, same as Nomega=1.

Complex-plane picture (see the rendered artifact for the actual
diagram): the 4 sampled nodes for the H2O/DZ run were 0.075i, 0.493i,
2.030i, 13.40i. Continuing to a valence MO (eps_HF~0.6, close to the
fitted domain) is a modest, reliable continuation. Continuing to the O
1s core (eps_HF~-20.6, ~14x past the fitted [0,1.5] domain) is a real
extrapolation-domain violation -- this is exactly the Section 11
finding.

## 8. Pade continuation to the real axis

Thiele's reciprocal-difference continued fraction (Vidberg-Serene
1977) -- standard, not new. Given N sampled points (z_i, f_i =
Sigma_c(p,p;z_i)):

    g_1(z_i) = f_i
    g_k(z_i) = [g_(k-1)(z_(k-1)) - g_(k-1)(z_i)] /
               [(z_i - z_(k-1)) * g_(k-1)(z_i)]        for k=2..N       (8)
    a_k = g_k(z_k)

giving f(z) ~ a1/(1 + a2(z-z1)/(1 + a3(z-z2)/(1 + ...))), evaluated by
the standard stable forward two-term recurrence:

    A_n = A_(n-1) + (z-z_(n-1))*a_n*A_(n-2)
    B_n = B_(n-1) + (z-z_(n-1))*a_n*B_(n-2)                            (9)
    f(z) ~ A_N/B_N        [A_0=0, A_1=a_1, B_0=B_1=1]

Complex arithmetic throughout (Double Complex, explicit -- always
needs explicit declaration regardless of name, not subject to the
Section 2 implicit-typing trap).

Validated first, in isolation: sampled f(z)=1/(z-2) at 8 imaginary
points, continued back to 4 real points. Matched the exact closed form
to ~1D-15 at every point, BEFORE trusting it on real Sigma_c data.

## 9. The quasiparticle equation

Solved on-shell, real part only, plain fixed-point iteration from
eps_QP^(0) = eps_HF(p):

    eps_QP(p) = eps_HF(p) + Re[Sigma_c(p,p; eps_QP(p))]              (10)

No separate bare-exchange Sigma_x term. Deliberate: the usual G0W0 "Sigma_x
- v_xc" double-counting correction exists because a DFT reference's
v_xc is not exact exchange. This project's reference is UHF -- and for
CANONICAL HF orbitals, by definition,

    <p|h + J - Sigma_x_hat|q> = eps_HF(p) * delta_pq                 (11)

the Fock operator's static-exchange piece Sigma_x_hat IS the same
operator GW's own static-exchange piece would be. So Sigma_x - v_xc^HF
= 0 exactly, already folded into eps_HF(p) via Eq. 11. Standard result
for HF-based G0W0, not new physics, but it's the reason this codebase
builds no Sigma_x machinery at all.

## 10. Validation status

Three systems, same driver settings throughout (N=4, Omega_c=1.0,
Nomega=6):

| system        | Nao | bare-v asym. | AO2MO check      | QP converged | notes                          |
|---------------|-----|--------------|-------------------|--------------|---------------------------------|
| H2 / STO-3G   | 2   | 0.0          | 2.2D-16           | 2/2          | first end-to-end pass          |
| H2O / STO-3G  | 7   | 0.0          | 4D-15 - 1.1D-14   | 7/7          | lower symmetry, core+valence   |
| H2O / DZ      | 14  | 0.0          | 3.6D-14 - 6.5D-14 | 12/14        | 2 non-convergent -- Section 11 |

Bare-v symmetry and the AO2MO brute-force check stay clean at machine
precision across every system/basis tested -- Sections 2-4's machinery
shows no basis-size sensitivity. The QP solve is where a real
limitation shows up.

## 11. Open items for discussion

### Pade domain-of-validity (the headline finding) -- NOT silently patched

H2O/DZ's O 1s core (eps_HF=-20.559) failed Eq. 10 (hit the 50-iteration
cap): a ~1 Ha shift (vs 0.01-0.06 Ha for every other orbital) and a
NEGATIVE Im[Sigma_c] (every converged orbital came out positive).
Sigma_c(1,1;i*omega) on the fitted grid is small and nearly flat
(-9.40D-3 to -9.28D-3 across all 6 points) -- evaluating the Pade
approximant at z=-20.6, ~14x past the fitted omega' in [0,1.5] domain,
is a genuine extrapolation-domain violation, not a bug in Eqs. 8-9
themselves (independently exact to ~1D-15 in Section 8's own isolated
test).

A second orbital (p=7, eps_HF=0.311, WELL INSIDE the fitted domain)
also failed -- different mechanism. Its Sigma_c(i*omega) drops ~8x
between the first two sample points then flattens; a sharp early
feature in only 6 points can seat a spurious near-pole of the fitted
rational function close to the real axis, destabilizing the iteration
even in-domain.

**Deliberately not patched with iteration damping** -- would likely
make both cases LOOK converged without fixing the underlying Pade
issue, converging smoothly to a still-untrustworthy number. Candidate
real fixes, none yet chosen:
(a) scale Omega_c / the sampled omega' range to each orbital's own
    energy scale, rather than one fixed grid for every orbital;
(b) more quadrature/Pade points, especially near features like p=7's
    early kink;
(c) an explicit domain-of-validity check that flags/excludes any QP
    request landing far outside the fitted range.
This needs a real decision, not a solo one.

### Other flagged, unchanged since earlier sessions

- O(Nao**4) storage for the AO-pair W matrix -- fine at toy-system
  scale (largest so far: 196x196), not intended to scale to a
  realistic basis as-is.
- Alpha-spin only (SCFEVECA/SCFEVALA) -- inherited scope, never
  revisited for beta-spin or spin-summed treatment.
- Off-diagonal Sigma_c(p,q), p!=q is computed but unused -- Eq. 10
  only reads the diagonal; fine for every non-degenerate case tested,
  would need real degenerate-orbital treatment otherwise.
- Single geometry, minimal-basis territory only -- H2 and H2O at
  equilibrium, STO-3G/DZ. No stress test yet against anything a
  manuscript reviewer would call "a real molecule."

---

Source: `aces2/scrnc/scrnc_sigma_prep.F`, `scrnc_sigma_driver.F`.
Log: `SCRNC_DEBUG_LOG.md`, "continued 10" through "continued 13".
Commits: 548bed40, 30e1bf1c, ab6774ab, a4714a87, 742e0d66.
