#!/usr/bin/env python3
"""
Standalone pilot probe for the MO-5 Pade/QP non-convergence
(H2O/STO-3G, post fold-fix, sym=off, eps_HF(5)=-0.3912734).

Reuses the driver's OWN 6 (z_i, f_i) sample points verbatim (copied by
hand from sigma_c_output.dat) and re-implements the Thiele
reciprocal-difference algorithm + forward recurrence evaluator EXACTLY
as in scrnc_sigma_driver.F's Scrnc_sigma_pade_coeff/Scrnc_sigma_pade_eval,
independently in Python (different language/arithmetic path -- same
spirit as [[feedback_pilot_code_debugging]]/[[feedback_consistency_check_debugging]]).

Purpose: reason about the MATH of the continuation/QP equation --
where are F(z)=z-eps_HF-Re[Sigma_c(z)]'s roots, where are the fitted
Pade function's poles -- not tune iteration damping to force a
"looks converged" number.
"""
import cmath

eps_HF = -0.3912734

# (omega, Re, Im) exactly as printed by the production code for p=q=5
data = [
    (0.0, -3.2821379466e-02, 0.0000000000e+00),
    (0.3, -3.3246090356e-02, 6.3682154327e-03),
    (0.6, -2.9838275414e-02, 1.8357143882e-02),
    (0.9, -1.7842779951e-02, 2.0131173731e-02),
    (1.2, -1.2581053023e-02, 1.4708933736e-02),
    (1.5, -1.3051371593e-02, 9.1639683764e-03),
]

N = len(data)
Zpt = [complex(0.0, om) for om, re, im in data]
Fpt = [complex(re, im) for om, re, im in data]


def pade_coeff(Zpt, Fpt):
    N = len(Zpt)
    G = [[0j] * N for _ in range(N)]
    for i in range(N):
        G[0][i] = Fpt[i]
    for k in range(1, N):
        for i in range(k, N):
            G[k][i] = (G[k - 1][k - 1] - G[k - 1][i]) / \
                      ((Zpt[i] - Zpt[k - 1]) * G[k - 1][i])
    return [G[k][k] for k in range(N)]


def pade_eval(Zpt, Acoef, Z):
    N = len(Zpt)
    A = [0j] * (N + 1)
    B = [0j] * (N + 1)
    A[0] = 0j
    A[1] = Acoef[0]
    B[0] = 1 + 0j
    B[1] = 1 + 0j
    for n in range(2, N + 1):
        A[n] = A[n - 1] + (Z - Zpt[n - 2]) * Acoef[n - 1] * A[n - 2]
        B[n] = B[n - 1] + (Z - Zpt[n - 2]) * Acoef[n - 1] * B[n - 2]
    return A[N] / B[N], B[N]


Acoef = pade_coeff(Zpt, Fpt)
print("Thiele coefficients a_k:")
for k, a in enumerate(Acoef, 1):
    print(f"  a{k} = {a}")

print()
print("=== Fixed-point iteration trace (production algorithm, exact) ===")
eps_qp = eps_HF
for it in range(1, 51):
    Fval, Bval = pade_eval(Zpt, Acoef, complex(eps_qp, 0.0))
    eps_new = eps_HF + Fval.real
    print(f"iter {it:2d}: eps_qp={eps_qp: .6f} -> "
          f"Sigma_c(z)={Fval.real: .6f}{Fval.imag:+.6f}i  "
          f"eps_new={eps_new: .6f}  |B(N)|={abs(Bval):.3e}")
    if abs(eps_new - eps_qp) < 1e-8:
        print("CONVERGED")
        eps_qp = eps_new
        break
    eps_qp = eps_new

print()
print("=== F(z) = z - eps_HF - Re[Sigma_c(z)] scan over real z ===")
print(f"{'z':>8} {'ReSigma':>12} {'ImSigma':>12} {'F(z)':>12} {'|B(N)|':>10}")
prev_F = None
prev_z = None
for i in range(-200, 201):
    z = i * 0.01
    Fval, Bval = pade_eval(Zpt, Acoef, complex(z, 0.0))
    Fz = z - eps_HF - Fval.real
    if prev_F is not None and (prev_F < 0) != (Fz < 0):
        print(f"  >>> sign change in F(z) between z={prev_z:.3f} and z={z:.3f}")
    if abs(Bval) < 0.05:
        print(f"  >>> near-pole: z={z:.3f}  |B(N)|={abs(Bval):.4e}")
    prev_F, prev_z = Fz, z
    if i % 20 == 0:
        print(f"{z:8.3f} {Fval.real:12.6f} {Fval.imag:12.6f} {Fz:12.6f} {abs(Bval):10.3e}")
