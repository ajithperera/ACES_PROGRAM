#!/usr/bin/env python3
"""Pilot check on DZ p=1 (core, extrapolation-domain concern) and
p=7 (in-domain kink concern), post fold-fix + bracketed-rootfind.
Same independent Thiele reimplementation as pade_mo5_probe.py."""

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


def analyze(name, eps_HF, data, eps_qp_found):
    print(f"\n===== {name}: eps_HF={eps_HF:.6f}, driver found eps_QP={eps_qp_found:.6f} =====")
    Zpt = [complex(0.0, om) for om, re, im in data]
    Fpt = [complex(re, im) for om, re, im in data]
    Acoef = pade_coeff(Zpt, Fpt)

    def F(z):
        Fval, B = pade_eval(Zpt, Acoef, complex(z, 0.0))
        return z - eps_HF - Fval.real, Fval, B

    # Scan a wide window and report every sign change + min |B(N)| along the
    # path between eps_HF and the driver's reported root, to see whether a
    # pole sits between them (which would make the "found" root suspect).
    zlo = min(eps_HF, eps_qp_found) - 1.0
    zhi = max(eps_HF, eps_qp_found) + 1.0
    n = 4000
    dz = (zhi - zlo) / n
    prev_F = None
    prev_z = None
    minB = 1e30
    minB_z = None
    for i in range(n + 1):
        z = zlo + i * dz
        Fz, Fval, B = F(z)
        if abs(B) < minB:
            minB = abs(B)
            minB_z = z
        if prev_F is not None and (prev_F < 0) != (Fz < 0):
            print(f"  sign change: z in [{prev_z:.4f},{z:.4f}]  |B(N)| there ~ {abs(B):.4e}")
        prev_F, prev_z = Fz, z

    print(f"  min|B(N)| over [{zlo:.3f},{zhi:.3f}] = {minB:.4e} at z={minB_z:.4f}")
    Fz_at_root, Fval_at_root, B_at_root = F(eps_qp_found)
    print(f"  at driver's root z={eps_qp_found:.6f}: F(z)={Fz_at_root:.3e}, "
          f"Sigma_c={Fval_at_root.real:.6f}{Fval_at_root.imag:+.6f}i, |B(N)|={abs(B_at_root):.4e}")
    h = 1e-4
    _, Fp, _ = F(eps_qp_found + h)
    _, Fm, _ = F(eps_qp_found - h)
    gprime = (Fp.real - Fm.real) / (2 * h)
    print(f"  g'(root) = {gprime:.4f}")
    # sampled-domain edge
    om_max = max(om for om, re, im in data)
    print(f"  sampled imaginary-axis domain: |omega| in [0,{om_max}];  "
          f"|eps_HF|={abs(eps_HF):.3f}, |eps_QP|={abs(eps_qp_found):.3f} "
          f"({'INSIDE' if abs(eps_qp_found) <= om_max else 'OUTSIDE'} sampled |omega| range)")


data_p1 = [
    (0.0, -1.8792895223e-02, 0.0000000000e+00),
    (0.3, -1.8795276404e-02, 5.7711879512e-04),
    (0.6, -1.8756979404e-02, 1.2401949954e-03),
    (0.9, -1.8755427134e-02, 1.6151946284e-03),
    (1.2, -1.8724756439e-02, 1.7979896198e-03),
    (1.5, -1.8715312829e-02, 1.9682074384e-03),
]
analyze("DZ p=1 (core)", -20.55918, data_p1, -20.55232)

data_p7 = [
    (0.0, 2.8869579546e-03, 0.0000000000e+00),
    (0.3, 2.7181758865e-03, 3.3002449278e-03),
    (0.6, 3.3582308127e-03, 9.9280831507e-03),
    (0.9, 4.9250347128e-04, 1.0516704449e-02),
    (1.2, 1.8387892774e-04, 8.6628226105e-03),
    (1.5, 2.6611607862e-04, 6.8918390299e-03),
]
analyze("DZ p=7 (kink)", 0.310978, data_p7, 0.318185)
