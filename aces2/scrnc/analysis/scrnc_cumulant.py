#!/usr/bin/env python3
"""
Post-processing of an xscrnc sigma run: Pade continuation of the
diagonal Sigma_c(p,p;iw) + GW+cumulant spectral function, reading ONLY
the on-disk sigma_c_output.dat / qp_energies.dat (no pipeline rerun).

The Pade part is a line-by-line port of scrnc_sigma_driver.F
(Scrnc_sigma_pade_coeff/_pade_eval, the inline node selection of
Scrnc_sigma_qp_solve, Scrnc_sigma_qp_bracket).  `--check` re-solves the
diagonal QP equation and compares against qp_energies.dat -- run it
before trusting anything else this script prints.

Cumulant (hole state; Aryasetiawan-Hedin-Karlsson 1996, Kas-Rehr-
Reining PRB 90, 085112 (2014)):
  beta(w') = (1/pi)|Im Sigma_c(eps_QP - w')|,  w' > 0
  C(t)     = Int_0^inf dw' beta(w')/w'^2 [exp(-i w' t) + i w' t - 1]
  A(w)     = (1/pi) Re Int_0^inf dt exp[-i(w-eps_QP)t - eta t] exp[C(t)]
Satellite weight is reported two ways:
  * window: 1 - Int_{|w-eps_QP|<win} A(w) dw   (eta-dependent, the
    2026-08-23 procedure)
  * split:  beta(w') = beta(0) + dbeta(w').  The constant beta(0)
    = |Im Sigma_c(eps_QP)|/pi part of C(t) integrates to a pure
    exponential decay (-pi beta(0) t/2 + log-phase), i.e. a LIFETIME
    broadening of the QP line, not satellite weight.  The genuinely
    frequency-dependent part gives an eta-independent incoherent weight
    1 - exp(-Int dbeta/w'^2), finite when dbeta ~ w' near 0.
"""
import argparse
import numpy as np

NSEL0, NSELSTEP, NSELMAX, TOL_STABLE = 6, 2, 20, 2.0e-2


def read_sigma(path):
    with open(path) as f:
        nmo, nom = (int(x) for x in f.readline().replace("=", " ").split()[1::2])
        omega = np.zeros(nom)
        sig = np.zeros((nmo, nmo, nom), complex)
        iw = -1
        for line in f:
            if line.startswith("Omega_ext"):
                iw += 1
                omega[iw] = float(line.split("=")[1].replace("D", "E"))
            elif line.strip() and line.split()[0].isdigit():
                p, q, re, im = line.replace("D", "E").split()
                sig[int(p) - 1, int(q) - 1, iw] = complex(float(re), float(im))
    return omega, sig


def read_qp(path):
    rows = []
    with open(path) as f:
        f.readline()
        for line in f:
            t = line.replace("D", "E").split()
            rows.append(dict(p=int(t[0]), eps_hf=float(t[1]), eps_qp=float(t[2]),
                             im=float(t[3]), conv=int(t[5]), stable=int(t[6]),
                             npade=int(t[7])))
    return rows


def select_nodes(omega, nuse, anchor):
    """Static point always kept + nuse-1 nearest to anchor, sorted ascending."""
    istat = int(np.argmin(omega))
    rest = [j for j in range(len(omega)) if j != istat]
    rest.sort(key=lambda j: abs(omega[j] - anchor))  # stable, like the insertion sort
    return sorted([istat] + rest[:nuse - 1], key=lambda j: omega[j])


class Pade:
    def __init__(self, z, f):
        n = len(z)
        g = np.zeros((n, n), complex)
        g[0] = f
        for k in range(1, n):
            g[k, k:] = (g[k - 1, k - 1] - g[k - 1, k:]) / ((z[k:] - z[k - 1]) * g[k - 1, k:])
        self.z, self.a = np.asarray(z), np.diag(g).copy()
        self.fscale = np.max(np.abs(f))

    def __call__(self, zz):
        zz = np.asarray(zz, complex)
        a_m2, a_m1 = np.zeros_like(zz), np.full_like(zz, self.a[0])
        b_m2, b_m1 = np.ones_like(zz), np.ones_like(zz)
        for n in range(1, len(self.a)):
            c = (zz - self.z[n - 1]) * self.a[n]
            a_m2, a_m1 = a_m1, a_m1 + c * a_m2
            b_m2, b_m1 = b_m1, b_m1 + c * b_m2
        return a_m1 / b_m1


def qp_bracket(pade, eps_hf):
    fmax = max(20.0 * pade.fscale, 1e-2)
    whalf = 2.0
    for _ in range(8):
        z = eps_hf - whalf + 2 * whalf / 200 * np.arange(201)
        fv = pade(z)
        F = z - eps_hf - fv.real
        best = None
        for i in range(200):
            if (F[i] < 0) != (F[i + 1] < 0) and abs(fv[i]) < fmax and abs(fv[i + 1]) < fmax:
                d = min(abs(z[i] - eps_hf), abs(z[i + 1] - eps_hf))
                if best is None or d < best[0]:
                    best = (d, z[i], z[i + 1])
        if best:
            lo, hi = best[1:]
            flo = lo - eps_hf - pade(lo).real
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                fm = mid - eps_hf - pade(mid).real
                if (fm < 0) == (flo < 0):
                    lo, flo = mid, fm
                else:
                    hi = mid
                if abs(hi - lo) < 1e-8:
                    break
            e = 0.5 * (lo + hi)
            return e, complex(pade(e))
        whalf *= 2
    return None, None


def diag_pade(omega, sig, p, nuse, anchor):
    idx = select_nodes(omega, nuse, anchor)
    return Pade(1j * omega[idx], sig[p, p, idx])


def solve_diag(omega, sig, p, eps_hf):
    """Port of Scrnc_sigma_qp_solve's escalation/stability loop for one MO."""
    nsel, prev = NSEL0, None
    while nsel <= min(NSELMAX, len(omega)):
        pade = diag_pade(omega, sig, p, nsel, abs(eps_hf))
        e, fv = qp_bracket(pade, eps_hf)
        if e is None:
            prev = None
        elif prev is not None and abs(e - prev) <= TOL_STABLE:
            return e, fv, nsel, pade
        else:
            prev = e
        nsel += NSELSTEP
    return e, fv, nsel - NSELSTEP, pade


def cumulant(pade, eps_qp, wcut, nw, tmax, nt, eta, split):
    wp = np.linspace(0, wcut, nw + 1)[1:]
    beta = np.abs(pade(eps_qp - wp).imag) / np.pi
    beta0 = abs(pade(eps_qp).imag) / np.pi
    if split:
        beta = beta - beta0
    t = np.linspace(0, tmax, nt)
    dw = wp[1] - wp[0]
    C = np.zeros(nt, complex)
    for k in range(0, nt, 2000):
        tt = t[k:k + 2000, None]
        x = wp[None, :] * tt
        C[k:k + 2000] = (beta / wp**2 * (np.exp(-1j * x) + 1j * x - 1)).sum(1) * dw
    a_inc = (beta / wp**2).sum() * dw
    return t, C, a_inc, beta0


def spectrum(t, C, w_rel, eta):
    dt = t[1] - t[0]
    g = np.exp(C - eta * t)
    wts = np.full(len(t), dt)
    wts[0] = wts[-1] = dt / 2
    A = np.array([(wts * np.exp(-1j * w * t) * g).sum().real for w in w_rel]) / np.pi
    return A


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rundir")
    ap.add_argument("--mo", type=int, nargs="*", help="1-based MO list (default: occupied)")
    ap.add_argument("--check", action="store_true", help="re-solve QP and compare to qp_energies.dat")
    ap.add_argument("--wcut", type=float, default=6.0)
    ap.add_argument("--nw", type=int, default=6000)
    ap.add_argument("--eta", type=float, nargs="*", default=[0.005, 0.0025])
    ap.add_argument("--win", type=float, default=3.0 / 27.211386, help="main-line half window (Ha)")
    ap.add_argument("--profile", action="store_true", help="print beta(w') profile")
    args = ap.parse_args()

    omega, sig = read_sigma(f"{args.rundir}/sigma_c_output.dat")
    qp = read_qp(f"{args.rundir}/qp_energies.dat")

    if args.check:
        print(" p   eps_QP(file)   eps_QP(py)    diff      Im(file)     Im(py)   npade file/py")
        for r in qp:
            p = r["p"] - 1
            e, fv, n, _ = solve_diag(omega, sig, p, r["eps_hf"])
            if e is None:
                print(f"{r['p']:2d}  no root (file conv={r['conv']})")
                continue
            print(f"{r['p']:2d} {r['eps_qp']:13.6f} {e:13.6f} {e - r['eps_qp']:9.1e} "
                  f"{r['im']:11.4e} {fv.imag:11.4e}   {r['npade']}/{n}")
        return

    mos = args.mo or [r["p"] for r in qp if r["eps_hf"] < 0]
    for m in mos:
        r = qp[m - 1]
        if not r["conv"]:
            print(f"MO{m}: QP not converged in this run, skipped")
            continue
        pade = diag_pade(omega, sig, m - 1, r["npade"], abs(r["eps_hf"]))
        e = r["eps_qp"]
        print(f"\nMO{m}: eps_HF={r['eps_hf']:.5f}  eps_QP={e:.5f} Ha  npade={r['npade']}  "
              f"Im Sigma(eps_QP)={pade(e).imag:.4e}")
        if args.profile:
            for x in (0.0, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0, 6.0):
                print(f"   w'={x:5.2f}  |ImSigma(eps_QP-w')|={abs(pade(e - x).imag):.4e}"
                      f"  ReSigma={pade(e - x).real:+.4e}")
        for split in (False, True):
            for eta in args.eta:
                tmax = 12.0 / eta
                nt = int(tmax / 0.2) + 1
                t, C, a_inc, beta0 = cumulant(pade, e, args.wcut, args.nw, tmax, nt, eta, split)
                w = np.linspace(-args.wcut - 0.5, 1.0, 6000)
                A = spectrum(t, C, w, eta)
                dw = w[1] - w[0]
                tot = A.sum() * dw
                main_ = A[np.abs(w) < args.win].sum() * dw
                lab = "split(dbeta) " if split else "full beta   "
                extra = f"  analytic Z-loss=1-exp(-a)={1 - np.exp(-a_inc):.4f}" if split else ""
                print(f"  {lab} eta={eta:.4f}: sum={tot:.4f}  outside+-win={(tot - main_) / tot:.4f}{extra}")
        print(f"  lifetime part: beta(0)={beta0:.4e} -> QP HWHM ~ {np.pi * beta0 / 2 * 27.2114:.3f} eV")


if __name__ == "__main__":
    main()
