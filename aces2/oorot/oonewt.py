#!/usr/bin/env python3
"""Pilot driver: OO-pCCD with exact Newton steps from a finite-difference
Hessian of ACES II's analytic orbital gradient (OBRTGRDA).

All electronic structure is done by ACES II (xvscf/xvtran/xintprc/xpccd);
this script only displaces orbitals (C = C0 * expm(K), K antisymmetric,
applied through xoorot) and does the Newton linear algebra.

usage: oonewt.py check  <workdir>   # base point, U=I reproducibility,
                                    # energy-FD vs gradient convention
       oonewt.py newton <workdir>   # Newton iterations
"""
import os, sys, shutil, subprocess, glob
import numpy as np


def expm(K):
    """exp of a real antisymmetric matrix via eigh of the Hermitian iK."""
    w, V = np.linalg.eigh(1j * K)
    return (V @ np.diag(np.exp(-1j * w)) @ V.conj().T).real

ACES_STEPS = ["xvscf", "xvtran", "xintprc", "xpccd"]
DELTA_E = 1.0e-3      # energy finite-difference step (check mode)
DELTA_H = 1.0e-3      # gradient finite-difference step (Hessian)
MAXIT = 15
GTOL = 1.0e-6


def sh(cmd, cwd, log):
    with open(log, "a") as f:
        f.write(f"\n===== {cmd}\n"); f.flush()
        r = subprocess.run(cmd, shell=True, cwd=cwd, stdout=f, stderr=f)
    if r.returncode != 0:
        raise RuntimeError(f"{cmd} failed in {cwd} (rc={r.returncode})")


def read_mat(lines):
    n1, n2 = map(int, lines[0].split())
    vals = np.array(" ".join(lines[1:]).split(), dtype=float)
    return vals.reshape((n2, n1)).T          # column-major on disk


def write_mat(path, m):
    with open(path, "w") as f:
        f.write(f"{m.shape[0]} {m.shape[1]}\n")
        for v in m.T.ravel():
            f.write(f"{v:25.17E}\n")


def setup_base(wd):
    """HF SCF, save C0, then one xpccd at the HF orbitals. Snapshot."""
    base = os.path.join(wd, "base")
    log = os.path.join(base, "base.log")
    for f in os.listdir(base):
        if f not in ("ZMAT", "GENBAS"):
            os.remove(os.path.join(base, f))
    for s in ["xjoda", "xvmol", "xvmol2ja", "xvscf"]:
        sh(s, base, log)
    open(os.path.join(base, "oorot.in"), "w").write("SAVE\n")
    sh("xoorot", base, log)
    for s in ["xvtran", "xintprc", "xpccd"]:
        sh(s, base, log)
    open(os.path.join(base, "oorot.in"), "w").write("DUMP\n")
    sh("xoorot", base, log)
    snap = os.path.join(wd, "snap")
    if os.path.isdir(snap):
        shutil.rmtree(snap)
    shutil.copytree(base, snap)
    if os.environ.get("OONEWT_NHF0"):
        open(os.path.join(snap, "oorot.in"), "w").write("NHF0\n")
        sh("xoorot", snap, os.path.join(snap, "nhf0.log"))
    return read_point(base)


def read_point(d):
    lines = open(os.path.join(d, "grad.dat")).read().split("\n")
    lines = [l for l in lines if l.strip()]
    return float(lines[0]), read_mat(lines[1:])


def run_point(wd, tag, U):
    d = os.path.join(wd, "pts", tag)
    if os.path.isdir(d):
        shutil.rmtree(d)
    shutil.copytree(os.path.join(wd, "snap"), d)
    for junk in ["Energy", "base.log", "grad.dat"]:
        if os.path.exists(os.path.join(d, junk)):
            os.remove(os.path.join(d, junk))
    write_mat(os.path.join(d, "u.dat"), U)
    log = os.path.join(d, "pt.log")
    open(os.path.join(d, "oorot.in"), "w").write("APPLY\n")
    sh("xoorot", d, log)
    for s in ACES_STEPS:
        sh(s, d, log)
    open(os.path.join(d, "oorot.in"), "w").write("DUMP\n")
    sh("xoorot", d, log)
    open(os.path.join(d, "oorot.in"), "w").write("DUMPF\n")
    sh("xoorot", d, log)
    E, G = read_point(d)
    scf = [l for l in open(log, errors="replace") if "E(SCF)=" in l]
    # keep only small records; drop the big integral files
    for f in glob.glob(os.path.join(d, "*")):
        if os.path.basename(f) not in ("pt.log", "grad.dat", "u.dat",
                                       "fock.dat"):
            if os.path.isfile(f):
                os.remove(f)
    return E, G, (scf[-1].strip() if scf else "")


def kmat(n, p, q, x):
    K = np.zeros((n, n)); K[p, q] = x; K[q, p] = -x
    return K


def gvec(G, pairs):
    """Antisymmetric part of OBRTGRDA on the independent pairs (p<q)."""
    return np.array([G[p, q] - G[q, p] for p, q in pairs])


def check(wd):
    E0, G0 = setup_base(wd)
    n = G0.shape[0]
    print(f"base: E = {E0:.12f}  nmo = {n}")
    print(f"OBRTGRDA symmetric-part max = {np.abs(G0+G0.T).max()/2:.3e}, "
          f"antisym-part max = {np.abs(G0-G0.T).max()/2:.3e}")
    np.set_printoptions(linewidth=160, precision=6, suppress=False)
    print("OBRTGRDA (base):\n", G0)
    E1, G1, s1 = run_point(wd, "identity", np.eye(n))
    print(f"U=I rerun: E = {E1:.12f}  dE = {E1-E0:.2e}  "
          f"max|dG| = {np.abs(G1-G0).max():.2e}   {s1}")
    # pick largest-|antisym grad| pair in each block, plus one zero pair
    with open(os.path.join(wd, "nocc")) as f:
        nocc = int(f.read())
    A = G0 - G0.T
    blocks = {"oo": [(p, q) for p in range(nocc) for q in range(p+1, nocc)],
              "vv": [(p, q) for p in range(nocc, n) for q in range(p+1, n)],
              "ov": [(p, q) for p in range(nocc) for q in range(nocc, n)]}
    for b, prs in blocks.items():
        if not prs:
            continue
        p, q = max(prs, key=lambda t: abs(A[t]))
        Ep, _, sp = run_point(wd, f"{b}_{p}_{q}_p", expm(kmat(n, p, q, DELTA_E)))
        Em, _, sm = run_point(wd, f"{b}_{p}_{q}_m", expm(kmat(n, p, q, -DELTA_E)))
        dE = (Ep - Em) / (2 * DELTA_E)
        ratio = dE / A[p, q] if abs(A[p, q]) > 1e-12 else float("nan")
        print(f"{b} ({p},{q}): FD dE/dk = {dE: .8e}  "
              f"G[p,q]-G[q,p] = {A[p,q]: .8e}  ratio = {ratio: .6f}")
        print(f"      E(+)-E0 = {Ep-E0: .3e}, E(-)-E0 = {Em-E0: .3e}   {sp}")


def newton(wd, ratio):
    E0, G0 = setup_base(wd)
    n = G0.shape[0]
    pairs = [(p, q) for p in range(n) for q in range(p+1, n)]
    npar = len(pairs)
    Uacc = np.eye(n)
    E, G = E0, G0
    print(f"Newton: nmo={n}, {npar} rotation parameters, ratio={ratio}")
    for it in range(MAXIT):
        g = ratio * gvec(G, pairs)
        print(f"it {it:2d}  E = {E:.12f}  |g|max = {np.abs(g).max():.3e} "
              f" rms = {np.sqrt((g**2).mean()):.3e}", flush=True)
        if np.abs(g).max() < GTOL:
            break
        H = np.zeros((npar, npar))
        for k, (p, q) in enumerate(pairs):
            _, Gp, _ = run_point(wd, f"h{it}_{k}p",
                                 Uacc @ expm(kmat(n, p, q, DELTA_H)))
            _, Gm, _ = run_point(wd, f"h{it}_{k}m",
                                 Uacc @ expm(kmat(n, p, q, -DELTA_H)))
            H[:, k] = ratio * (gvec(Gp, pairs) - gvec(Gm, pairs)) / (2*DELTA_H)
        asym = np.abs(H - H.T).max()
        H = 0.5 * (H + H.T)
        w, V = np.linalg.eigh(H)
        np.save(os.path.join(wd, f"hess_{it}.npy"), H)
        print(f"      Hessian asym = {asym:.2e}; lowest eigenvalues: "
              f"{np.array2string(w[:6], precision=5)}", flush=True)
        step = -np.linalg.solve(H, g)
        smax = np.abs(step).max()
        if smax > 0.2:
            step *= 0.2 / smax
        K = np.zeros((n, n))
        for (p, q), x in zip(pairs, step):
            K[p, q] = x; K[q, p] = -x
        Uacc = Uacc @ expm(K)
        E, G, s = run_point(wd, f"it{it+1}", Uacc)
        print(f"      step max = {smax:.3e}   {s}", flush=True)
    np.save(os.path.join(wd, "U_final.npy"), Uacc)
    print("final E = %.12f" % E)


def fd_hessian(wd, it, Uacc, pairs, ratio, n):
    npar = len(pairs)
    H = np.zeros((npar, npar))
    for k, (p, q) in enumerate(pairs):
        _, Gp, _ = run_point(wd, f"h{it}_{k}p", Uacc @ expm(kmat(n, p, q, DELTA_H)))
        _, Gm, _ = run_point(wd, f"h{it}_{k}m", Uacc @ expm(kmat(n, p, q, -DELTA_H)))
        H[:, k] = ratio * (gvec(Gp, pairs) - gvec(Gm, pairs)) / (2 * DELTA_H)
    return H


def mode_report(w, V, pairs, nocc, nneg=3):
    for m in range(min(nneg, len(w))):
        if w[m] >= 0:
            break
        top = np.argsort(-np.abs(V[:, m]))[:4]
        desc = ", ".join(f"({pairs[t][0]+1},{pairs[t][1]+1}){V[t, m]:+.2f}"
                         for t in top)
        print(f"      neg mode {m}: {w[m]: .5f}  [{desc}]")


def rfo(wd, ratio, trust=0.1):
    """Trust-region augmented-Hessian (RFO) minimization, exact FD Hessian."""
    E0, G0 = setup_base(wd)
    n = G0.shape[0]
    nocc = int(open(os.path.join(wd, "nocc")).read())
    pairs = [(p, q) for p in range(n) for q in range(p+1, n)]
    Uacc, E, G = np.eye(n), E0, G0
    print(f"RFO: nmo={n}, {len(pairs)} parameters, start E = {E0:.12f}")
    for it in range(40):
        g = ratio * gvec(G, pairs)
        print(f"it {it:2d}  E = {E:.12f}  |g|max = {np.abs(g).max():.3e}"
              f"  trust = {trust:.3f}", flush=True)
        if np.abs(g).max() < GTOL:
            break
        H = fd_hessian(wd, it, Uacc, pairs, ratio, n)
        H = 0.5 * (H + H.T)
        w, V = np.linalg.eigh(H)
        np.save(os.path.join(wd, f"hess_{it}.npy"), H)
        print(f"      Hessian lowest eig: {np.array2string(w[:5], precision=5)}"
              f"  (#neg = {(w < 0).sum()})", flush=True)
        mode_report(w, V, pairs, nocc)
        while True:
            npar = len(g)
            AH = np.zeros((npar + 1, npar + 1))
            AH[:npar, :npar] = H; AH[:npar, npar] = g; AH[npar, :npar] = g
            a, B = np.linalg.eigh(AH)
            step = B[:npar, 0] / B[npar, 0]
            if np.linalg.norm(step) > trust:
                step *= trust / np.linalg.norm(step)
            pred = g @ step + 0.5 * step @ H @ step
            K = np.zeros((n, n))
            for (p, q), x in zip(pairs, step):
                K[p, q] = x; K[q, p] = -x
            Unew = Uacc @ expm(K)
            En, Gn, s = run_point(wd, f"it{it+1}_t{trust:.4f}", Unew)
            rho = (En - E) / pred if pred != 0 else 0.0
            print(f"      |step| = {np.linalg.norm(step):.3e}  dE = {En-E: .3e}"
                  f"  pred = {pred: .3e}  rho = {rho:.2f}   {s}", flush=True)
            if En < E or trust < 1e-4:
                break
            trust *= 0.5
        Uacc, E, G = Unew, En, Gn
        if rho > 0.75 and np.linalg.norm(step) > 0.8 * trust:
            trust = min(2 * trust, 0.5)
        elif rho < 0.25:
            trust *= 0.5
    np.save(os.path.join(wd, "U_final.npy"), Uacc)
    H = fd_hessian(wd, "final", Uacc, pairs, ratio, n)
    w, V = np.linalg.eigh(0.5 * (H + H.T))
    print(f"final E = {E:.12f}; final Hessian lowest eig: "
          f"{np.array2string(w[:5], precision=5)}  (#neg = {(w < 0).sum()})")
    mode_report(w, V, pairs, nocc)


def gradcheck(wd, ratio, scale=0.05, seed=7, block="all", tag="gc"):
    """Full FD-energy vs analytic-gradient check at a rotated point;
    the displacement is restricted to one block (oo, vv, ov) or all."""
    E0, G0 = setup_base(wd) if tag == "gc" or not os.path.isdir(
        os.path.join(wd, "snap")) else read_point(os.path.join(wd, "base"))
    n = G0.shape[0]
    nocc = int(open(os.path.join(wd, "nocc")).read())
    pairs = [(p, q) for p in range(n) for q in range(p+1, n)]
    def blk(p, q):
        return "oo" if q < nocc else ("vv" if p >= nocc else "ov")
    rng = np.random.default_rng(seed)
    K = np.zeros((n, n))
    for p, q in pairs:
        x = scale * rng.standard_normal()
        if block in ("all", blk(p, q)):
            K[p, q] = x; K[q, p] = -x
    U0 = expm(K)
    Ec, Gc, s = run_point(wd, f"{tag}_center", U0)
    g = ratio * gvec(Gc, pairs)
    print(f"center: E = {Ec:.12f}  {s}")
    print(" pair      analytic          FD(E)           diff")
    worst = 0.0
    for k, (p, q) in enumerate(pairs):
        Ep, _, _ = run_point(wd, f"{tag}_{k}p", U0 @ expm(kmat(n, p, q, DELTA_E)))
        Em, _, _ = run_point(wd, f"{tag}_{k}m", U0 @ expm(kmat(n, p, q, -DELTA_E)))
        fd = (Ep - Em) / (2 * DELTA_E)
        worst = max(worst, abs(fd - g[k]))
        print(f"({p},{q}) {blk(p,q)}  {g[k]: .8e}  {fd: .8e}  {fd-g[k]: .2e}", flush=True)
    print(f"max |FD - analytic| = {worst:.3e}")


if __name__ == "__main__":
    mode, wd = sys.argv[1], os.path.abspath(sys.argv[2])
    if mode == "check":
        check(wd)
    elif mode == "rfo":
        rfo(wd, float(sys.argv[3]))
    elif mode == "gradcheck":
        gradcheck(wd, float(sys.argv[3]))
    elif mode == "gradov":
        gradcheck(wd, float(sys.argv[3]), block="ov", tag="gc")
    elif mode == "gradblocks":
        r = float(sys.argv[3])
        for i, b in enumerate(["oo", "vv", "ov"]):
            print(f"\n######## displacement restricted to {b} block")
            gradcheck(wd, r, block=b, tag=("gc" if i == 0 else f"gc{b}"))
    else:
        newton(wd, float(sys.argv[3]) if len(sys.argv) > 3 else 1.0)
