#!/usr/bin/env python3
"""Broken-symmetry UHF guess for a homonuclear diatomic from RHF MOs.

usage: mkbsguess.py <c0.dat> <nocc> <npair> <nao_atom1> <outdir>

c0.dat: RHF MOs (AO x MO, column-major, "nrow ncol" header; xoorot SAVE).
W = span(npair highest occupied, npair lowest virtual MOs), i.e. the
bonding/antibonding orbitals of the bonds to break. The Loewdin projector
onto atom 1 (AOs 1..nao_atom1) is diagonalized in W; this is invariant
to rotations within degenerate (pi) sets. Its npair largest-eigenvalue
vectors (on atom 1) replace the alpha occupied part of W, the npair
smallest (on atom 2) the beta occupied part; the remaining vectors fill
the virtual slots. Writes ca.dat, cb.dat.
"""
import sys
import numpy as np


def read_mat(path):
    lines = [l for l in open(path).read().split("\n") if l.strip()]
    n1, n2 = map(int, lines[0].split())
    vals = np.array(" ".join(lines[1:]).split(), dtype=float)
    return vals.reshape((n2, n1)).T


def write_mat(path, m):
    with open(path, "w") as f:
        f.write(f"{m.shape[0]} {m.shape[1]}\n")
        for v in m.T.ravel():
            f.write(f"{v:25.17E}\n")


def main():
    c = read_mat(sys.argv[1])
    nocc, npair, na1, out = (int(sys.argv[2]), int(sys.argv[3]),
                             int(sys.argv[4]), sys.argv[5])
    s = np.linalg.inv(c @ c.T)
    w, v = np.linalg.eigh(s)
    sh = v @ np.diag(np.sqrt(w)) @ v.T
    slots = list(range(nocc - npair, nocc + npair))
    cw = c[:, slots]
    xw = sh @ cw                                  # Loewdin-basis coefficients
    m = xw[:na1, :].T @ xw[:na1, :]
    lam, u = np.linalg.eigh(m)                    # ascending
    print("atom-1 Loewdin populations in W:",
          " ".join(f"{x:.3f}" for x in lam))
    on1, on2 = cw @ u[:, ::-1][:, :npair], cw @ u[:, :npair]
    occ_slots, vir_slots = slots[:npair], slots[npair:]
    ca, cb = c.copy(), c.copy()
    ca[:, occ_slots], ca[:, vir_slots] = on1, on2
    cb[:, occ_slots], cb[:, vir_slots] = on2, on1
    write_mat(f"{out}/ca.dat", ca)
    write_mat(f"{out}/cb.dat", cb)


if __name__ == "__main__":
    main()
