#!/usr/bin/env python3
"""Restricted guess orbitals for OO-pCCD from a (broken-symmetry) UHF.

usage: mkguess.py <ca.dat> <cb.dat> <nocc> <outdir>

ca.dat/cb.dat: alpha/beta MOs (AO x MO, column-major, "nrow ncol" header)
from `xoorot` mode SAVEAB. Writes
  <outdir>/c0_alpha.dat : the UHF alpha MOs as they are
  <outdir>/c0_uno.dat   : UHF natural orbitals, sorted by occupation
and prints the UNO occupation numbers. The AO overlap is recovered from
the S-orthonormality of the alpha MOs (square C): S^-1 = Ca Ca^T.
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
    ca, cb = read_mat(sys.argv[1]), read_mat(sys.argv[2])
    nocc, out = int(sys.argv[3]), sys.argv[4]
    if ca.shape[0] != ca.shape[1]:
        sys.exit("mkguess: needs Nao == Nmo")
    s = np.linalg.inv(ca @ ca.T)
    for c, lab in ((ca, "alpha"), (cb, "beta")):
        err = np.abs(c.T @ s @ c - np.eye(c.shape[1])).max()
        print(f"{lab} orthonormality error {err:.2e}")
    p = ca[:, :nocc] @ ca[:, :nocc].T + cb[:, :nocc] @ cb[:, :nocc].T
    w, v = np.linalg.eigh(s)
    sh, smh = v @ np.diag(np.sqrt(w)) @ v.T, v @ np.diag(w ** -0.5) @ v.T
    n, u = np.linalg.eigh(sh @ p @ sh)
    order = np.argsort(-n)
    n, cno = n[order], smh @ u[:, order]
    print("UNO occupations (first 2*nocc):")
    print(" ".join(f"{x:.4f}" for x in n[:2 * nocc]))
    print(f"<S^2>-like overlap: sum_i<occ |<a_i|b_j>|^2 = "
          f"{np.sum((ca[:, :nocc].T @ s @ cb[:, :nocc]) ** 2):.4f} of {nocc}")
    write_mat(f"{out}/c0_alpha.dat", ca)
    write_mat(f"{out}/c0_uno.dat", cno)


if __name__ == "__main__":
    main()
