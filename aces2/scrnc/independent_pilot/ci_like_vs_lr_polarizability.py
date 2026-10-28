#!/usr/bin/env python3
"""
Independent PySCF test of Ajith's claim (2026-08-16, session after the
reciprocity_pilot.py entry): the alpha_xy/alpha_yx asymmetry seen in
ACES II's "CI-like" (IEOMPROP=0) CCSD polarizability is NOT a scrnc/ACES2
coding bug -- it is a generic consequence of CCSD being a NON-VARIATIONAL
method evaluated with a "CI-like" (expectation-value, not fully-relaxed
Lagrangian) property construction. Exact symmetry of a response tensor
only follows automatically (Clairaut/Schwarz symmetry of mixed partial
derivatives) when the quantity is a genuine second derivative of a
stationary Lagrangian/energy functional w.r.t. BOTH perturbation
directions -- true for variational methods (trivially), and true for CC
only when the full 2n+1-rule double-perturbation machinery is used (i.e.
"LR" proper -- which is what IEOMPROP=2/QUADRATIC is meant to be,
chronologically the last step after CI-like and LINEAR).

This script builds BOTH constructions, independently, in PySCF, on the
SAME distorted (C2v-broken) H2O/STO-3G/UHF/CCSD system, using only
PySCF's own trusted CCSD/Lambda/AO-integral machinery (no hand-rolled
response code):

  1. "CI-like" construction: freeze the ground-state Lambda amplitudes at
     their UNPERTURBED (zero-field) values; only relax T under the field.
     mu_i(F_j) = Tr[ D(T(F_j), Lambda(F=0)) * r_i ]
     alpha_ij = -d mu_i / dF_j  (central finite difference)
     This is NOT a genuine mixed second derivative of the energy (Lambda
     is not stationary w.r.t. F_j at the perturbed T), so nothing
     guarantees alpha_ij == alpha_ji.

  2. "LR" (fully relaxed) construction: re-solve Lambda self-consistently
     at EACH field value, same as reciprocity_pilot.py's proven method.
     mu_i(F_j) = Tr[ D(T(F_j), Lambda(F_j)) * r_i ]
     This IS (via Hellmann-Feynman on the fully relaxed density) a
     genuine dE/dF_i, so alpha_ij = -d^2E/dF_i dF_j is guaranteed
     symmetric by Schwarz's theorem, independent of any CC-specific
     physics.

Expectation: (1) shows a real, eps-independent (does NOT shrink toward
zero as the finite-difference step shrinks) alpha_xy/alpha_yx asymmetry.
(2) shows alpha_xy == alpha_yx to finite-difference precision, same
eps-scaling signature (shrinks ~eps^2) as reciprocity_pilot.py already
established for a different pair of operators.
"""
import numpy as np
from pyscf import gto, scf, cc

# Same distorted-water idea as the earlier ACES2-native EOM_NLO check:
# break C2v so alpha_xy/alpha_yx are not trivially zero by symmetry.
O = np.array([0., 0., 0.])
H1 = np.array([0.957, 0., 0.])
ANGLE = np.radians(100.0)
H2 = np.array([1.05 * np.cos(ANGLE), 1.05 * np.sin(ANGLE), 0.])

mol = gto.M(atom=[['O', O], ['H', H1], ['H', H2]], basis='sto-3g',
            symmetry=False, verbose=0, unit='Angstrom')

mf0 = scf.UHF(mol)
mf0.kernel()
mo_coeff = mf0.mo_coeff
mo_occ = mf0.mo_occ
hcore0 = mf0.get_hcore()
nao = mol.nao

# Dipole integrals in AO basis, electron contribution only (nuclear part
# is field-independent, drops out of any derivative w.r.t. field).
r_ao = mol.intor('int1e_r')  # shape (3, nao, nao); r_ao[0]=x, [1]=y, [2]=z
IX, IY = 0, 1


def base_cc(field_dir, eps):
    """UHF+UCCSD under a shifted core Hamiltonian H0 - eps*r_{field_dir},
    frozen orbitals (no SCF reoptimization). Returns the converged CC
    object (T amplitudes set, Lambda NOT yet solved)."""
    mf = scf.UHF(mol)
    mf.mo_coeff = mo_coeff
    mf.mo_occ = mo_occ
    shifted = hcore0 - eps * r_ao[field_dir]
    mf.get_hcore = lambda *args, **kwargs: shifted
    mycc = cc.UCCSD(mf)
    mycc.verbose = 0
    mycc.kernel()
    return mycc


def mu_ci_like(field_dir, eps, measure_dir, l1l2_unperturbed):
    """CI-like: T relaxed under the field, Lambda frozen at its
    unperturbed (F=0) value."""
    mycc = base_cc(field_dir, eps)
    mycc.l1, mycc.l2 = l1l2_unperturbed
    dm1a, dm1b = mycc.make_rdm1(t1=mycc.t1, t2=mycc.t2,
                                 l1=mycc.l1, l2=mycc.l2, ao_repr=True)
    dm_ao = dm1a + dm1b
    return np.einsum('pq,pq->', dm_ao, r_ao[measure_dir])


def mu_relaxed(field_dir, eps, measure_dir):
    """LR: both T and Lambda relaxed self-consistently under the field."""
    mycc = base_cc(field_dir, eps)
    mycc.solve_lambda()
    dm1a, dm1b = mycc.make_rdm1(ao_repr=True)
    dm_ao = dm1a + dm1b
    return np.einsum('pq,pq->', dm_ao, r_ao[measure_dir])


def unperturbed_lambda():
    mycc = base_cc(IX, 0.0)  # eps=0.0 => unperturbed
    mycc.solve_lambda()
    return mycc.l1, mycc.l2


def run(eps):
    l1l2_0 = unperturbed_lambda()

    # --- CI-like (frozen Lambda) ---
    muX_p = mu_ci_like(IY, +eps, IX, l1l2_0)
    muX_m = mu_ci_like(IY, -eps, IX, l1l2_0)
    alpha_xy_cilike = -(muX_p - muX_m) / (2 * eps)

    muY_p = mu_ci_like(IX, +eps, IY, l1l2_0)
    muY_m = mu_ci_like(IX, -eps, IY, l1l2_0)
    alpha_yx_cilike = -(muY_p - muY_m) / (2 * eps)

    # --- Fully relaxed (LR) ---
    muX_p_r = mu_relaxed(IY, +eps, IX)
    muX_m_r = mu_relaxed(IY, -eps, IX)
    alpha_xy_lr = -(muX_p_r - muX_m_r) / (2 * eps)

    muY_p_r = mu_relaxed(IX, +eps, IY)
    muY_m_r = mu_relaxed(IX, -eps, IY)
    alpha_yx_lr = -(muY_p_r - muY_m_r) / (2 * eps)

    print(f"eps = {eps:.1e}")
    print("  CI-like (Lambda frozen at F=0):")
    print(f"    alpha_xy = {alpha_xy_cilike:.8f}")
    print(f"    alpha_yx = {alpha_yx_cilike:.8f}")
    print(f"    diff     = {alpha_xy_cilike - alpha_yx_cilike:.3e}"
          f"   (rel {(alpha_xy_cilike - alpha_yx_cilike) / alpha_xy_cilike:.3e})")
    print("  Fully relaxed (LR, Lambda re-solved at each F):")
    print(f"    alpha_xy = {alpha_xy_lr:.8f}")
    print(f"    alpha_yx = {alpha_yx_lr:.8f}")
    print(f"    diff     = {alpha_xy_lr - alpha_yx_lr:.3e}"
          f"   (rel {(alpha_xy_lr - alpha_yx_lr) / alpha_xy_lr:.3e})")
    print()
    return (alpha_xy_cilike - alpha_yx_cilike,
            alpha_xy_lr - alpha_yx_lr)


if __name__ == '__main__':
    d1_cilike, d1_lr = run(2.0e-3)
    d2_cilike, d2_lr = run(1.0e-3)
    print("Ratio of asymmetry when eps is halved (expect ~1 for a real")
    print("effect, ~0.25 for O(eps^2) finite-difference truncation noise):")
    print(f"  CI-like: {d2_cilike / d1_cilike:.3f}")
    print(f"  LR:      {d2_lr / d1_lr:.3f}")
