#!/usr/bin/env python3
"""
Independent numerical cross-check for scrnc's CCSD response-function
symmetry bug (see SCRNC_DEBUG_LOG.md, 2026-08-16, "H2O/STO-3G/CCSD
asymmetric even with sym=off"). Not a re-implementation of scrnc's own
code, and not a claim that PySCF is more authoritative than ACES II --
ACES II's own SCF/CCSD/Lambda machinery is the established, decades-
validated reference this project already trusts. This pilot's only
value is that it is a genuinely SEPARATE codebase, so it cannot share
whatever bug is specific to scrnc's own response-construction Fortran
(Scrnc_form_mubar_s.F/_d_uhf.F/its Lineqy invocation).

Physics: for two one-electron operators A, B, the static CCSD linear
response function satisfies the reciprocity relation
    <<A;B>>(0) == <<B;A>>(0)
This is a direct consequence of the CCSD energy Lagrangian being a
smooth function of the perturbation strengths (Clairaut/Schwarz
symmetry of mixed partial derivatives, via the standard 2n+1-rule
Lagrangian derivation -- see Koch & Jorgensen, JCP 93, 3333 (1990)) --
not something specific to any one implementation's formulas.

Method: rather than hand-deriving and re-coding the CC linear-response
equations (Hbar solve + Lambda-dressed property-gradient vector) --
which would risk introducing an independent bug of its own, defeating
the point -- this pilot computes the SAME reciprocity via finite
difference of the CCSD+Lambda RELAXED one-particle density matrix,
using ONLY PySCF's own trusted CCSD/Lambda/make_rdm1 machinery
end-to-end:
    <<A;B>>(0) = d<A>/dlambda_B |_{lambda_B=0}
with a SHIFTED one-electron (core) Hamiltonian H_0 + lambda_B*B but
FROZEN HF orbitals (no SCF re-optimization) -- matching scrnc's own
"pure T/Lambda amplitude response, no orbital relaxation" formulation
(Lineqy solves only for the T-amplitude response; there is no CPHF/
Z-vector orbital-relaxation piece anywhere in scrnc).

Operators used: diagonal AO-pair "self" operators A_mu,mu matching the
SAME physical AOs behind scrnc's own worst-offender AO-pair indices
for H2O/STO-3G/UHF (ACES2 AO-pair index 9 = AO(2,2) = O 2s self-pair;
AO-pair index 41 = AO(6,6) = one H's 1s self-pair -- confirmed via
mol.ao_labels() that PySCF's own 0-based AO ordering for this same
STO-3G basis places these at indices 1 and 5 respectively, i.e. the
identical physical orbitals in the identical order).

Result (2026-08-16, eps=1e-3): d<A>/dlambda_B = -0.39660964,
d<B>/dlambda_A = -0.39660751, relative difference 5.4e-6. Halving eps
to 5e-4 shrinks the relative difference to 2.0e-6 (consistent with
O(eps^2) finite-difference truncation, not a real effect). scrnc's own
JJ(9,41) vs JJ(41,9) relative discrepancy is ~1.4e-3 -- about 1000x
larger than anything finite-difference truncation could produce here.
Conclusion: the reciprocity theorem genuinely holds for this system and
these operators (confirmed independently); scrnc's asymmetry is a real
implementation bug, not a missed physics effect (this also directly
refutes the earlier IEOMPROP=2 "missing quadratic correction"
hypothesis, which was tested and found to make things worse, not
better -- see the debug log).

Not yet done: this pilot only tests the SYMMETRY of the response
function, not scrnc's exact numerical VALUE (raw magnitudes differ,
sign and ~10% scale, most likely just a perturbation-strength/Fact=2
convention difference between "d<A>/dlambda_B for H=H0+lambda*B" and
whatever exact convention scrnc's Ddot(373,374)*Fact construction
uses -- not yet reconciled, and not needed for the symmetry conclusion
above).
"""
import numpy as np
from pyscf import gto, scf, cc

R = 0.957
ANGLE = np.radians(104.5)
O = np.array([0., 0., 0.])
H1 = np.array([R, 0., 0.])
H2 = np.array([R * np.cos(ANGLE), R * np.sin(ANGLE), 0.])
mol = gto.M(atom=[['O', O], ['H', H1], ['H', H2]], basis='sto-3g',
            symmetry=False, verbose=0, unit='Angstrom')

mf0 = scf.UHF(mol)
mf0.kernel()
mo_coeff = mf0.mo_coeff
mo_occ = mf0.mo_occ
hcore0 = mf0.get_hcore()
nao = mol.nao


def ao_pair_op(idx):
    m = np.zeros((nao, nao))
    m[idx, idx] = 1.0
    return m


A_ao = ao_pair_op(1)  # O 2s self-pair  (ACES2 AO-pair index 9)
B_ao = ao_pair_op(5)  # H1 1s self-pair (ACES2 AO-pair index 41)


def expectation_at(perturb_ao, eps, measure_ao):
    """Frozen orbitals, shifted core Hamiltonian, fresh CCSD+Lambda,
    return Tr[D_ao * measure_ao] and the alpha T1 amplitudes."""
    mf = scf.UHF(mol)
    mf.mo_coeff = mo_coeff
    mf.mo_occ = mo_occ
    shifted = hcore0 + eps * perturb_ao
    mf.get_hcore = lambda *args, **kwargs: shifted
    mycc = cc.UCCSD(mf)
    mycc.verbose = 0
    mycc.kernel()
    mycc.solve_lambda()
    dm1a, dm1b = mycc.make_rdm1(ao_repr=True)
    dm_ao = dm1a + dm1b
    return np.einsum('pq,pq->', dm_ao, measure_ao), mycc.t1


def run(eps):
    eA_p, t1_p = expectation_at(B_ao, +eps, A_ao)
    eA_m, _ = expectation_at(B_ao, -eps, A_ao)
    dA_dB = (eA_p - eA_m) / (2 * eps)

    eB_p, _ = expectation_at(A_ao, +eps, B_ao)
    eB_m, _ = expectation_at(A_ao, -eps, B_ao)
    dB_dA = (eB_p - eB_m) / (2 * eps)

    print(f"eps={eps:.1e}  T1(alpha) norm at +eps = "
          f"{np.linalg.norm(t1_p[0]):.6f}")
    print(f"  d<A>/dlambda_B = {dA_dB:.10f}")
    print(f"  d<B>/dlambda_A = {dB_dA:.10f}")
    print(f"  difference     = {dA_dB - dB_dA:.3e}")
    print(f"  relative diff  = {(dA_dB - dB_dA) / dA_dB:.3e}")
    print()


if __name__ == '__main__':
    run(1.0e-3)
    run(5.0e-4)
