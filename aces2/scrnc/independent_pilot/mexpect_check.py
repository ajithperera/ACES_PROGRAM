import numpy as np
from pyscf import gto, scf, cc

R = 0.957
Ang = np.radians(104.5)
O = np.array([0.,0.,0.])
H1 = np.array([R,0.,0.])
H2 = np.array([R*np.cos(Ang), R*np.sin(Ang), 0.])
mol = gto.M(atom=[['O',O],['H',H1],['H',H2]], basis='sto-3g', symmetry=False, verbose=0, unit='Angstrom')

mf = scf.UHF(mol)
mf.kernel()
mycc = cc.UCCSD(mf)
mycc.kernel()
l1, l2 = mycc.solve_lambda()

dm1a, dm1b = mycc.make_rdm1(ao_repr=True)
dm_ao = dm1a + dm1b

nao = mol.nao
def ao_pair_op(idx):
    M = np.zeros((nao, nao))
    M[idx, idx] = 1.0
    return M

for label, idx, scrnc_mexpect, scrnc_mbar00 in [
    ("AO-pair 9 (O 2s self)", 1, -0.021827, -0.008283),
    ("AO-pair 41 (H1 1s self)", 5, -0.016931, -0.005255),
]:
    M_ao = ao_pair_op(idx)
    mexpect = np.einsum('pq,pq->', dm_ao, M_ao)
    print(f"{label}: pyscf Tr[D*M] = {mexpect:.6f}   scrnc Mexpect = {scrnc_mexpect}   ratio = {mexpect/scrnc_mexpect:.4f}")

print()
print("--- testing correlation-only density hypothesis ---")
dm_hf_a, dm_hf_b = mf.make_rdm1()
dm_hf_ao = dm_hf_a + dm_hf_b
dm_corr_ao = dm_ao - dm_hf_ao

for label, idx, scrnc_mexpect, scrnc_mbar00 in [
    ("AO-pair 9 (O 2s self)", 1, -0.021827, -0.008283),
    ("AO-pair 41 (H1 1s self)", 5, -0.016931, -0.005255),
]:
    M_ao = ao_pair_op(idx)
    mexpect_corr = np.einsum('pq,pq->', dm_corr_ao, M_ao)
    print(f"{label}: pyscf Tr[(D_cc-D_hf)*M] = {mexpect_corr:.6f}   scrnc Mexpect = {scrnc_mexpect}   ratio = {mexpect_corr/scrnc_mexpect:.4f}")

# Result (2026-08-16): raw Tr[D_relaxed*M] is wrong order of magnitude
# vs scrnc's Mexpect. Tr[(D_cc-D_hf)*M] (correlation-only density) gets
# to the right order of magnitude (0.04 vs scrnc's 0.02-ish) but the
# ratio isn't a clean constant across AO pairs 9 vs 41 (-1.96 vs -2.43)
# -- a genuine remaining convention gap, most likely in how the AO-pair
# operator should be applied per spin block (this script sums dm1a+dm1b
# before contracting; Scrnc_form_mubar_00.F contracts per-ISPIN inside
# its own loop -- try replicating that structure exactly on resume).
# See SCRNC_DEBUG_LOG.md, "Mexpect convention reconciliation attempted".
