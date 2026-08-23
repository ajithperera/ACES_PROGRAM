#!/usr/bin/env python3
"""
From-scratch analytic/symbolic RCCD & DRCCD derivation for the minimal
H2/STO-3G single-amplitude system, per Ajith's approved next step (see
SCRNC_DEBUG_LOG.md, "SESSION 2026-08-22 (autonomous fork #2, minimal-case
pilot...)" -- this is the flagged, not-yet-started "actual from-scratch
analytic or symbolic derivation" follow-up).

Goal: get a THIRD, independent number for the RHF-vs-UHF JJ-matrix
response gap (5.18x RCCD / 1.15x DRCCD, observed 2026-08-22), and pin
down the mechanism precisely rather than just re-observing the gap.

All numeric inputs below are REAL ACES II output, not invented --
pulled from run_scrnc.log in:
    /blue/ufhpc/perera/scrnc_test/h2_scrnc_rccd_rhf_sto3g_verify3
    /blue/ufhpc/perera/scrnc_test/h2_scrnc_rccd_uhf_sto3g_reverify
    /blue/ufhpc/perera/scrnc_test/h2_scrnc_drccd_rhf_sto3g_pilot0822
H2/STO-3G, R=0.742 Angstrom, sym=off, CALC=RCCD or DRCCD, DENSITY=RESPONSE.

This is a genuinely minimal case: 1 occupied + 1 virtual spatial MO for
RHF (H2 has 2 electrons total), so RCCD/DRCCD have exactly ONE non-zero
doubles amplitude (RHF: two spin-adapted "channel" versions of it,
singlet + triplet; UHF: the single ABAB spin-orbital amplitude, since
AAAA/BBBB are identically zero -- Pauli exclusion, only 1 alpha + 1
beta electron in the system).

RESULT SUMMARY (see script output and SCRNC_DEBUG_LOG.md entry for the
full writeup): the standard closed-shell spin-adaptation identity
    T(ABAB)_spin-orbital = [T_singlet - T_triplet] / 2
is verified EXACTLY (to ~1e-9) against ACES's own two independently-
computed numbers for T(ABAB) -- RHF's own internal spin-orbital
reconstruction check (Rcl_check_uhf) and UHF's own directly-converged
amplitude. This is a clean, from-first-principles, non-ACES-derived
confirmation that RHF's spin-adapted singlet/triplet channel amplitudes
and UHF's spin-orbital amplitude are exactly consistent representations
of the SAME physical ground state -- i.e. nothing is wrong with the
underlying rCCD amplitude solve itself, RHF vs UHF, at the amplitude
level.

A from-scratch reproduction of the numeric VALUES of T_singlet/
T_triplet/T_dRPA via the standard textbook closed-shell RPA/TDHF
(Casida-form) A,B matrices built from these same real MO integrals gets
the right sign/order of magnitude/qualitative pattern but only agrees
with ACES's own converged amplitudes to ~1-3% (not to machine
precision) -- see the "Riccati sanity check" section below. This is
good enough to confirm the basic single-amplitude Riccati-equation
PHYSICS is right, but not tight enough to trust a high-precision
downstream numeric prediction of the exact 5.18x/1.15x response-
magnitude ratio from these coefficients alone. The residual mismatch
most likely reflects ACES's actual internal RCCD/DRCCD contraction
(rcc_newt2.F et al., DPD-packed) using some normalization/intermediate
convention not identical to the plain textbook Casida formula -- tracing
that exactly was out of scope for this session (see log entry).
"""

import math

Ha_eV = 27.2113957

# ---------------------------------------------------------------------
# Real ACES MO data (RHF and UHF give BIT-IDENTICAL orbitals/integrals
# for this system -- confirmed: alpha and beta orbital eigenvalues in
# the UHF run are identical to each other and to the RHF run to all 10
# printed digits, and <S**2>=0.0000000 exactly. I.e. UHF does NOT break
# spin symmetry for H2/STO-3G at R=0.742A -- this is well inside the
# Coulson-Fischer stability radius. So the RHF/UHF response gap studied
# here is NOT a "different mean-field solution" effect -- it is purely
# about how the two pipelines construct/diagonalize the response Hbar
# from the SAME orbitals/integrals.)
# ---------------------------------------------------------------------
eps_i = -0.5777270364   # occupied MO energy, Hartree
eps_a = 0.6690811064    # virtual MO energy, Hartree
Delta = eps_a - eps_i

Jia = 0.66336391034632  # (ii|aa) Coulomb integral, chemist notation
Kia = 0.18132242419990  # (ia|ia) exchange integral, chemist notation

# ACES's own converged amplitudes/eigenvalues (ground truth to check
# the analytic model against):
T_S_aces      = -0.0976373126652   # RHF RCCD singlet-channel amplitude (list 63)
T_T_aces      =  0.163485120840    # RHF RCCD triplet-channel amplitude (list 61)
T_ABAB_uhf    = -0.130561216326    # UHF RCCD, directly converged ABAB amplitude
T_ABAB_rhfchk = -0.130561216375    # RHF's OWN internal Rcl_check_uhf reconstruction
T_dRPA_aces   = -0.0578172560853   # RHF DRCCD converged amplitude (single channel)

E_singlet_eV  = 25.50578519   # ringcl's own Hbar diagonalization, Imult=1 (singlet)
E_triplet_eV  = 15.48325067   # ringcl's own Hbar diagonalization, Imult=2 (triplet)
E_dRPA_eV     = 43.23234947   # ringcl's own Hbar diagonalization, DRCCD

JJresp_rhf_rccd  = 0.01161070026110285   # jjlrespn_resp_packed.dat diagonal, RHF RCCD
JJresp_uhf_rccd  = 0.06020595989104796   # ... UHF RCCD
JJresp_rhf_drccd = 0.02612050914120012   # ... RHF DRCCD
JJresp_uhf_drccd = 0.03008849944184275   # ... UHF DRCCD

print("=" * 70)
print("1. EXACT check: spin-orbital/spin-adapted amplitude identity")
print("=" * 70)
T_ABAB_predicted = (T_S_aces - T_T_aces) / 2.0
print(f"Predicted T(ABAB) = (T_S - T_T)/2 = {T_ABAB_predicted:.12f}")
print(f"ACES UHF (direct)                 = {T_ABAB_uhf:.12f}")
print(f"ACES RHF-internal reconstruction  = {T_ABAB_rhfchk:.12f}")
print(f"|predicted - UHF direct| = {abs(T_ABAB_predicted - T_ABAB_uhf):.3e}")
print("-> EXACT match (to ~4e-10), confirming RHF's singlet+triplet")
print("   channel amplitudes and UHF's spin-orbital amplitude are")
print("   exactly consistent representations of the same ground state.")

T_AAAA_predicted = (T_S_aces + T_T_aces) / 2.0
print(f"\n(For reference) predicted T(AAAA) = (T_S+T_T)/2 = {T_AAAA_predicted:.10f}")
print("ACES RHF-internal check T(AAAA) = 0.0329239037  (not a physically")
print("populated channel for this 2-electron system -- Rcl_check_uhf's")
print("own generic consistency diagnostic, not a real amplitude here)")

print()
print("=" * 70)
print("2. Riccati sanity check: textbook closed-shell RPA/TDHF (Casida)")
print("   A,B built from real ACES integrals, vs ACES's own converged")
print("   amplitudes and Hbar eigenvalues")
print("=" * 70)


def solve_riccati(A, B):
    """B t^2 + 2 A t + B = 0 (scalar rCCD/RPA Riccati equation, T=YX^-1
    reduced to a 1-amplitude system). Returns both roots."""
    disc = 4 * A * A - 4 * B * B
    sq = math.sqrt(disc)
    return (-2 * A + sq) / (2 * B), (-2 * A - sq) / (2 * B)


channels = {
    "singlet (RCCD, A=D+2K-J,B=K)": (Delta + 2 * Kia - Jia, Kia, T_S_aces, E_singlet_eV),
    "triplet (RCCD, A=D-J,  B=-K)": (Delta - Jia, -Kia, T_T_aces, E_triplet_eV),
    "dRPA    (DRCCD,A=D+2K, B=K)": (Delta + 2 * Kia, Kia, T_dRPA_aces, E_dRPA_eV),
}

for name, (A, B, T_ref, E_ref_eV) in channels.items():
    r1, r2 = solve_riccati(A, B)
    E_pred = (A + B * r1) * Ha_eV
    pct_T = 100 * (r1 - T_ref) / T_ref
    pct_E = 100 * (E_pred - E_ref_eV) / E_ref_eV
    print(f"\n{name}")
    print(f"  predicted T = {r1:.8f}   ACES T = {T_ref:.8f}   ({pct_T:+.2f}%)")
    print(f"  predicted E = {E_pred:.6f} eV  ACES E = {E_ref_eV:.6f} eV  ({pct_E:+.2f}%)")

print()
print("-> All three channels: right sign, right order of magnitude,")
print("   qualitatively correct pattern (dRPA < singlet < triplet in T")
print("   magnitude order observed also holds in the analytic model),")
print("   but only ~1-3% quantitative agreement -- NOT machine precision.")
print("   This is good enough to validate the basic single-amplitude")
print("   Riccati/A+BT structure, but not tight enough to trust a")
print("   precise downstream numeric prediction of the 5.18x/1.15x")
print("   response ratio built from these coefficients alone.")

print()
print("=" * 70)
print("3. Response-magnitude ratios actually observed (ACES, real runs)")
print("=" * 70)
r_rccd = JJresp_uhf_rccd / JJresp_rhf_rccd
r_drccd = JJresp_uhf_drccd / JJresp_rhf_drccd
print(f"RCCD  UHF/RHF response-only diagonal ratio  = {r_rccd:.4f}")
print(f"DRCCD UHF/RHF response-only diagonal ratio  = {r_drccd:.4f}")

print()
print("=" * 70)
print("4. Suggestive (UNCONFIRMED) mechanism probe: scrnc_form_mubar_s.F's")
print("   RHF-only '2*list37 - list34' dressing vs UHF's raw list37/T2AB")
print("=" * 70)
print("scrnc_form_mubar_s.F lines 286-306: for RHF (IUHF=0) the T2AB")
print("dressing fed into the transition-moment (mubar) construction is")
print("  SCR = 2*T2AB(list 37) - T2(list 34)")
print("while for UHF the raw T2AB (list 37) alone is used, no combination.")
print("scrnc's own Nmult=1 restriction (scrnc_driver.F line 251) means")
print("list 37, as read during the actual JJ-matrix run, IS the singlet")
print("channel amplitude T_S (Imult=1 is the only pass that runs).")
print()
print("run_scrnc.log shows list 34 is built by RCC_DWMBEJ specifically")
print("during the IMULT=2 (triplet) Hbar-construction pass inside xringcc")
print("(grep 'RCC_DWMBEJ RHF IMULT= 2 ... T2 list=  34') -- suggestive")
print("that list 34 may equal (or be closely derived from) T_T, the")
print("triplet-channel amplitude, but this was NOT independently verified")
print("this session (would need an instrumented rerun with an added")
print("print of list 34's actual numeric content, or full RNABIJ-resort")
print("tracing of rcc_dwmbej.F -- flagged as the concrete next step).")
print()
print("IF list34 == T_T (unconfirmed hypothesis), the RHF dressing factor")
print("would be 2*T_S - T_T, vs UHF's raw T_ABAB:")
combined_rhf = 2 * T_S_aces - T_T_aces
print(f"  2*T_S - T_T (hypothetical RHF dressing) = {combined_rhf:.6f}")
print(f"  T_ABAB (UHF dressing)                   = {T_ABAB_uhf:.6f}")
mubar_ratio_sq = ((1 + T_ABAB_uhf) / (1 + combined_rhf)) ** 2
print(f"  -> [(1+T_ABAB)/(1+2T_S-T_T)]^2 = {mubar_ratio_sq:.4f}")
print(f"     (observed RCCD ratio: {r_rccd:.4f})")
print()
print("Right DIRECTION (UHF bigger) and same order of magnitude as the")
print("real 5.18x, but undershoots substantially -- consistent with this")
print("being A real contributing mechanism, not necessarily the WHOLE")
print("story (the Nmult=1 singlet-only restriction itself, mechanism 1")
print("above, plausibly contributes independently on top of this).")
