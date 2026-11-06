#!/usr/bin/env python3
"""EOM-(D)RCCD excited-state regression suite -- reproduces Table I of
Rishi, Perera, and Bartlett, "A route to improving RPA excitation energies
through its connection to equation-of-motion coupled-cluster theory,"
J. Chem. Phys. 153, 234101 (2020) [arXiv:2008.00336], H2O row (1B1/1A2/
21A1/1B2 states, modified aug-cc-pVDZ* basis = CC-PVDZ+DIFFUSE(3) below).

Covers the xvee excited-state pipeline (xjoda...xringcc/xringcl/xvee),
which is a SEPARATE code path from the xringcl ground-state Lambda solver
fixed in tests/ringcl/ -- this suite only pins its current, verified-good
behavior, it doesn't exercise anything touched by that fix.

Provenance: ZMAT/GENBAS recovered 2026-08-14 from Varun Rishi's own
HiPerGator directory (co-author on the paper), the original run that
produced this data:
    /blue/bartlett/varunrishi/scratch/check/eom_drccd_dec9_18/sauer_1995/H2O
(ZMAT geometry/basis there cites Christiansen/Koch/Jorgensen, JCP 105, 1451
(1996) and Christiansen/Bak/Koch/Sauer, CPL 284, 47 (1998) for the
geometry; that directory's own name "sauer_1995" and its extraction script
extractExcitationEnergies.py -- replicated below via the same "grep lines
containing 'Singlet'" logic -- confirm this is the actual source of the
paper's Table I numbers, not a lookalike). Rerun on the current (2026-08-14
ringcc/ringcl-fixed) build and found to reproduce Varun's original archived
output BIT-FOR-BIT for all four ZMAT variants below -- i.e. this week's
ringcc/ringcl fixes had zero effect on this pathway, as expected. Three of
the four Table I columns (EOM(Sf)-rCCD, EOM(Sf)-drCCD, EOM(S)-drCCD) also
match the PUBLISHED paper to the two decimals it prints.

RESOLVED DISCREPANCY (2026-08-14, kept for history -- do not re-open
without new evidence): H2O's EOM(S)-rCCD (Fock-dressed, WITH exchange)
column reproduces Varun's original 2018 archive exactly, but is
consistently ~0.28-0.33 eV LOWER than the published Table I value on every
H2O state (e.g. 1B1: 10.007 eV here vs 10.33 eV published). Investigated
2026-08-14 by extending this suite to N2, Ne, CH2, BH (below) using the
same recovered archive: N2 and Ne's EOM(S)-rCCD columns match the
published paper closely (~0.01 eV), so this is NOT a universal EOM(S)-rCCD
formula bug -- it is H2O-specific. Confirmed the gap is present bit-for-
bit in Varun's own original 2018 run (`eom_s_rccd.out`), predating any
2026 work, with no post-processing between his raw output and the
published-comparison numbers (his `extractExcitationEnergies.py` just
`grep`s "Singlet" lines verbatim). User's working conclusion: likely a
manuscript-side transcription error by Varun when Table I was assembled,
not an ACES II code bug -- pending his own confirmation, not indepedently
verified. See `project_jcp2020_eom_rccd_recovery.md` (memory) Checkpoints
4-5 for the full comparison tables and reasoning.

KNOWN SEPARATE, UNRELATED ISSUES (not fixed here, do not let them block
this suite's PASS status): (1) N2's EOM summary output mislabels some
near-degenerate singlet/triplet state pairs (swapped Singlet/Triplet tags,
or unresolved "??????") -- the underlying eigenvalues are correct, only
the multiplicity classification is wrong, so this suite's N2 REFERENCE
values below reflect the (buggy) CURRENT labeling, matching this suite's
existing philosophy of pinning current-build behavior rather than a
presumed-correct one. (2) CH2 and BH's RCCD/UHF ground state (`xringcc`)
fails to converge in 50 cycles (oscillating residual, not decreasing) --
blocks `eomsfrccd`/`eomsrccd` for both molecules entirely (never reaches
`xringcl`/`xvee`). Likely the same already-known DRCCD/UHF bug tracked
separately in the scrnc/GW-manuscript project memory. Their ZMATs are
present under `cases/ch2_eomsfrccd/`, `cases/ch2_eomsrccd/`,
`cases/bh_eomsfrccd/`, `cases/bh_eomsrccd/` (recovered from the same
Varun Rishi archive) but deliberately NOT registered in REFERENCE below --
add them once the convergence bug is fixed.

Usage:
    module load intel/2025.1.0 openmpi/5.0.7   # must be loaded first
    python3 tests/jcp2020_eom_rccd/run_tests.py [--jobs N] [--case NAME ...]
"""
import argparse
import concurrent.futures
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

SUITE_DIR = os.path.dirname(os.path.abspath(__file__))
CASES_DIR = os.path.join(SUITE_DIR, "cases")
# One GENBAS for all suites (ACES_PROGRAM/basis).
GENBAS = os.path.join(SUITE_DIR, os.pardir, os.pardir, os.pardir, "basis", "GENBAS")

TOLERANCE_EV = 1e-4

PIPELINE = ["xjoda", "xvmol", "xvmol2ja", "xvprops", "xvscf", "xvtran",
            "xintprc", "xringcc", "xringcl", "xvee"]

_EXTRA_LD_PATH = "/apps/compilers/intel/2025/1.0.666/compiler/2025.1/lib"

_SINGLET_RE = re.compile(
    r"^\s*(A1|A2|B1|B2|Ag|B1g|B2g|B3g|Au|B1u|B2u|B3u)"
    r"\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+Singlet\s*$"
)

# Reference: the "Singlet"-tagged rows of each case's EOM summary table,
# as produced by the current (2026-08-14) build -- confirmed bit-identical
# to Varun Rishi's original archived run (see module docstring).
# 2026-11-01: the doubles corrections that changed by more than 1e-5 eV
# when vee started to biorthonormalize the left and right eigenvectors
# (vee/biorth_rcc_vecs.F) were updated. Before that, the correction of a
# root was scaled by <L|R> (up to 1.6e-3 eV off in N2), and in a block of
# degenerate roots it depended on how DGEEV happened to pick the vectors,
# which made the Ne Ag rows change with the compiler (0.02-0.24 eV). See
# docs/programmer/issues/2026-10-29-ring-ne-ag-correction.md.
# 2026-11-03: every doubles correction (and Exc. Eng.) regenerated after
# vee/form_dbls_correctns_2rcc.F started to use each root's own omega in
# the denominator (omega_k - D2), as in eqs. 37-41 of the paper, instead
# of omega of the first root of the irrep (314 rows changed; no base value
# changed). The first root of each irrep is unchanged. The published
# Table 4 (D) values for N2 and Ne were made with the old omega; see
# docs/programmer/issues/2026-11-01-notes-for-upstream.md. Rows of high
# roots that sit near a doubles pole (omega_k ~ D2) now carry large
# corrections of either sign; that is the CIS(D)-type intruder behavior,
# not noise, but those rows are the most sensitive ones in the suite.
# N2, Ne, CH2, BH added 2026-08-14 from the same Varun Rishi archive (see
# module docstring). N2's Singlet/Triplet labeling issue (see docstring)
# means these lists reflect the current labeling. ch2_eomsfrccd,
# ch2_eomsrccd, bh_eomsfrccd and bh_eomsrccd are deliberately NOT
# included: their RCCD/UHF ground state does not converge (ZMATs kept
# under cases/).
# (irrep, base_value_eV, doubles_correction_eV, exc_energy_with_D_eV)
REFERENCE = {
    "h2o_eomsfrccd": [
        ("A1", 10.954420, -3.356265, 7.598155),
        ("A1", 12.830829, -3.393248, 9.437581),
        ("B1", 8.627536, -3.456701, 5.170835),
        ("B2", 12.605486, -3.377088, 9.228398),
        ("A2", 10.315738, -3.452184, 6.863553),
    ],
    "h2o_eomsrccd": [
        ("A1", 12.337861, -3.505824, 8.832037),
        ("A1", 14.217284, -3.523466, 10.693819),
        ("B1", 10.007320, -3.622603, 6.384717),
        ("B2", 13.995828, -3.567981, 10.427847),
        ("A2", 11.718210, -3.679667, 8.038542),
    ],
    "h2o_eomsfdrccd": [
        ("A1", 17.008345, -4.293124, 12.715221),
        ("A1", 20.242838, -4.190455, 16.052383),
        ("A1", 21.208067, -3.504780, 17.703287),
        ("A1", 21.750608, -4.353393, 17.397214),
        ("B1", 14.926477, -4.592328, 10.334149),
        ("B1", 19.049443, -4.169264, 14.880178),
        ("B1", 21.708752, -3.697878, 18.010874),
        ("B2", 17.535867, -4.260578, 13.275289),
        ("B2", 20.673838, -3.401086, 17.272753),
        ("B2", 22.929277, -3.918013, 19.011264),
        ("A2", 15.474952, -4.584040, 10.890911),
        ("A2", 20.826857, -4.232607, 16.594250),
    ],
    "h2o_eomsdrccd": [
        ("A1", 17.582465, -4.372148, 13.210317),
        ("A1", 20.848223, -4.274551, 16.573673),
        ("A1", 21.807949, -3.569379, 18.238569),
        ("A1", 22.337549, -4.443325, 17.894223),
        ("B1", 15.484287, -4.670568, 10.813719),
        ("B1", 19.633844, -4.247783, 15.386061),
        ("B1", 22.334563, -3.774423, 18.560140),
        ("B2", 18.106337, -4.338738, 13.767599),
        ("B2", 21.278075, -3.462629, 17.815446),
        ("B2", 23.511683, -3.992882, 19.518800),
        ("A2", 16.028104, -4.661479, 11.366625),
        ("A2", 21.387651, -4.304583, 17.083068),
    ],
    "n2_eomsfrccd": [
        ("Ag", 27.321470, -2.562695, 24.758775),
        ("Ag", 28.709746, -2.144312, 26.565434),
        ("B3u", 15.787973, -4.070344, 11.717629),
        ("B3u", 28.370025, -1.472330, 26.897696),
        ("B3u", 29.157583, -2.996557, 26.161025),
        ("B2u", 15.787973, -4.070344, 11.717629),
        ("B2u", 28.370025, -1.472330, 26.897696),
        ("B2u", 29.157583, -2.996557, 26.161025),
        ("B1g", 28.709746, -2.144312, 26.565434),
        ("B1g", 51.720898, -32.011038, 19.709859),
        ("B1u", 8.955490, -1.845626, 7.109865),
        ("B1u", 15.619110, -2.051796, 13.567314),
        ("B1u", 23.831848, -3.090990, 20.740858),
        ("B2g", 9.825162, -3.166627, 6.658535),
        ("B2g", 23.653722, -1.992529, 21.661193),
        ("B3g", 9.825162, -3.166627, 6.658535),
        ("B3g", 23.653722, -1.992529, 21.661193),
        ("Au", 8.955490, -1.845626, 7.109865),
    ],
    "n2_eomsrccd": [
        ("Ag", 28.770206, -2.709538, 26.060669),
        ("Ag", 31.123703, -2.269301, 28.854402),
        ("B3u", 18.448374, -4.298551, 14.149824),
        ("B3u", 30.402217, -2.803635, 27.598583),
        ("B3u", 30.894367, -1.809280, 29.085086),
        ("B2u", 18.448374, -4.298551, 14.149824),
        ("B2u", 30.402217, -2.803635, 27.598583),
        ("B2u", 30.894367, -1.809280, 29.085086),
        ("B1g", 31.123703, -2.269301, 28.854402),
        ("B1g", 53.129268, -2.994953, 50.134315),
        ("B1u", 12.859721, -1.917223, 10.942498),
        ("B1u", 19.199484, -2.119573, 17.079911),
        ("B1u", 25.263156, -3.173472, 22.089684),
        ("B2g", 12.744858, -3.342339, 9.402519),
        ("B2g", 26.158156, -2.090482, 24.067673),
        ("B3g", 12.744858, -3.342339, 9.402519),
        ("B3g", 26.158156, -2.090482, 24.067673),
        ("Au", 12.859721, -1.917223, 10.942498),
    ],
    "n2_eomsfdrccd": [
        ("Ag", 38.003303, -7.327844, 30.675458),
        ("Ag", 40.660063, -3.689593, 36.970471),
        ("Ag", 41.172367, -3.742288, 37.430078),
        ("Ag", 44.487004, -4.869785, 39.617219),
        ("Ag", 48.503447, -6.417820, 42.085627),
        ("Ag", 55.870471, -5.437332, 50.433140),
        ("B3u", 28.631537, -6.391077, 22.240460),
        ("B3u", 39.763654, -2.249099, 37.514555),
        ("B3u", 41.268096, -6.650830, 34.617266),
        ("B3u", 44.706753, -3.584363, 41.122390),
        ("B3u", 50.799677, -7.983921, 42.815757),
        ("B3u", 64.875881, -16.533286, 48.342595),
        ("B2u", 28.631537, -6.391077, 22.240460),
        ("B2u", 39.763654, -2.249099, 37.514555),
        ("B2u", 41.268096, -6.650830, 34.617266),
        ("B2u", 44.706753, -3.584363, 41.122390),
        ("B2u", 50.799677, -7.983921, 42.815757),
        ("B2u", 64.875881, -16.533286, 48.342595),
        ("B1g", 40.295821, -3.045091, 37.250730),
        ("B1g", 41.172367, -3.742288, 37.430078),
        ("B1g", 65.360863, -5.942933, 59.417930),
        ("B1g", 67.894327, -2.583343, 65.310984),
        ("B1g", 69.109108, -12.870278, 56.238830),
        ("B1g", 85.001796, -14.695183, 70.306613),
        ("B1u", 22.735565, -2.130530, 20.605035),
        ("B1u", 28.852097, -2.366443, 26.485654),
        ("B1u", 34.415276, -5.193567, 29.221709),
        ("B1u", 44.105401, -5.286050, 38.819351),
        ("B1u", 46.011451, -4.136875, 41.874577),
        ("B1u", 47.469159, -6.290609, 41.178550),
        ("B2g", 23.429307, -4.541947, 18.887360),
        ("B2g", 33.362924, -2.520967, 30.841957),
        ("B2g", 45.062055, -12.288605, 32.773450),
        ("B2g", 45.932150, -16.381314, 29.550837),
        ("B2g", 46.682549, -11.983956, 34.698592),
        ("B2g", 49.168568, -2.527981, 46.640587),
        ("B3g", 23.429307, -4.541947, 18.887360),
        ("B3g", 33.362924, -2.520967, 30.841957),
        ("B3g", 45.062055, -12.288605, 32.773450),
        ("B3g", 45.932150, -16.381314, 29.550837),
        ("B3g", 46.682549, -11.983956, 34.698592),
        ("B3g", 49.168568, -2.527981, 46.640587),
        ("Au", 21.415936, -1.929757, 19.486180),
        ("Au", 22.735565, -2.130530, 20.605035),
        ("Au", 45.209016, -3.733725, 41.475291),
        ("Au", 46.011451, -4.136875, 41.874577),
        ("Au", 69.993765, -8.512042, 61.481723),
        ("Au", 80.553300, -9.608395, 70.944904),
    ],
    "n2_eomsdrccd": [
        ("Ag", 38.570267, -7.503165, 31.067102),
        ("Ag", 41.313324, -3.848009, 37.465314),
        ("Ag", 42.004768, -4.889890, 37.114877),
        ("Ag", 45.194044, -5.396664, 39.797380),
        ("Ag", 49.244848, -6.677084, 42.567764),
        ("Ag", 56.623981, -5.829369, 50.794612),
        ("B3u", 29.512302, -6.622458, 22.889844),
        ("B3u", 40.594232, -2.320766, 38.273466),
        ("B3u", 41.891982, -7.376750, 34.515231),
        ("B3u", 45.589604, -3.743705, 41.845899),
        ("B3u", 51.409127, -8.234994, 43.174133),
        ("B3u", 65.353975, -23.483938, 41.870037),
        ("B2u", 29.512302, -6.622458, 22.889844),
        ("B2u", 40.594232, -2.320766, 38.273466),
        ("B2u", 41.891982, -7.376750, 34.515231),
        ("B2u", 45.589604, -3.743705, 41.845899),
        ("B2u", 51.409127, -8.234994, 43.174133),
        ("B2u", 65.353975, -23.483938, 41.870037),
        ("B1g", 41.128256, -3.398178, 37.730077),
        ("B1g", 42.004768, -4.889890, 37.114877),
        ("B1g", 66.002192, -6.832722, 59.169469),
        ("B1g", 68.741434, -7.496922, 61.244512),
        ("B1g", 69.950303, 25.083281, 95.033585),
        ("B1g", 85.512377, -15.780440, 69.731937),
        ("B1u", 23.887088, -2.160931, 21.726157),
        ("B1u", 29.913882, -2.404656, 27.509226),
        ("B1u", 35.026538, -5.300087, 29.726451),
        ("B1u", 44.673005, -5.411901, 39.261104),
        ("B1u", 46.878073, -4.323185, 42.554888),
        ("B1u", 48.270026, -6.729121, 41.540905),
        ("B2g", 24.384212, -4.663812, 19.720400),
        ("B2g", 34.181325, -2.569582, 31.611742),
        ("B2g", 45.650827, -13.072402, 32.578425),
        ("B2g", 46.672992, -29.539824, 17.133168),
        ("B2g", 47.362217, 77.636004, 124.998221),
        ("B2g", 50.148972, -3.023856, 47.125116),
        ("B3g", 24.384212, -4.663812, 19.720400),
        ("B3g", 34.181325, -2.569582, 31.611742),
        ("B3g", 45.650827, -13.072402, 32.578425),
        ("B3g", 46.672992, -29.539824, 17.133168),
        ("B3g", 47.362217, 77.636004, 124.998221),
        ("B3g", 50.148972, -3.023856, 47.125116),
        ("Au", 22.579574, -1.952362, 20.627212),
        ("Au", 23.887088, -2.160931, 21.726157),
        ("Au", 46.064912, -3.880447, 42.184465),
        ("Au", 46.878073, -4.323185, 42.554888),
        ("Au", 70.572592, -15.072540, 55.500052),
        ("Au", 81.126084, -10.374559, 70.751525),
    ],
    "ne_eomsfrccd": [
        ("Ag", 19.943569, -4.273877, 15.669693),
        ("Ag", 19.943569, -4.273877, 15.669693),
        ("B3u", 18.084460, -4.089664, 13.994796),
        ("B3u", 49.556465, -5.722165, 43.834300),
        ("B2u", 18.084460, -4.089664, 13.994796),
        ("B2u", 49.556465, -5.722165, 43.834300),
        ("B1g", 19.943569, -4.273877, 15.669693),
        ("B1u", 18.084460, -4.089664, 13.994796),
        ("B1u", 49.556465, -5.722165, 43.834300),
        ("B2g", 19.943569, -4.273877, 15.669693),
        ("B3g", 19.943569, -4.273877, 15.669693),
    ],
    "ne_eomsrccd": [
        ("Ag", 21.026280, -4.340012, 16.686268),
        ("Ag", 21.026280, -4.340012, 16.686268),
        ("B3u", 19.177309, -4.149346, 15.027963),
        ("B3u", 50.071509, -6.553852, 43.517657),
        ("B2u", 19.177309, -4.149346, 15.027963),
        ("B2u", 50.071509, -6.553852, 43.517657),
        ("B1g", 21.026280, -4.340012, 16.686268),
        ("B1u", 19.177309, -4.149346, 15.027963),
        ("B1u", 50.071509, -6.553852, 43.517657),
        ("B2g", 21.026280, -4.340012, 16.686268),
        ("B3g", 21.026280, -4.340012, 16.686268),
    ],
    "ne_eomsfdrccd": [
        ("Ag", 25.030008, -4.625626, 20.404381),
        ("Ag", 25.030008, -4.625626, 20.404381),
        ("Ag", 25.384598, -4.949522, 20.435077),
        ("Ag", 55.079028, -5.507676, 49.571352),
        ("Ag", 71.840892, -1.299051, 70.541840),
        ("B3u", 25.514205, -4.569976, 20.944229),
        ("B3u", 54.611536, -5.649528, 48.962009),
        ("B3u", 88.198177, -2.277810, 85.920367),
        ("B3u", 102.052077, -2.119276, 99.932801),
        ("B3u", 164.024774, -1.101723, 162.923051),
        ("B2u", 25.514205, -4.569976, 20.944229),
        ("B2u", 54.611536, -5.649528, 48.962009),
        ("B2u", 88.198177, -2.277810, 85.920367),
        ("B2u", 102.052077, -2.119276, 99.932801),
        ("B2u", 164.024774, -1.101723, 162.923051),
        ("B1g", 24.991224, -4.602894, 20.388330),
        ("B1g", 25.030008, -4.625626, 20.404381),
        ("B1g", 69.609172, -1.262315, 68.346857),
        ("B1g", 71.840892, -1.299051, 70.541840),
        ("B1g", 199.239213, 1.279173, 200.518386),
        ("B1u", 25.514205, -4.569976, 20.944229),
        ("B1u", 54.611536, -5.649528, 48.962009),
        ("B1u", 88.198177, -2.277810, 85.920367),
        ("B1u", 102.052077, -2.119276, 99.932801),
        ("B1u", 164.024774, -1.101723, 162.923051),
        ("B2g", 24.991224, -4.602894, 20.388330),
        ("B2g", 25.030008, -4.625626, 20.404381),
        ("B2g", 69.609172, -1.262315, 68.346857),
        ("B2g", 71.840892, -1.299051, 70.541840),
        ("B2g", 199.239213, 1.279173, 200.518386),
        ("B3g", 24.991224, -4.602894, 20.388330),
        ("B3g", 25.030008, -4.625626, 20.404381),
        ("B3g", 69.609172, -1.262315, 68.346857),
        ("B3g", 71.840892, -1.299051, 70.541840),
        ("B3g", 199.239213, 1.279173, 200.518386),
    ],
    "ne_eomsdrccd": [
        ("Ag", 25.590032, -4.664711, 20.925321),
        ("Ag", 25.590032, -4.664711, 20.925321),
        ("Ag", 25.941437, -4.996764, 20.944673),
        ("Ag", 55.373786, -5.541688, 49.832098),
        ("Ag", 72.687813, -1.309464, 71.378349),
        ("B3u", 26.080400, -4.608932, 21.471468),
        ("B3u", 54.897676, -5.678597, 49.219078),
        ("B3u", 88.942412, -2.326246, 86.616167),
        ("B3u", 102.632253, -2.132594, 100.499659),
        ("B3u", 164.740265, -1.170890, 163.569375),
        ("B2u", 26.080400, -4.608932, 21.471468),
        ("B2u", 54.897676, -5.678597, 49.219078),
        ("B2u", 88.942412, -2.326246, 86.616167),
        ("B2u", 102.632253, -2.132594, 100.499659),
        ("B2u", 164.740265, -1.170890, 163.569375),
        ("B1g", 25.551712, -4.641473, 20.910239),
        ("B1g", 25.590032, -4.664711, 20.925321),
        ("B1g", 70.455694, -1.271357, 69.184337),
        ("B1g", 72.687813, -1.309464, 71.378349),
        ("B1g", 199.680895, 0.226525, 199.907420),
        ("B1u", 26.080400, -4.608932, 21.471468),
        ("B1u", 54.897676, -5.678597, 49.219078),
        ("B1u", 88.942412, -2.326246, 86.616167),
        ("B1u", 102.632253, -2.132594, 100.499659),
        ("B1u", 164.740265, -1.170890, 163.569375),
        ("B2g", 25.551712, -4.641473, 20.910239),
        ("B2g", 25.590032, -4.664711, 20.925321),
        ("B2g", 70.455694, -1.271357, 69.184337),
        ("B2g", 72.687813, -1.309464, 71.378349),
        ("B2g", 199.680895, 0.226525, 199.907420),
        ("B3g", 25.551712, -4.641473, 20.910239),
        ("B3g", 25.590032, -4.664711, 20.925321),
        ("B3g", 70.455694, -1.271357, 69.184337),
        ("B3g", 72.687813, -1.309464, 71.378349),
        ("B3g", 199.680895, 0.226525, 199.907420),
    ],
    "ch2_eomsfdrccd": [
        ("A1", 11.513752, -1.840017, 9.673735),
        ("A1", 12.831347, -1.964111, 10.867236),
        ("A1", 17.055453, -1.710121, 15.345331),
        ("A1", 19.085036, -2.159289, 16.925747),
        ("A1", 24.526772, -2.067102, 22.459670),
        ("A1", 25.068516, -2.423108, 22.645408),
        ("A1", 26.470242, -2.530760, 23.939482),
        ("A1", 31.843926, -1.961435, 29.882491),
        ("B1", 13.744809, -1.600541, 12.144268),
        ("B1", 28.360733, -2.804925, 25.555808),
        ("B1", 31.675855, -2.496268, 29.179587),
        ("B1", 45.482729, -3.335743, 42.146986),
        ("B1", 47.878942, -1.842472, 46.036469),
        ("B1", 53.340187, -1.716137, 51.624051),
        ("B1", 61.906111, -2.187537, 59.718574),
        ("B1", 66.580107, -6.847486, 59.732621),
        ("B2", 12.365332, -1.906897, 10.458435),
        ("B2", 16.039345, -1.647700, 14.391645),
        ("B2", 17.412190, -1.668007, 15.744183),
        ("B2", 19.505745, -2.075594, 17.430151),
        ("B2", 24.269488, -1.880806, 22.388681),
        ("B2", 25.913563, -2.573302, 23.340261),
        ("B2", 28.964353, -2.013103, 26.951250),
        ("B2", 32.753844, -2.981225, 29.772619),
        ("A2", 17.749626, -1.868016, 15.881610),
        ("A2", 35.763381, -1.212157, 34.551224),
        ("A2", 44.203749, -7.986479, 36.217270),
        ("A2", 51.231681, -1.381864, 49.849816),
        ("A2", 59.491848, -3.053377, 56.438472),
        ("A2", 66.334660, -7.984981, 58.349679),
        ("A2", 71.873234, 0.624785, 72.498020),
        ("A2", 79.997889, -27.634678, 52.363211),
    ],
    "ch2_eomsdrccd": [
        ("A1", 11.923914, -1.866434, 10.057480),
        ("A1", 13.246221, -1.992251, 11.253970),
        ("A1", 17.485731, -1.734435, 15.751297),
        ("A1", 19.553736, -2.191840, 17.361896),
        ("A1", 24.997186, -2.150415, 22.846771),
        ("A1", 25.451245, -2.427560, 23.023685),
        ("A1", 26.834458, -2.568512, 24.265946),
        ("A1", 32.318577, -2.021556, 30.297021),
        ("B1", 14.261260, -1.616990, 12.644270),
        ("B1", 28.819233, -2.874320, 25.944913),
        ("B1", 32.113882, -2.645518, 29.468364),
        ("B1", 45.876635, -3.774440, 42.102195),
        ("B1", 48.359879, -2.026018, 46.333861),
        ("B1", 53.822261, -1.456976, 52.365285),
        ("B1", 62.334210, -2.332112, 60.002098),
        ("B1", 67.005759, 1.396339, 68.402098),
        ("B2", 12.777379, -1.933886, 10.843493),
        ("B2", 16.468971, -1.669933, 14.799038),
        ("B2", 17.846719, -1.691029, 16.155689),
        ("B2", 19.992889, -2.109948, 17.882941),
        ("B2", 24.752194, -1.909965, 22.842229),
        ("B2", 26.272506, -2.614888, 23.657619),
        ("B2", 29.437112, -1.724453, 27.712659),
        ("B2", 33.190730, -3.043658, 30.147072),
        ("A2", 18.287264, -1.887025, 16.400239),
        ("A2", 36.219435, -2.389719, 33.829715),
        ("A2", 44.677834, -1.644987, 43.032847),
        ("A2", 51.727538, -1.506517, 50.221021),
        ("A2", 59.914795, -4.289308, 55.625487),
        ("A2", 66.757276, 11.285475, 78.042752),
        ("A2", 72.313781, -1.713686, 70.600095),
        ("A2", 80.366851, 9.090327, 89.457178),
    ],
    "bh_eomsfdrccd": [
        ("A1", 9.666880, -1.013520, 8.653360),
        ("A1", 9.800959, -0.983012, 8.817947),
        ("A1", 10.459147, -1.034761, 9.424386),
        ("A1", 11.187447, -1.093535, 10.093912),
        ("A1", 11.595249, -0.998907, 10.596343),
        ("A1", 11.728388, -1.013810, 10.714578),
        ("B1", 9.769641, -0.993065, 8.776575),
        ("B1", 10.795332, -1.130684, 9.664648),
        ("B1", 11.683946, -1.024538, 10.659408),
        ("B1", 12.467622, -1.111272, 11.356350),
        ("B1", 13.984088, -1.243861, 12.740228),
        ("B1", 17.994275, -1.252493, 16.741782),
        ("B2", 9.769641, -0.993065, 8.776575),
        ("B2", 10.795332, -1.130684, 9.664648),
        ("B2", 11.683946, -1.024538, 10.659408),
        ("B2", 12.467622, -1.111272, 11.356350),
        ("B2", 13.984088, -1.243861, 12.740228),
        ("B2", 17.994275, -1.252493, 16.741782),
        ("A2", 11.595249, -0.998907, 10.596343),
        ("A2", 18.401156, -1.460472, 16.940684),
        ("A2", 19.814572, -1.254643, 18.559928),
        ("A2", 26.530594, -1.527320, 25.003274),
        ("A2", 39.881972, 0.074328, 39.956301),
        ("A2", 48.257949, -1.107831, 47.150119),
    ],
    "bh_eomsdrccd": [
        ("A1", 10.014796, -1.030485, 8.984311),
        ("A1", 10.148635, -0.999648, 9.148986),
        ("A1", 10.808169, -1.052180, 9.755989),
        ("A1", 11.537284, -1.112488, 10.424796),
        ("A1", 11.942859, -1.016032, 10.926826),
        ("A1", 12.076368, -1.031318, 11.045050),
        ("B1", 10.117714, -1.009571, 9.108142),
        ("B1", 11.157845, -1.146221, 10.011624),
        ("B1", 12.033845, -1.042320, 10.991525),
        ("B1", 12.830186, -1.128312, 11.701874),
        ("B1", 14.378153, -1.267693, 13.110460),
        ("B1", 18.331755, -1.275862, 17.055893),
        ("B2", 10.117714, -1.009571, 9.108142),
        ("B2", 11.157845, -1.146221, 10.011624),
        ("B2", 12.033845, -1.042320, 10.991525),
        ("B2", 12.830186, -1.128312, 11.701874),
        ("B2", 14.378153, -1.267693, 13.110460),
        ("B2", 18.331755, -1.275862, 17.055893),
        ("A2", 11.942859, -1.016032, 10.926826),
        ("A2", 18.750185, -1.534836, 17.215349),
        ("A2", 20.151401, -1.278282, 18.873119),
        ("A2", 26.868839, -1.578205, 25.290634),
        ("A2", 40.235579, 1.031270, 41.266849),
        ("A2", 48.600887, -0.483841, 48.117045),
    ],
}


def _parse_singlets(out_text):
    rows = []
    for line in out_text.splitlines():
        m = _SINGLET_RE.match(line)
        if m:
            rows.append((m.group(1), float(m.group(2)), float(m.group(3)),
                         float(m.group(4))))
    return rows


def run_case(case_name, env):
    case_dir = os.path.join(CASES_DIR, case_name)
    zmat = os.path.join(case_dir, "ZMAT")
    if not os.path.isfile(zmat):
        return case_name, "ERROR", 0.0, f"no ZMAT found in {case_dir}"

    t0 = time.time()
    with tempfile.TemporaryDirectory(prefix=f"acesii_jcp2020_{case_name}_") as scratch:
        shutil.copy(GENBAS, os.path.join(scratch, "GENBAS"))
        shutil.copy(zmat, os.path.join(scratch, "ZMAT"))

        out_chunks = []
        for exe in PIPELINE:
            # errors="replace": RCCD/UHF xvee output (N2/Ne eomsfrccd/
            # eomsrccd) has been observed to contain a handful of stray
            # NUL/0x01 bytes; harmless to the parsed Singlet table but
            # fatal to strict text decoding.
            try:
                run = subprocess.run(
                    [exe], cwd=scratch, env=env,
                    capture_output=True, text=True, errors="replace",
                    timeout=300,
                )
            except subprocess.TimeoutExpired:
                return case_name, "ERROR", time.time() - t0, f"{exe} timed out (300s)"
            out_chunks.append(run.stdout)
            if run.returncode != 0:
                tail = "\n".join(run.stdout.splitlines()[-15:])
                return case_name, "ERROR", time.time() - t0, f"{exe} exited {run.returncode}\n{tail}"

        got = _parse_singlets("\n".join(out_chunks))
        want = REFERENCE[case_name]
        elapsed = time.time() - t0

        if len(got) != len(want):
            return case_name, "FAIL", elapsed, (
                f"expected {len(want)} Singlet rows, got {len(got)}\n"
                f"got: {got}"
            )

        mismatches = []
        for i, (g, w) in enumerate(zip(got, want)):
            if g[0] != w[0] or any(abs(g[j] - w[j]) > TOLERANCE_EV for j in (1, 2, 3)):
                mismatches.append(f"row {i}: got {g}, want {w}")
        if mismatches:
            return case_name, "FAIL", elapsed, "\n".join(mismatches)

        return case_name, "PASS", elapsed, f"{len(got)} Singlet state(s) matched (tol {TOLERANCE_EV:.0e} eV)"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--jobs", type=int, default=4,
                         help="cases to run concurrently (default 4)")
    parser.add_argument("--case", action="append", dest="cases",
                         help="run only these case names (repeatable); "
                              "default is all cases")
    parser.add_argument("--xaces2", action="store_true",
                         help="run each case with xaces2 instead of the "
                              "program chain above, against the same "
                              "references (checks the driver's ring-CC route)")
    args = parser.parse_args()
    if args.xaces2:
        PIPELINE[:] = ["xaces2"]

    if shutil.which("xvee") is None:
        sys.exit("xvee not found on PATH -- module load intel/2025.1.0 "
                  "openmpi/5.0.7 first, and make sure ~/Develop/bin is on PATH")

    env = os.environ.copy()
    ld = env.get("LD_LIBRARY_PATH", "")
    if _EXTRA_LD_PATH not in ld:
        env["LD_LIBRARY_PATH"] = _EXTRA_LD_PATH + (":" + ld if ld else "")

    all_case_names = sorted(REFERENCE.keys())
    case_names = args.cases if args.cases else all_case_names
    unknown = set(case_names) - set(all_case_names)
    if unknown:
        sys.exit(f"unknown case(s): {', '.join(sorted(unknown))}")

    print(f"Running {len(case_names)} case(s), --jobs {args.jobs}\n")
    print(f"{'CASE':<20} {'TIME':>8}  {'STATUS':<6}  DETAIL")
    print("-" * 90)

    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futures = {ex.submit(run_case, name, env): name for name in case_names}
        for fut in concurrent.futures.as_completed(futures):
            results.append(fut.result())

    results.sort(key=lambda r: r[0])
    n_pass = n_fail = n_error = 0
    total_time = 0.0
    for name, status, elapsed, detail in results:
        total_time += elapsed
        if status == "PASS":
            n_pass += 1
        elif status == "FAIL":
            n_fail += 1
        else:
            n_error += 1
        detail_first_line = detail.splitlines()[0] if detail else ""
        print(f"{name:<20} {elapsed:>6.2f}s  {status:<6}  {detail_first_line}")
        if status in ("FAIL", "ERROR") and "\n" in detail:
            for line in detail.splitlines()[1:]:
                print(f"    {line}")

    print("-" * 90)
    print(f"{len(case_names)} case(s), {n_pass} passed, {n_fail} failed, "
          f"{n_error} errored, {total_time:.1f}s total CPU-time-of-runs")

    sys.exit(0 if (n_fail == 0 and n_error == 0) else 1)


if __name__ == "__main__":
    main()
