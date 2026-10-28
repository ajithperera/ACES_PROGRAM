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
# (irrep, base_value_eV, doubles_correction_eV, exc_energy_with_D_eV)
REFERENCE = {
    "h2o_eomsfrccd": [
        ("A1", 10.954420, -3.285141, 7.669279),
        ("A1", 12.830829, -3.186331, 9.644498),
        ("B1", 8.627536, -3.403779, 5.223757),
        ("B2", 12.605486, -3.305809, 9.299677),
        ("A2", 10.315738, -3.422541, 6.893197),
    ],
    "h2o_eomsrccd": [
        ("A1", 12.337861, -3.437060, 8.900801),
        ("A1", 14.217284, -3.307989, 10.909295),
        ("B1", 10.007320, -3.568195, 6.439125),
        ("B2", 13.995828, -3.504418, 10.491411),
        ("A2", 11.718210, -3.651711, 8.066499),
    ],
    "h2o_eomsfdrccd": [
        ("A1", 17.008345, -4.293118, 12.715226),
        ("A1", 20.242838, -3.828783, 16.414055),
        ("A1", 21.208067, -3.134362, 18.073706),
        ("A1", 21.750608, -3.787776, 17.962832),
        ("B1", 14.926477, -4.592328, 10.334149),
        ("B1", 19.049443, -3.746733, 15.302710),
        ("B1", 21.708752, -3.147022, 18.561730),
        ("B2", 17.535867, -4.260578, 13.275289),
        ("B2", 20.673838, -3.125580, 17.548259),
        ("B2", 22.929277, -3.413511, 19.515765),
        ("A2", 15.474952, -4.584040, 10.890911),
        ("A2", 20.826857, -3.691448, 17.135409),
    ],
    "h2o_eomsdrccd": [
        ("A1", 17.582465, -4.372143, 13.210322),
        ("A1", 20.848223, -3.892845, 16.955378),
        ("A1", 21.807949, -3.181302, 18.626647),
        ("A1", 22.337549, -3.850622, 18.486926),
        ("B1", 15.484287, -4.670568, 10.813720),
        ("B1", 19.633844, -3.805191, 15.828653),
        ("B1", 22.334563, -3.195678, 19.138885),
        ("B2", 18.106337, -4.338737, 13.767600),
        ("B2", 21.278075, -3.172533, 18.105542),
        ("B2", 23.511683, -3.466704, 20.044979),
        ("A2", 16.028104, -4.661478, 11.366626),
        ("A2", 21.387651, -3.742449, 17.645203),
    ],
    # N2, Ne, CH2, BH added 2026-08-14 from the same Varun Rishi archive
    # (see module docstring) while investigating the H2O EOM(S)-rCCD gap.
    # N2's Singlet/Triplet labeling bug (see docstring) means these
    # REFERENCE lists reflect current (imperfect) labeling, not a
    # presumed-correct one -- consistent with this suite's existing
    # regression-pinning philosophy.
    "n2_eomsfrccd": [
        ("Ag", 27.321470, -2.404701, 24.916768),
        ("Ag", 28.709746, -1.958657, 26.751090),
        ("B3u", 15.787973, -3.739206, 12.048767),
        ("B3u", 28.370025, -1.150552, 27.219473),
        ("B3u", 29.157583, -2.130759, 27.026824),
        ("B2u", 15.787973, -3.739206, 12.048767),
        ("B2u", 28.370025, -1.150552, 27.219473),
        ("B2u", 29.157583, -2.130759, 27.026824),
        ("B1g", 28.709746, -2.101374, 26.608372),
        ("B1g", 51.720898, -2.080125, 49.640772),
        ("B1u", 8.955490, -1.737318, 7.218172),
        ("B1u", 15.619110, -1.777276, 13.841833),
        ("B1u", 23.831848, -2.223538, 21.608310),
        ("B2g", 9.825162, -3.058975, 6.766187),
        ("B2g", 23.653722, -1.582947, 22.070776),
        ("B3g", 9.825162, -3.058975, 6.766187),
        ("B3g", 23.653722, -1.582947, 22.070776),
        ("Au", 8.955490, -1.784342, 7.171148),
    ],
    "n2_eomsrccd": [
        ("Ag", 28.770206, -2.529684, 26.240523),
        ("Ag", 31.123703, -2.013343, 29.110360),
        ("B3u", 18.448374, -3.919586, 14.528788),
        ("B3u", 30.402217, -2.044968, 28.357249),
        ("B3u", 30.894367, -1.344962, 29.549405),
        ("B2u", 18.448374, -3.919586, 14.528788),
        ("B2u", 30.402217, -2.044968, 28.357249),
        ("B2u", 30.894367, -1.344962, 29.549405),
        ("B1g", 31.123703, -2.218918, 28.904784),
        ("B1g", 53.129268, -2.167959, 50.961309),
        ("B1u", 12.859721, -1.798803, 11.060917),
        ("B1u", 19.199484, -1.827225, 17.372259),
        ("B1u", 25.263156, -2.339674, 22.923482),
        ("B2g", 12.744858, -3.218043, 9.526815),
        ("B2g", 26.158156, -1.647621, 24.510535),
        ("B3g", 12.744858, -3.218043, 9.526815),
        ("B3g", 26.158156, -1.647621, 24.510535),
        ("Au", 12.859721, -1.850150, 11.009571),
    ],
    "n2_eomsfdrccd": [
        ("Ag", 38.003303, -7.327840, 30.675463),
        ("Ag", 40.660063, -3.298721, 37.361342),
        ("Ag", 41.172367, -2.833357, 38.339010),
        ("Ag", 44.487004, -5.378699, 39.108305),
        ("Ag", 48.503447, -5.304740, 43.198707),
        ("Ag", 55.870471, -3.732318, 52.138153),
        ("B3u", 28.631537, -6.390969, 22.240568),
        ("B3u", 39.763654, -1.688311, 38.075343),
        ("B3u", 41.268096, -4.091848, 37.176248),
        ("B3u", 44.706753, -2.474821, 42.231932),
        ("B3u", 50.799677, -4.133065, 46.666613),
        ("B3u", 64.875881, -2.601559, 62.274322),
        ("B2u", 28.631537, -6.390969, 22.240568),
        ("B2u", 39.763654, -1.688311, 38.075343),
        ("B2u", 41.268096, -4.091848, 37.176248),
        ("B2u", 44.706753, -2.474821, 42.231932),
        ("B2u", 50.799677, -4.133065, 46.666613),
        ("B2u", 64.875881, -2.601559, 62.274322),
        ("B1g", 40.295821, -3.045091, 37.250730),
        ("B1g", 41.172367, -3.311598, 37.860769),
        ("B1g", 65.360863, -3.712105, 61.648757),
        ("B1g", 67.894327, -4.916882, 62.977445),
        ("B1g", 69.109108, -5.459945, 63.649163),
        ("B1g", 85.001796, -4.116200, 80.885595),
        ("B1u", 22.735565, -2.130529, 20.605036),
        ("B1u", 28.852097, -2.133753, 26.718344),
        ("B1u", 34.415276, -3.890563, 30.524713),
        ("B1u", 44.105401, -2.999871, 41.105531),
        ("B1u", 46.011451, -2.213253, 43.798198),
        ("B1u", 47.469159, -3.113980, 44.355179),
        ("B2g", 23.429307, -4.541929, 18.887378),
        ("B2g", 33.362924, -2.067011, 31.295913),
        ("B2g", 45.062055, -5.536443, 39.525612),
        ("B2g", 45.932150, -3.393206, 42.538944),
        ("B2g", 46.682549, -2.324916, 44.357633),
        ("B2g", 49.168568, -2.062301, 47.106267),
        ("B3g", 23.429307, -4.541929, 18.887378),
        ("B3g", 33.362924, -2.067011, 31.295913),
        ("B3g", 45.062055, -5.536443, 39.525612),
        ("B3g", 45.932150, -3.393206, 42.538944),
        ("B3g", 46.682549, -2.324916, 44.357633),
        ("B3g", 49.168568, -2.062301, 47.106267),
        ("Au", 21.415936, -1.929757, 19.486180),
        ("Au", 22.735565, -2.087389, 20.648176),
        ("Au", 45.209016, -2.055941, 43.153075),
        ("Au", 46.011451, -2.162214, 43.849238),
        ("Au", 69.993765, -2.828057, 67.165708),
        ("Au", 80.553300, -2.128944, 78.424356),
    ],
    "n2_eomsdrccd": [
        ("Ag", 38.570267, -7.503157, 31.067111),
        ("Ag", 41.313324, -3.357783, 37.955540),
        ("Ag", 42.004768, -2.919773, 39.084995),
        ("Ag", 45.194044, -5.523971, 39.670073),
        ("Ag", 49.244848, -5.435887, 43.808960),
        ("Ag", 56.623981, -3.780673, 52.843308),
        ("B3u", 29.512302, -6.622344, 22.889958),
        ("B3u", 40.594232, -1.716319, 38.877914),
        ("B3u", 41.891982, -4.196048, 37.695933),
        ("B3u", 45.589604, -2.518098, 43.071505),
        ("B3u", 51.409127, -4.218520, 47.190607),
        ("B3u", 65.353975, -2.641841, 62.712135),
        ("B2u", 29.512302, -6.622344, 22.889958),
        ("B2u", 40.594232, -1.716319, 38.877914),
        ("B2u", 41.891982, -4.196048, 37.695933),
        ("B2u", 45.589604, -2.518098, 43.071505),
        ("B2u", 51.409127, -4.218520, 47.190607),
        ("B2u", 65.353975, -2.641841, 62.712135),
        ("B1g", 41.128256, -3.398178, 37.730077),
        ("B1g", 42.004768, -3.717667, 38.287101),
        ("B1g", 66.002192, -3.963678, 62.038513),
        ("B1g", 68.741434, -6.540515, 62.200919),
        ("B1g", 69.950303, -7.340810, 62.609493),
        ("B1g", 85.512377, -4.209224, 81.303153),
        ("B1u", 23.887088, -2.160926, 21.726162),
        ("B1u", 29.913882, -2.163990, 27.749892),
        ("B1u", 35.026538, -3.991921, 31.034617),
        ("B1u", 44.673005, -3.065854, 41.607151),
        ("B1u", 46.878073, -2.269663, 44.608409),
        ("B1u", 48.270026, -3.313411, 44.956615),
        ("B2g", 24.384212, -4.663774, 19.720438),
        ("B2g", 34.181325, -2.102141, 32.079184),
        ("B2g", 45.650827, -5.593755, 40.057072),
        ("B2g", 46.672992, -3.346647, 43.326345),
        ("B2g", 47.362217, -2.541333, 44.820884),
        ("B2g", 50.148972, -2.079941, 48.069031),
        ("B3g", 24.384212, -4.663774, 19.720438),
        ("B3g", 34.181325, -2.102141, 32.079184),
        ("B3g", 45.650827, -5.593755, 40.057072),
        ("B3g", 46.672992, -3.346647, 43.326345),
        ("B3g", 47.362217, -2.541333, 44.820884),
        ("B3g", 50.148972, -2.079941, 48.069031),
        ("Au", 22.579574, -1.952357, 20.627218),
        ("Au", 23.887088, -2.116322, 21.770766),
        ("Au", 46.064912, -2.105738, 43.959174),
        ("Au", 46.878073, -2.216468, 44.661604),
        ("Au", 70.572592, -2.870635, 67.701958),
        ("Au", 81.126084, -2.164964, 78.961120),
    ],
    "ne_eomsfrccd": [
        ("Ag", 19.943569, -4.172985, 15.770585),
        ("Ag", 19.943569, -4.172985, 15.770585),
        ("B3u", 18.084460, -4.074080, 14.010380),
        ("B3u", 49.556465, -3.829686, 45.726779),
        ("B2u", 18.084460, -4.074080, 14.010380),
        ("B2u", 49.556465, -3.829686, 45.726779),
        ("B1g", 19.943569, -4.268984, 15.674586),
        ("B1u", 18.084460, -4.074080, 14.010380),
        ("B1u", 49.556465, -3.829686, 45.726779),
        ("B2g", 19.943569, -4.268984, 15.674586),
        ("B3g", 19.943569, -4.268984, 15.674586),
    ],
    "ne_eomsrccd": [
        ("Ag", 21.026280, -4.298957, 16.727322),
        ("Ag", 21.026280, -4.298957, 16.727322),
        ("B3u", 19.177309, -4.133467, 15.043842),
        ("B3u", 50.071509, -3.868722, 46.202787),
        ("B2u", 19.177309, -4.133467, 15.043842),
        ("B2u", 50.071509, -3.868722, 46.202787),
        ("B1g", 21.026280, -4.335160, 16.691120),
        ("B1u", 19.177309, -4.133467, 15.043842),
        ("B1u", 50.071509, -3.868722, 46.202787),
        ("B2g", 21.026280, -4.335160, 16.691120),
        ("B3g", 21.026280, -4.335160, 16.691120),
    ],
    "ne_eomsfdrccd": [
        ("Ag", 25.030008, -4.606569, 20.423438),
        ("Ag", 25.030008, -4.606569, 20.423438),
        ("Ag", 25.384598, -4.920036, 20.464563),
        ("Ag", 55.079028, -4.069861, 51.009167),
        ("Ag", 71.840892, -0.682103, 71.158789),
        ("B3u", 25.514205, -4.569976, 20.944230),
        ("B3u", 54.611536, -4.128032, 50.483504),
        ("B3u", 88.198177, -1.356077, 86.842100),
        ("B3u", 102.052077, -1.257053, 100.795024),
        ("B3u", 164.024774, -0.517069, 163.507705),
        ("B2u", 25.514205, -4.569976, 20.944230),
        ("B2u", 54.611536, -4.128032, 50.483504),
        ("B2u", 88.198177, -1.356077, 86.842100),
        ("B2u", 102.052077, -1.257053, 100.795024),
        ("B2u", 164.024774, -0.517069, 163.507705),
        ("B1g", 24.991224, -4.602894, 20.388330),
        ("B1g", 25.030008, -4.622987, 20.407020),
        ("B1g", 69.609172, -0.926065, 68.683107),
        ("B1g", 71.840892, -0.924126, 70.916766),
        ("B1g", 199.239213, -0.850471, 198.388742),
        ("B1u", 25.514205, -4.569976, 20.944230),
        ("B1u", 54.611536, -4.128032, 50.483504),
        ("B1u", 88.198177, -1.356077, 86.842100),
        ("B1u", 102.052077, -1.257053, 100.795024),
        ("B1u", 164.024774, -0.517069, 163.507705),
        ("B2g", 24.991224, -4.602894, 20.388330),
        ("B2g", 25.030008, -4.622987, 20.407020),
        ("B2g", 69.609172, -0.926065, 68.683107),
        ("B2g", 71.840892, -0.924126, 70.916766),
        ("B2g", 199.239213, -0.850471, 198.388742),
        ("B3g", 24.991224, -4.602894, 20.388330),
        ("B3g", 25.030008, -4.622987, 20.407020),
        ("B3g", 69.609172, -0.926065, 68.683107),
        ("B3g", 71.840892, -0.924126, 70.916766),
        ("B3g", 199.239213, -0.850471, 198.388742),
    ],
    "ne_eomsdrccd": [
        ("Ag", 25.590032, -4.500696, 21.089335),
        ("Ag", 25.590032, -4.500696, 21.089335),
        ("Ag", 25.941437, -4.966773, 20.974664),
        ("Ag", 55.373786, -4.091760, 51.282026),
        ("Ag", 72.687813, -0.881774, 71.806039),
        ("B3u", 26.080400, -4.608932, 21.471469),
        ("B3u", 54.897676, -4.150484, 50.747192),
        ("B3u", 88.942412, -1.361667, 87.580746),
        ("B3u", 102.632253, -1.258415, 101.373838),
        ("B3u", 164.740265, -0.518088, 164.222177),
        ("B2u", 26.080400, -4.608932, 21.471469),
        ("B2u", 54.897676, -4.150484, 50.747192),
        ("B2u", 88.942412, -1.361667, 87.580746),
        ("B2u", 102.632253, -1.258415, 101.373838),
        ("B2u", 164.740265, -0.518088, 164.222177),
        ("B1g", 25.551712, -4.641473, 20.910239),
        ("B1g", 25.590032, -4.662047, 20.927985),
        ("B1g", 70.455694, -0.928296, 69.527398),
        ("B1g", 72.687813, -0.926564, 71.761248),
        ("B1g", 199.680895, -0.852091, 198.828805),
        ("B1u", 26.080400, -4.608932, 21.471469),
        ("B1u", 54.897676, -4.150484, 50.747192),
        ("B1u", 88.942412, -1.361667, 87.580746),
        ("B1u", 102.632253, -1.258415, 101.373838),
        ("B1u", 164.740265, -0.518088, 164.222177),
        ("B2g", 25.551712, -4.641473, 20.910239),
        ("B2g", 25.590032, -4.662047, 20.927985),
        ("B2g", 70.455694, -0.928296, 69.527398),
        ("B2g", 72.687813, -0.926564, 71.761248),
        ("B2g", 199.680895, -0.852091, 198.828805),
        ("B3g", 25.551712, -4.641473, 20.910239),
        ("B3g", 25.590032, -4.662047, 20.927985),
        ("B3g", 70.455694, -0.928296, 69.527398),
        ("B3g", 72.687813, -0.926564, 71.761248),
        ("B3g", 199.680895, -0.852091, 198.828805),
    ],
    # ch2_eomsfrccd / ch2_eomsrccd / bh_eomsfrccd / bh_eomsrccd
    # deliberately NOT included -- RCCD/UHF ground state doesn't converge
    # for these two molecules (see docstring). ZMATs kept under cases/
    # for when that's fixed.
    "ch2_eomsfdrccd": [
        ("A1", 11.513752, -1.840015, 9.673736),
        ("A1", 12.831347, -1.881607, 10.949739),
        ("A1", 17.055453, -1.453217, 15.602236),
        ("A1", 19.085036, -1.766654, 17.318383),
        ("A1", 24.526772, -1.483006, 23.043766),
        ("A1", 25.068516, -1.605298, 23.463217),
        ("A1", 26.470242, -1.594387, 24.875855),
        ("A1", 31.843926, -1.242304, 30.601622),
        ("B1", 13.744809, -1.600513, 12.144296),
        ("B1", 28.360733, -1.822241, 26.538493),
        ("B1", 31.675855, -1.322329, 30.353526),
        ("B1", 45.482729, -1.483146, 43.999583),
        ("B1", 47.878942, -0.989837, 46.889105),
        ("B1", 53.340187, -0.886161, 52.454026),
        ("B1", 61.906111, -0.804412, 61.101699),
        ("B1", 66.580107, -1.031450, 65.548657),
        ("B2", 12.365332, -1.906897, 10.458435),
        ("B2", 16.039345, -1.480872, 14.558473),
        ("B2", 17.412190, -1.443704, 15.968487),
        ("B2", 19.505745, -1.714570, 17.791175),
        ("B2", 24.269488, -1.401767, 22.867720),
        ("B2", 25.913563, -1.678716, 24.234847),
        ("B2", 28.964353, -1.246719, 27.717633),
        ("B2", 32.753844, -1.655913, 31.097931),
        ("A2", 17.749626, -1.868014, 15.881612),
        ("A2", 35.763381, -1.662349, 34.101032),
        ("A2", 44.203749, -1.317154, 42.886596),
        ("A2", 51.231681, -0.988541, 50.243140),
        ("A2", 59.491848, -1.069005, 58.422843),
        ("A2", 66.334660, -1.164116, 65.170544),
        ("A2", 71.873234, -0.970804, 70.902430),
        ("A2", 79.997889, -1.002404, 78.995485),
    ],
    "ch2_eomsdrccd": [
        ("A1", 11.923914, -1.866433, 10.057481),
        ("A1", 13.246221, -1.906753, 11.339468),
        ("A1", 17.485731, -1.469080, 16.016651),
        ("A1", 19.553736, -1.783946, 17.769790),
        ("A1", 24.997186, -1.519513, 23.477673),
        ("A1", 25.451245, -1.598116, 23.853128),
        ("A1", 26.834458, -1.607275, 25.227183),
        ("A1", 32.318577, -1.255844, 31.062733),
        ("B1", 14.261260, -1.616966, 12.644295),
        ("B1", 28.819233, -1.840215, 26.979018),
        ("B1", 32.113882, -1.339905, 30.773977),
        ("B1", 45.876635, -1.500380, 44.376254),
        ("B1", 48.359879, -0.997128, 47.362751),
        ("B1", 53.822261, -0.894290, 52.927972),
        ("B1", 62.334210, -0.808752, 61.525458),
        ("B1", 67.005759, -1.038122, 65.967637),
        ("B2", 12.777379, -1.933886, 10.843493),
        ("B2", 16.468971, -1.497540, 14.971432),
        ("B2", 17.846719, -1.459281, 16.387437),
        ("B2", 19.992889, -1.732916, 18.259973),
        ("B2", 24.752194, -1.413655, 23.338539),
        ("B2", 26.272506, -1.695187, 24.577319),
        ("B2", 29.437112, -1.254428, 28.182684),
        ("B2", 33.190730, -1.666960, 31.523770),
        ("A2", 18.287264, -1.887023, 16.400241),
        ("A2", 36.219435, -1.686616, 34.532819),
        ("A2", 44.677834, -1.333148, 43.344687),
        ("A2", 51.727538, -0.998120, 50.729418),
        ("A2", 59.914795, -1.076394, 58.838402),
        ("A2", 66.757276, -1.175683, 65.581593),
        ("A2", 72.313781, -0.978578, 71.335203),
        ("A2", 80.366851, -1.009029, 79.357822),
    ],
    "bh_eomsfdrccd": [
        ("A1", 9.666880, -1.013520, 8.653360),
        ("A1", 9.800959, -0.976769, 8.824190),
        ("A1", 10.459147, -0.997156, 9.461991),
        ("A1", 11.187447, -1.018382, 10.169065),
        ("A1", 11.595249, -0.915022, 10.680227),
        ("A1", 11.728388, -0.922484, 10.805904),
        ("B1", 9.769641, -0.993065, 8.776576),
        ("B1", 10.795332, -1.080792, 9.714540),
        ("B1", 11.683946, -0.938437, 10.745508),
        ("B1", 12.467622, -0.995332, 11.472290),
        ("B1", 13.984088, -1.072245, 12.911843),
        ("B1", 17.994275, -0.891804, 17.102471),
        ("B2", 9.769641, -0.993065, 8.776576),
        ("B2", 10.795332, -1.080792, 9.714540),
        ("B2", 11.683946, -0.938437, 10.745508),
        ("B2", 12.467622, -0.995332, 11.472290),
        ("B2", 13.984088, -1.072245, 12.911843),
        ("B2", 17.994275, -0.891804, 17.102471),
        ("A2", 11.595249, -0.998907, 10.596343),
        ("A2", 18.401156, -0.942521, 17.458635),
        ("A2", 19.814572, -0.892036, 18.922536),
        ("A2", 26.530594, -0.808854, 25.721740),
        ("A2", 39.881972, -0.918235, 38.963737),
        ("A2", 48.257949, -0.802020, 47.455929),
    ],
    "bh_eomsdrccd": [
        ("A1", 10.014796, -1.030485, 8.984311),
        ("A1", 10.148635, -0.993164, 9.155471),
        ("A1", 10.808169, -1.013032, 9.795138),
        ("A1", 11.537284, -1.034246, 10.503038),
        ("A1", 11.942859, -0.928935, 11.013923),
        ("A1", 12.076368, -0.936457, 11.139910),
        ("B1", 10.117714, -1.009571, 9.108143),
        ("B1", 11.157845, -1.093682, 10.064163),
        ("B1", 12.033845, -0.952711, 11.081134),
        ("B1", 12.830186, -1.007460, 11.822726),
        ("B1", 14.378153, -1.088719, 13.289434),
        ("B1", 18.331755, -0.901916, 17.429840),
        ("B2", 10.117714, -1.009571, 9.108143),
        ("B2", 11.157845, -1.093682, 10.064163),
        ("B2", 12.033845, -0.952711, 11.081134),
        ("B2", 12.830186, -1.007460, 11.822726),
        ("B2", 14.378153, -1.088719, 13.289434),
        ("B2", 18.331755, -0.901916, 17.429840),
        ("A2", 11.942859, -1.016032, 10.926827),
        ("A2", 18.750185, -0.956506, 17.793679),
        ("A2", 20.151401, -0.902462, 19.248939),
        ("A2", 26.868839, -0.816619, 26.052220),
        ("A2", 40.235579, -0.933159, 39.302420),
        ("A2", 48.600887, -0.809544, 47.791343),
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
    args = parser.parse_args()

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
