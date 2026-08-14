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

KNOWN OPEN DISCREPANCY (not a regression -- do not "fix" this suite to
hide it): EOM(S)-rCCD (Fock-dressed, WITH exchange -- the h2o_eomsrccd
case) reproduces Varun's original archive exactly, but is consistently
~0.3 eV LOWER than the published Table I value on every state (e.g. 1B1:
10.007 eV here vs 10.33 eV published). Since the archive and this rebuild
agree bit-for-bit, this predates and is unrelated to any 2026-08 work --
it lives in whatever Fock-matrix-dressing (similarity transformation) code
underlies EOM(S) specifically. Flagged for a separate investigation.

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
GENBAS = os.path.join(SUITE_DIR, "GENBAS")

TOLERANCE_EV = 1e-4

PIPELINE = ["xjoda", "xvmol", "xvmol2ja", "xvprops", "xvscf", "xvtran",
            "xintprc", "xringcc", "xringcl", "xvee"]

_EXTRA_LD_PATH = "/apps/compilers/intel/2025/1.0.666/compiler/2025.1/lib"

_SINGLET_RE = re.compile(
    r"^\s*(A1|A2|B1|B2)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+Singlet\s*$"
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
            try:
                run = subprocess.run(
                    [exe], cwd=scratch, env=env,
                    capture_output=True, text=True, timeout=300,
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
                  "openmpi/5.0.7 first, and make sure aces2/bin is on PATH")

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
