#!/usr/bin/env python3
"""Native RPA (aces2/rpa/, executable xrpa) regression suite.

This is a DIFFERENT RPA implementation from the ring-CC-based "EOM(SF)-
RCCD (RPA)" path covered by tests/jcp2020_eom_rccd/ -- this one is the
standalone module in aces2/rpa/ (drive_rpa.F, do_rhf_rpa.F, etc.),
CALC=SCF + EXCITE=RPA, running its own short pipeline (xjoda, xvmol,
xvmol2ja, xvscf, xvtran, xintprc, xrpa -- no xvprops/xvee/xringcc/
xringcl). It prints excitation energies in a different, per-state format
("The excitation energy = X eV" / "The spin-state = Singlet"), grouped
under "The excited states of symmetry block: N" headers (N is a raw
symmetry-block index, not mapped to an irrep letter here).

Provenance: ZMATs recovered 2026-08-14 from the same Varun Rishi archive
used for tests/jcp2020_eom_rccd/ (Rishi/Perera/Bartlett JCP 153, 234101
(2020) Table I RPA column):
    /blue/bartlett/varunrishi/scratch/check/eom_drccd_dec9_18/sauer_1995/{H2O,N2,Ne_ccpvdz,CH2,BH}/rpa/ZMAT
Varun's own rpa.csv extraction for ALL FIVE molecules came back EMPTY
(header only) -- his extractExcitationEnergies.py regex expects a single
line containing "  Singlet" with the value as line.split()[1], but this
module's output splits the energy and spin-state across two separate
lines, so his script never matched anything here. Not a code bug, just an
extraction-script mismatch with this module's differently-shaped output;
this suite parses it directly instead.

N2, CH2, BH ran cleanly and are pinned in REFERENCE below (values
confirmed self-consistent within this session, not yet cross-checked
against the published paper's RPA column -- do that before treating them
as physics-validated, not just regression-safe).

h2o_rpa and ne_rpa (DROPMO=1) were broken until 2026-10-28: vtran ignored
DROPMO for CALC=SCF even with EXCITE set, leaving SYMPOPVA one orbital
larger than NVRTORB, and xintprc (DS4AA) ran off its arrays (H2O) or
xrpa's NTO step read ORBMAP_A with the full basis length (Ne). Fixed in
vtran.F and rpa_dens_main.F (no NTOs with frozen orbitals). Their
singlets agree to 1e-6 eV with the base column of tests/jcp2020_eom_rccd's
h2o_eomsfrccd / ne_eomsfrccd (EOM(SF)-RCCD through xringcc/xringcl/xvee,
the same RPA by another code path); the references below are from this
module's own run.

Usage:
    module load intel/2025.1.0 openmpi/5.0.7   # must be loaded first
    python3 tests/rpa/run_tests.py [--jobs N] [--case NAME ...]
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

PIPELINE = ["xjoda", "xvmol", "xvmol2ja", "xvscf", "xvtran", "xintprc", "xrpa"]

_EXTRA_LD_PATH = "/apps/compilers/intel/2025/1.0.666/compiler/2025.1/lib"

_BLOCK_RE = re.compile(r"The excited states of symmetry block:\s*(\d+)")
_ENERGY_RE = re.compile(r"The excitation energy\s*=\s*([-\d.]+)eV")
_SPIN_RE = re.compile(r"The spin-state\s*=\s*(\w+)")

# Reference: (symmetry_block, excitation_energy_eV, spin_state), in the
# order printed, as produced by the current (2026-08-14) build.
REFERENCE = {
    "n2_rpa": [
        (1, 27.321219, "Singlet"), (1, 28.709359, "Singlet"),
        (1, 29.732999, "Singlet"), (1, 32.952161, "Singlet"),
        (1, 36.181166, "Singlet"), (1, 40.907820, "Singlet"),
        (2, 15.779776, "Singlet"), (2, 28.369658, "Singlet"),
        (2, 29.157194, "Singlet"), (2, 32.658648, "Singlet"),
        (2, 39.449668, "Singlet"), (2, 51.065325, "Singlet"),
        (3, 15.779776, "Singlet"), (3, 28.369658, "Singlet"),
        (3, 29.157194, "Singlet"), (3, 32.658648, "Singlet"),
        (3, 39.449668, "Singlet"), (3, 51.065325, "Singlet"),
        (4, 28.654641, "Singlet"), (4, 28.709359, "Singlet"),
        (4, 51.720746, "Singlet"), (4, 54.283146, "Singlet"),
        (4, 54.372059, "Singlet"), (4, 71.826047, "Singlet"),
        (5, 8.954847, "Singlet"), (5, 15.617690, "Singlet"),
        (5, 23.831108, "Singlet"), (5, 32.894226, "Singlet"),
        (5, 33.832220, "Singlet"), (5, 34.644455, "Singlet"),
        (6, 9.823709, "Singlet"), (6, 23.651360, "Singlet"),
        (6, 31.423590, "Singlet"), (6, 31.925437, "Singlet"),
        (6, 34.374780, "Singlet"), (6, 36.756482, "Singlet"),
        (7, 9.823709, "Singlet"), (7, 23.651360, "Singlet"),
        (7, 31.423590, "Singlet"), (7, 31.925437, "Singlet"),
        (7, 34.374780, "Singlet"), (7, 36.756482, "Singlet"),
        (8, 8.085119, "Singlet"), (8, 8.954847, "Singlet"),
        (8, 33.819991, "Singlet"), (8, 33.832220, "Singlet"),
        (8, 56.203859, "Singlet"), (8, 67.700865, "Singlet"),
    ],
    "ch2_rpa": [
        (1, 6.943607, "Singlet"), (1, 8.844964, "Singlet"),
        (1, 11.805845, "Singlet"), (1, 12.768728, "Singlet"),
        (1, 15.874301, "Singlet"), (1, 20.274647, "Singlet"),
        (1, 21.438322, "Singlet"), (1, 22.366784, "Singlet"),
        (2, 1.302540, "Singlet"), (2, 16.356653, "Singlet"),
        (2, 19.285264, "Singlet"), (2, 32.938242, "Singlet"),
        (2, 35.295975, "Singlet"), (2, 39.073030, "Singlet"),
        (2, 47.939140, "Singlet"), (2, 55.695811, "Singlet"),
        (3, 8.102412, "Singlet"), (3, 11.127202, "Singlet"),
        (3, 11.228635, "Singlet"), (3, 13.350147, "Singlet"),
        (3, 16.142558, "Singlet"), (3, 19.255665, "Singlet"),
        (3, 21.513376, "Singlet"), (3, 23.838209, "Singlet"),
        (4, 6.037432, "Singlet"), (4, 24.322564, "Singlet"),
        (4, 32.418451, "Singlet"), (4, 38.857335, "Singlet"),
        (4, 45.864261, "Singlet"), (4, 54.979447, "Singlet"),
        (4, 59.060517, "Singlet"), (4, 67.578930, "Singlet"),
    ],
    "bh_rpa": [
        (1, 6.363443, "Singlet"), (1, 7.380581, "Singlet"),
        (1, 8.064313, "Singlet"), (1, 8.069552, "Singlet"),
        (1, 8.078888, "Singlet"), (1, 8.417586, "Singlet"),
        (2, 2.644557, "Singlet"), (2, 7.366206, "Singlet"),
        (2, 8.108484, "Singlet"), (2, 8.420138, "Singlet"),
        (2, 9.598671, "Singlet"), (2, 11.292905, "Singlet"),
        (3, 2.644557, "Singlet"), (3, 7.366206, "Singlet"),
        (3, 8.108484, "Singlet"), (3, 8.420138, "Singlet"),
        (3, 9.598671, "Singlet"), (3, 11.292905, "Singlet"),
        (4, 8.078888, "Singlet"), (4, 12.059077, "Singlet"),
        (4, 16.295507, "Singlet"), (4, 20.202479, "Singlet"),
        (4, 29.662806, "Singlet"), (4, 37.834261, "Singlet"),
    ],
    "h2o_rpa": [
        (1, 10.954420, "Singlet"), (1, 12.830829, "Singlet"),
        (1, 14.750139, "Singlet"), (1, 15.767523, "Singlet"),
        (2, 8.627536, "Singlet"), (2, 12.187912, "Singlet"),
        (2, 13.604293, "Singlet"), (3, 12.605485, "Singlet"),
        (3, 14.350839, "Singlet"), (3, 16.060907, "Singlet"),
        (4, 10.315737, "Singlet"), (4, 13.414325, "Singlet"),
    ],
    "ne_rpa": [
        (1, 19.943569, "Singlet"), (1, 19.943569, "Singlet"),
        (1, 20.350752, "Singlet"), (1, 47.760815, "Singlet"),
        (1, 50.452063, "Singlet"), (2, 18.084459, "Singlet"),
        (2, 49.556464, "Singlet"), (2, 67.177969, "Singlet"),
        (2, 80.516181, "Singlet"), (2, 139.608156, "Singlet"),
        (3, 18.084459, "Singlet"), (3, 49.556464, "Singlet"),
        (3, 67.177969, "Singlet"), (3, 80.516181, "Singlet"),
        (3, 139.608156, "Singlet"), (4, 19.943569, "Singlet"),
        (4, 19.972230, "Singlet"), (4, 49.931572, "Singlet"),
        (4, 50.452063, "Singlet"), (4, 171.308768, "Singlet"),
        (5, 18.084459, "Singlet"), (5, 49.556464, "Singlet"),
        (5, 67.177969, "Singlet"), (5, 80.516181, "Singlet"),
        (5, 139.608156, "Singlet"), (6, 19.943569, "Singlet"),
        (6, 19.972230, "Singlet"), (6, 49.931572, "Singlet"),
        (6, 50.452063, "Singlet"), (6, 171.308768, "Singlet"),
        (7, 19.943569, "Singlet"), (7, 19.972230, "Singlet"),
        (7, 49.931572, "Singlet"), (7, 50.452063, "Singlet"),
        (7, 171.308768, "Singlet"), (8, 139.608156, "Singlet"),
        (8, 139.608156, "Singlet"), (8, 139.770652, "Singlet"),
        (8, 50.452063, "Singlet"), (8, 171.308768, "Singlet"),
    ],
}

# Cases with a known, real, unfixed crash. Run every time (never silently
# skipped) so the suite catches it either way: reports XFAIL while the
# crash still reproduces exactly as documented, but flags loudly (UPASS or
# a differently-shaped ERROR) the moment that stops being true -- either
# because someone fixed it (move the case into REFERENCE and drop it from
# here) or because something regressed differently (needs investigation).
#
EXPECTED_CRASH = {}


def _parse_states(out_text):
    lines = out_text.splitlines()
    block = None
    rows = []
    for i, line in enumerate(lines):
        m = _BLOCK_RE.search(line)
        if m:
            block = int(m.group(1))
        m2 = _ENERGY_RE.search(line)
        if m2:
            energy = float(m2.group(1))
            spin = None
            for j in range(i + 1, min(i + 4, len(lines))):
                m3 = _SPIN_RE.search(lines[j])
                if m3:
                    spin = m3.group(1)
                    break
            rows.append((block, energy, spin))
    return rows


def run_case(case_name, env):
    case_dir = os.path.join(CASES_DIR, case_name)
    zmat = os.path.join(case_dir, "ZMAT")
    if not os.path.isfile(zmat):
        return case_name, "ERROR", 0.0, f"no ZMAT found in {case_dir}"

    expected = EXPECTED_CRASH.get(case_name)

    t0 = time.time()
    with tempfile.TemporaryDirectory(prefix=f"acesii_rpa_{case_name}_") as scratch:
        shutil.copy(GENBAS, os.path.join(scratch, "GENBAS"))
        shutil.copy(zmat, os.path.join(scratch, "ZMAT"))

        out_chunks = []
        for exe in PIPELINE:
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
                elapsed = time.time() - t0
                tail = "\n".join((run.stdout + run.stderr).splitlines()[-15:])
                if (expected and exe == expected["exe"]
                        and run.returncode == expected["returncode"]
                        and expected["signature"] in (run.stdout + run.stderr)):
                    return case_name, "XFAIL", elapsed, (
                        f"known bug reproduced: {exe} exited "
                        f"{run.returncode} ('{expected['signature']}')"
                    )
                return case_name, "ERROR", elapsed, f"{exe} exited {run.returncode}\n{tail}"

        elapsed = time.time() - t0
        if expected:
            return case_name, "UPASS", elapsed, (
                f"expected known crash in {expected['exe']} but the full "
                f"pipeline completed -- bug may be fixed; check output, "
                f"pin REFERENCE, and remove this case from EXPECTED_CRASH"
            )

        got = _parse_states("\n".join(out_chunks))
        want = REFERENCE[case_name]

        if len(got) != len(want):
            return case_name, "FAIL", elapsed, (
                f"expected {len(want)} state(s), got {len(got)}\n"
                f"got: {got}"
            )

        mismatches = []
        for i, (g, w) in enumerate(zip(got, want)):
            if g[0] != w[0] or g[2] != w[2] or abs(g[1] - w[1]) > TOLERANCE_EV:
                mismatches.append(f"row {i}: got {g}, want {w}")
        if mismatches:
            return case_name, "FAIL", elapsed, "\n".join(mismatches)

        return case_name, "PASS", elapsed, f"{len(got)} state(s) matched (tol {TOLERANCE_EV:.0e} eV)"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--jobs", type=int, default=4,
                         help="cases to run concurrently (default 4)")
    parser.add_argument("--case", action="append", dest="cases",
                         help="run only these case names (repeatable); "
                              "default is all cases")
    args = parser.parse_args()

    if shutil.which("xrpa") is None:
        sys.exit("xrpa not found on PATH -- module load intel/2025.1.0 "
                  "openmpi/5.0.7 first, and make sure ~/Develop/bin is on PATH")

    env = os.environ.copy()
    ld = env.get("LD_LIBRARY_PATH", "")
    if _EXTRA_LD_PATH not in ld:
        env["LD_LIBRARY_PATH"] = _EXTRA_LD_PATH + (":" + ld if ld else "")

    all_case_names = sorted(set(REFERENCE) | set(EXPECTED_CRASH))
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
    n_pass = n_fail = n_error = n_xfail = n_upass = 0
    total_time = 0.0
    for name, status, elapsed, detail in results:
        total_time += elapsed
        if status == "PASS":
            n_pass += 1
        elif status == "FAIL":
            n_fail += 1
        elif status == "XFAIL":
            n_xfail += 1
        elif status == "UPASS":
            n_upass += 1
        else:
            n_error += 1
        detail_first_line = detail.splitlines()[0] if detail else ""
        print(f"{name:<20} {elapsed:>6.2f}s  {status:<6}  {detail_first_line}")
        if status in ("FAIL", "ERROR", "UPASS") and "\n" in detail:
            for line in detail.splitlines()[1:]:
                print(f"    {line}")

    print("-" * 90)
    print(f"{len(case_names)} case(s), {n_pass} passed, {n_fail} failed, "
          f"{n_error} errored, {n_xfail} known-crash (xfail), "
          f"{n_upass} unexpectedly-fixed (upass), "
          f"{total_time:.1f}s total CPU-time-of-runs")

    sys.exit(0 if (n_fail == 0 and n_error == 0 and n_upass == 0) else 1)


if __name__ == "__main__":
    main()
