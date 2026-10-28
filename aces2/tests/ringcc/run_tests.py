#!/usr/bin/env python3
"""ACES II ring-CC (ringcc/ringcl) regression suite.

Covers CALC=RCCD and CALC=DRCCD (both REF=RHF and REF=UHF) -- the separate
ring-CC executables (xringcc/xringcl), NOT the standard CC driver. `xaces2`
cannot drive these: it routes CALC=RCCD/DRCCD through the ordinary xvcc CC
program instead, which has its own dead/broken RCCD/DRCCD branch (crashes
with "@GETLST: Error reading list (1,14)") -- confirmed empirically
2026-08-13, see project memory. So unlike tests/fast/run_tests.py (which
just runs `xaces2` once per case), this driver manually chains the same
executable sequence /blue/bartlett/perera/Testing/Testdir_rcc/run_rcc uses:
xjoda, xvmol, xvmol2ja, xvprops, xvscf, xvtran, xintprc, xringcc, xringcl.

Reference TOTENERG values (test_results) for drccd_rhf, drccd_uhf and
rccd_uhf come from re-running each case's exact ZMAT keyword block against
its archived legacy output in /blue/bartlett/perera/Testing/Testdir_rcc/
(drccd_rhf_sto3g.out, drccd_uhf_sto3g.out, rccd_uhf_sto3g.out) -- same
geometry (H2O, R=0.957 Ang, A=104.5 deg) throughout.

rccd_rhf (CALC=RCCD, REF=RHF, cc-pVDZ) is NOT the archive value. Its archive
(rccd_rhf_adz.out, -76.3039658461) was made with the old singlet/triplet
weighting of the RHF RCCD energy (0.75/0.25) in aces2/ringcc/rcc_cmpeng.F.
d80f2742 (2026-08-19) changed it to the equal weighting (0.25/0.25) of
Rishi, Perera and Bartlett, JCP 153, 234101 (2020), Eq. 49, after which RHF
and UHF RCCD agree for a closed shell. The reference is that value,
-76.1821730567, and rccd_uhf_ccpvdz (the same job with REF=UHF, through the
UHF code path) must give it too.

Usage:
    module load intel/2025.1.0 openmpi/5.0.7   # must be loaded first
    python3 tests/ringcc/run_tests.py [--jobs N] [--case NAME ...]
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

RINGCC_DIR = os.path.dirname(os.path.abspath(__file__))
CASES_DIR = os.path.join(RINGCC_DIR, "cases")
# One GENBAS for all suites (ACES_PROGRAM/basis).
GENBAS = os.path.join(RINGCC_DIR, os.pardir, os.pardir, os.pardir, "basis", "GENBAS")
TEST_RESULTS = os.path.join(RINGCC_DIR, "test_results")

DEFAULT_TOLERANCE = 1e-6  # a.u.

# Ground-state-only ring-CC pipeline -- deliberately stops at xringcc.
# xringcl (Lambda equations) is only needed for EOM/properties, not the
# correlation energy itself, and running it unconditionally here would
# actually corrupt TOTENERG: found 2026-08-13 that xringcl re-solves T2
# from scratch starting from xringcc's own (correct) converged amplitudes,
# then diverges away from them to a different, wrong answer -- a second,
# separate ringcl-side amplitude-convergence bug, not yet diagnosed/fixed
# (see project memory). xvee is also skipped -- no EXCITE requested by any
# case here, and it errors out with no excited states specified.
PIPELINE = ["xjoda", "xvmol", "xvmol2ja", "xvprops", "xvscf", "xvtran",
            "xintprc", "xringcc"]

_EXTRA_LD_PATH = "/apps/compilers/intel/2025/1.0.666/compiler/2025.1/lib"


def load_test_results():
    refs = {}
    with open(TEST_RESULTS) as f:
        lines = [l.rstrip("\n") for l in f]
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        name = lines[i].split()[0]
        record = lines[i + 1].split()[0]
        value = float(lines[i + 2].replace("D", "E").replace("d", "e"))
        refs[name] = (record, value)
        i += 3
    return refs


def _compare_record(scratch, env, name, record, ref_energy, elapsed):
    jr = subprocess.run(
        ["xa2proc", "jareq", "d", record, "1"],
        cwd=scratch, env=env, capture_output=True, text=True, timeout=60,
    )
    m = re.search(r"[-+]?\d*\.\d+[eEdD][-+]?\d+", jr.stdout)
    if not m:
        return name, "ERROR", elapsed, (
            f"xa2proc jareq produced no parseable value for {record}\n"
            f"stdout: {jr.stdout}\nstderr: {jr.stderr}"
        )
    got = float(m.group(0).replace("D", "E").replace("d", "e"))
    diff = abs(got - ref_energy)
    status = "PASS" if diff <= DEFAULT_TOLERANCE else "FAIL"
    detail = f"got {got:.10f}  ref {ref_energy:.10f}  diff {diff:.2e} (tol {DEFAULT_TOLERANCE:.1e})"
    return name, status, elapsed, detail


def run_case(name, record, ref_energy, env):
    case_dir = os.path.join(CASES_DIR, name)
    zmat = os.path.join(case_dir, "ZMAT")
    if not os.path.isfile(zmat):
        return name, "ERROR", 0.0, f"no ZMAT found in {case_dir}"

    with tempfile.TemporaryDirectory(prefix=f"acesii_ringcc_{name}_") as scratch:
        shutil.copy(GENBAS, os.path.join(scratch, "GENBAS"))
        shutil.copy(zmat, os.path.join(scratch, "ZMAT"))

        t0 = time.time()
        log = []
        for exe in PIPELINE:
            try:
                run = subprocess.run(
                    [exe], cwd=scratch, env=env,
                    capture_output=True, text=True, timeout=300,
                )
            except subprocess.TimeoutExpired:
                return name, "ERROR", time.time() - t0, f"{exe} timed out (300s)"
            log.append(run.stdout)
            if run.returncode != 0:
                elapsed = time.time() - t0
                tail = "\n".join(run.stdout.splitlines()[-15:])
                return name, "ERROR", elapsed, f"{exe} exited {run.returncode}\n{tail}"
        elapsed = time.time() - t0

        return _compare_record(scratch, env, name, record, ref_energy, elapsed)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=4,
                         help="cases to run concurrently (default 4)")
    parser.add_argument("--case", action="append", dest="cases",
                         help="run only these case names (repeatable); "
                              "default is all cases in tests/ringcc/cases/")
    args = parser.parse_args()

    if shutil.which("xringcc") is None:
        sys.exit("xringcc not found on PATH -- module load intel/2025.1.0 "
                  "openmpi/5.0.7 first, and make sure ~/Develop/bin is on PATH")

    env = os.environ.copy()
    ld = env.get("LD_LIBRARY_PATH", "")
    if _EXTRA_LD_PATH not in ld:
        env["LD_LIBRARY_PATH"] = _EXTRA_LD_PATH + (":" + ld if ld else "")

    refs = load_test_results()
    all_cases = sorted(
        d for d in os.listdir(CASES_DIR)
        if os.path.isfile(os.path.join(CASES_DIR, d, "ZMAT"))
    )
    cases = args.cases if args.cases else all_cases
    unknown = set(cases) - set(all_cases)
    if unknown:
        sys.exit(f"unknown case(s): {', '.join(sorted(unknown))}")
    missing_ref = [c for c in cases if c not in refs]
    if missing_ref:
        sys.exit(f"case(s) with no reference in test_results: {', '.join(missing_ref)}")

    print(f"Running {len(cases)} case(s), --jobs {args.jobs}\n")
    print(f"{'CASE':<28} {'TIME':>8}  {'STATUS':<6}  DETAIL")
    print("-" * 80)

    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futures = {
            ex.submit(run_case, name, refs[name][0], refs[name][1], env): name
            for name in cases
        }
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
        print(f"{name:<28} {elapsed:>6.2f}s  {status:<6}  {detail_first_line}")
        if status in ("FAIL", "ERROR") and "\n" in detail:
            for line in detail.splitlines()[1:]:
                print(f"    {line}")

    print("-" * 80)
    print(f"{len(cases)} case(s), {n_pass} passed, {n_fail} failed, "
          f"{n_error} errored, {total_time:.1f}s total CPU-time-of-runs")

    sys.exit(0 if (n_fail == 0 and n_error == 0) else 1)


if __name__ == "__main__":
    main()
