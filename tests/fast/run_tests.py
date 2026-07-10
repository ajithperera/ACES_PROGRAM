#!/usr/bin/env python3
"""ACES III fast regression suite.

Runs each small test case in tests/fast/cases/<name>/ against the current
build, compares the resulting JOBARC energy to test_results_fast via the
existing xtest_compare tool, and reports PASS/FAIL with a tolerance-based
gate (not just a printed RMS error for a human to eyeball).

Usage:
    module load intel/2025.1.0 openmpi/5.0.7   # must be loaded first
    python3 tests/fast/run_tests.py [--jobs N] [--np N] [--case NAME ...]
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

FAST_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(os.path.dirname(FAST_DIR))
CASES_DIR = os.path.join(FAST_DIR, "cases")
GENBAS = os.path.join(FAST_DIR, "GENBAS")
XACES3 = os.path.join(REPO_DIR, "bin", "xaces3")
XTEST_COMPARE = os.path.join(REPO_DIR, "bin", "xtest_compare")

DEFAULT_TOLERANCE = 1e-6  # a.u., RMS error gate

# Per-case tolerance overrides, for methods that may not reproduce to the
# same tightness every run (e.g. loosely-converging EOM/STEOM iterative
# roots). Empty for now -- every case reproduced to ~1e-10 or better when
# this suite was built, so the flat default is used unless proven otherwise.
TOLERANCE_OVERRIDES = {}

RMS_ERROR_RE = re.compile(r"RMS error\s+([-\d.eEdD+]+)")


def run_case(name, np):
    case_dir = os.path.join(CASES_DIR, name)
    zmat = os.path.join(case_dir, "ZMAT")
    if not os.path.isfile(zmat):
        return name, False, 0.0, f"no ZMAT found at {zmat}"

    with tempfile.TemporaryDirectory(prefix=f"acesiii_fast_{name}_") as scratch:
        shutil.copy(GENBAS, os.path.join(scratch, "GENBAS"))
        shutil.copy(zmat, os.path.join(scratch, "ZMAT"))

        env = os.environ.copy()
        env["ACES_EXE_PATH"] = os.path.join(REPO_DIR, "bin")

        t0 = time.time()
        run = subprocess.run(
            ["mpirun", "--oversubscribe", "-np", str(np), XACES3],
            cwd=scratch, env=env, capture_output=True, text=True, timeout=600,
        )
        elapsed = time.time() - t0

        if run.returncode != 0:
            tail = "\n".join(run.stdout.splitlines()[-15:])
            return name, False, elapsed, f"xaces3 exited {run.returncode}\n{tail}"

        env["ACES_EXE_PATH"] = FAST_DIR
        cmp_run = subprocess.run(
            ["mpirun", "--oversubscribe", "-np", "1", XTEST_COMPARE, name],
            cwd=scratch, env=env, capture_output=True, text=True, timeout=60,
        )

        m = RMS_ERROR_RE.search(cmp_run.stdout)
        if not m:
            return name, False, elapsed, (
                "xtest_compare produced no RMS error line\n"
                f"stdout: {cmp_run.stdout}\nstderr: {cmp_run.stderr}"
            )

        rms_error = float(m.group(1).replace("D", "E").replace("d", "e"))
        tol = TOLERANCE_OVERRIDES.get(name, DEFAULT_TOLERANCE)
        passed = rms_error <= tol
        detail = f"RMS error {rms_error:.3e} (tol {tol:.1e})"
        return name, passed, elapsed, detail


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=1,
                         help="cases to run concurrently (default 1 -- each "
                              "case already uses --np MPI ranks, so keep this "
                              "low on a 2-4 core allocation)")
    parser.add_argument("--np", type=int, default=4,
                         help="MPI ranks per case (default 4 -- 2 ranks was "
                              "found unreliable for gradient/EOM cases, "
                              "crashing inside the SIP block-server layer)")
    parser.add_argument("--case", action="append", dest="cases",
                         help="run only these case names (repeatable); "
                              "default is all cases in tests/fast/cases/")
    args = parser.parse_args()

    if not os.path.isfile(XACES3):
        sys.exit(f"bin/xaces3 not found at {XACES3} -- build it first")

    all_cases = sorted(
        d for d in os.listdir(CASES_DIR)
        if os.path.isfile(os.path.join(CASES_DIR, d, "ZMAT"))
    )
    cases = args.cases if args.cases else all_cases
    unknown = set(cases) - set(all_cases)
    if unknown:
        sys.exit(f"unknown case(s): {', '.join(sorted(unknown))}")

    print(f"Running {len(cases)} case(s), --np {args.np}, --jobs {args.jobs}\n")

    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(run_case, name, args.np): name for name in cases}
        for future in concurrent.futures.as_completed(futures):
            results.append(future.result())

    results.sort(key=lambda r: r[0])

    print(f"{'CASE':24s} {'TIME':>8s}  {'STATUS':6s}  DETAIL")
    print("-" * 80)
    n_failed = 0
    for name, passed, elapsed, detail in results:
        status = "PASS" if passed else "FAIL"
        if not passed:
            n_failed += 1
        first_line = detail.splitlines()[0] if detail else ""
        print(f"{name:24s} {elapsed:7.2f}s  {status:6s}  {first_line}")
        if not passed and len(detail.splitlines()) > 1:
            for line in detail.splitlines()[1:]:
                print(f"    {line}")

    total_time = sum(r[2] for r in results)
    print("-" * 80)
    print(f"{len(results)} case(s), {n_failed} failed, "
          f"{total_time:.1f}s total CPU-time-of-runs")

    sys.exit(1 if n_failed else 0)


if __name__ == "__main__":
    main()
