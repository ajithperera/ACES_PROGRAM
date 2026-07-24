#!/usr/bin/env python3
"""ACES II fast regression suite.

Runs each small test case in tests/fast/cases/<name>/ against the current
build, pulls the resulting value straight out of JOBARC via
`xa2proc jareq d <RECORD> 1` (a quiet, single-line JOBARC-record print --
no Fortran comparison module needed), compares it to test_results, and
reports PASS/FAIL/ERROR with a tolerance-based gate.

test_results normally names TOTENERG as the record to check. A few cases
(VIB=EXACT analytic-Hessian jobs at a correlated level) never get TOTENERG
written to JOBARC at all -- ACES II's formt2.F deliberately skips that
PUTREC for "SCF Hessian Polyrate" runs (see project memory) -- so their
test_results entry names SCFENEG instead, which is always written and still
verifies the underlying electronic structure.

Usage:
    module load intel/2025.1.0 openmpi/5.0.7   # must be loaded first
    python3 tests/fast/run_tests.py [--jobs N] [--case NAME ...]
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
CASES_DIR = os.path.join(FAST_DIR, "cases")
GENBAS = os.path.join(FAST_DIR, "GENBAS")
TEST_RESULTS = os.path.join(FAST_DIR, "test_results")

DEFAULT_TOLERANCE = 1e-6  # a.u.
TOLERANCE_OVERRIDES = {}

# ifx's runtime libs aren't always on LD_LIBRARY_PATH after `module load` --
# append the known location defensively so callers don't have to remember.
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
        # lines[i+1] is "<RECORD> 1", lines[i+2] is the value
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
    tol = TOLERANCE_OVERRIDES.get(name, DEFAULT_TOLERANCE)
    diff = abs(got - ref_energy)
    status = "PASS" if diff <= tol else "FAIL"
    detail = f"got {got:.10f}  ref {ref_energy:.10f}  diff {diff:.2e} (tol {tol:.1e})"
    return name, status, elapsed, detail


def run_case_zmat(name, zmat, record, ref_energy, env):
    with tempfile.TemporaryDirectory(prefix=f"acesii_fast_{name}_") as scratch:
        shutil.copy(GENBAS, os.path.join(scratch, "GENBAS"))
        shutil.copy(zmat, os.path.join(scratch, "ZMAT"))

        t0 = time.time()
        try:
            run = subprocess.run(
                ["xaces2"], cwd=scratch, env=env,
                capture_output=True, text=True, timeout=300,
            )
        except subprocess.TimeoutExpired:
            return name, "ERROR", time.time() - t0, "xaces2 timed out (300s)"
        elapsed = time.time() - t0

        if run.returncode != 0:
            tail = "\n".join(run.stdout.splitlines()[-15:])
            return name, "ERROR", elapsed, f"xaces2 exited {run.returncode}\n{tail}"

        return _compare_record(scratch, env, name, record, ref_energy, elapsed)


def run_case_script(name, script, record, ref_energy, env):
    # SCRIPT cases chain multiple ACES2 jobs in one shell script, saving an
    # intermediate file (FCMINT/OLDMOS/etc.) outside the working dir before
    # an internal `rm -rf *` cleanup between jobs -- USERDIR must therefore
    # be a directory the script's own cleanup can't reach, i.e. NOT a
    # subdirectory of the script's cwd.
    with tempfile.TemporaryDirectory(prefix=f"acesii_fast_{name}_") as scratch, \
         tempfile.TemporaryDirectory(prefix=f"acesii_fast_{name}_userdir_") as userdir:
        script_env = dict(env)
        script_env["GENBAS"] = GENBAS
        script_env["USERDIR"] = userdir

        t0 = time.time()
        try:
            run = subprocess.run(
                ["sh", script], cwd=scratch, env=script_env,
                capture_output=True, text=True, timeout=300,
            )
        except subprocess.TimeoutExpired:
            return name, "ERROR", time.time() - t0, "script timed out (300s)"
        elapsed = time.time() - t0

        if run.returncode != 0:
            tail = "\n".join(run.stdout.splitlines()[-15:])
            return name, "ERROR", elapsed, f"script exited {run.returncode}\n{tail}"

        return _compare_record(scratch, env, name, record, ref_energy, elapsed)


def run_case(name, record, ref_energy, env):
    case_dir = os.path.join(CASES_DIR, name)
    zmat = os.path.join(case_dir, "ZMAT")
    script = os.path.join(case_dir, "SCRIPT")
    if os.path.isfile(zmat):
        return run_case_zmat(name, zmat, record, ref_energy, env)
    if os.path.isfile(script):
        return run_case_script(name, script, record, ref_energy, env)
    return name, "ERROR", 0.0, f"no ZMAT or SCRIPT found in {case_dir}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=4,
                         help="cases to run concurrently (default 4 -- "
                              "ACES II is serial per-job, so this is a "
                              "straight core count, unlike ACES III's --np)")
    parser.add_argument("--case", action="append", dest="cases",
                         help="run only these case names (repeatable); "
                              "default is all cases in tests/fast/cases/")
    args = parser.parse_args()

    if shutil.which("xaces2") is None:
        sys.exit("xaces2 not found on PATH -- module load intel/2025.1.0 "
                  "openmpi/5.0.7 first, and make sure ~/Develop/bin is on PATH")

    env = os.environ.copy()
    ld = env.get("LD_LIBRARY_PATH", "")
    if _EXTRA_LD_PATH not in ld:
        env["LD_LIBRARY_PATH"] = _EXTRA_LD_PATH + (":" + ld if ld else "")

    refs = load_test_results()
    all_cases = sorted(
        d for d in os.listdir(CASES_DIR)
        if os.path.isfile(os.path.join(CASES_DIR, d, "ZMAT"))
        or os.path.isfile(os.path.join(CASES_DIR, d, "SCRIPT"))
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
