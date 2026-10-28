#!/usr/bin/env python3
"""ACES II ring-CC Lambda (xringcl) regression suite.

Covers CALC=DRCCD Lambda-equation solving (xringcl), for both REF=RHF and
REF=UHF. Runs the full manual pipeline (same as tests/ringcc/run_tests.py)
but continues through xringcl instead of stopping at xringcc.

Reference strategy: NOT matched against the 2018 archive
(/blue/bartlett/perera/Testing/Testdir_rcc/drccd_rhf_lambda.out). That
archive turned out to predate a deliberate refactor (the RHF/UHF-identity
factor of 2 in ringcl/drcl_dwmbej.F was moved into
ringcl/rcl_lrngdrv_r.F's RHF-collapsed contraction sometime after 2018 --
see the routine's own comment and project memory, 2026-08-13) and is not a
reliable reference for the current, fully-migrated source.

Instead, for DRCCD specifically, REF=RHF and REF=UHF are mathematically
required to agree (unlike RCCD-with-exchange, where the archive itself
shows RHF and UHF legitimately compute different quantities for this
truncated method -- see tests/ringcc's rccd_rhf/rccd_uhf notes). So each
case here is a matched RHF/UHF pair, and the test asserts the two
TOTENERG values agree to DEFAULT_TOLERANCE, rather than checking either
one against a fixed number. This was verified by hand 2026-08-13 for both
the STO-3G and DZ systems below (closed-shell UHF confirmed spin-pure via
identical AAAA/BBBB checksums, no symmetry breaking) before being encoded
here as a regression check.

Usage:
    module load intel/2025.1.0 openmpi/5.0.7   # must be loaded first
    python3 tests/ringcl/run_tests.py [--jobs N] [--pair NAME ...]
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

RINGCL_DIR = os.path.dirname(os.path.abspath(__file__))
CASES_DIR = os.path.join(RINGCL_DIR, "cases")
# One GENBAS for all suites (ACES_PROGRAM/basis).
GENBAS = os.path.join(RINGCL_DIR, os.pardir, os.pardir, os.pardir, "basis", "GENBAS")

DEFAULT_TOLERANCE = 1e-6  # a.u.

PIPELINE = ["xjoda", "xvmol", "xvmol2ja", "xvprops", "xvscf", "xvtran",
            "xintprc", "xringcc", "xringcl"]

_EXTRA_LD_PATH = "/apps/compilers/intel/2025/1.0.666/compiler/2025.1/lib"

# (pair name, rhf case dir, uhf case dir)
PAIRS = [
    ("drccd_sto3g", "drccd_sto3g_rhf", "drccd_sto3g_uhf"),
    ("drccd_dz", "drccd_dz_rhf", "drccd_dz_uhf"),
]


def _run_case(case_name, env):
    case_dir = os.path.join(CASES_DIR, case_name)
    zmat = os.path.join(case_dir, "ZMAT")
    if not os.path.isfile(zmat):
        return None, f"no ZMAT found in {case_dir}"

    with tempfile.TemporaryDirectory(prefix=f"acesii_ringcl_{case_name}_") as scratch:
        shutil.copy(GENBAS, os.path.join(scratch, "GENBAS"))
        shutil.copy(zmat, os.path.join(scratch, "ZMAT"))

        for exe in PIPELINE:
            try:
                run = subprocess.run(
                    [exe], cwd=scratch, env=env,
                    capture_output=True, text=True, timeout=300,
                )
            except subprocess.TimeoutExpired:
                return None, f"{exe} timed out (300s)"
            if run.returncode != 0:
                tail = "\n".join(run.stdout.splitlines()[-15:])
                return None, f"{exe} exited {run.returncode}\n{tail}"

        jr = subprocess.run(
            ["xa2proc", "jareq", "d", "TOTENERG", "1"],
            cwd=scratch, env=env, capture_output=True, text=True, timeout=60,
        )
        m = re.search(r"[-+]?\d*\.\d+[eEdD][-+]?\d+", jr.stdout)
        if not m:
            return None, (
                f"xa2proc jareq produced no parseable TOTENERG\n"
                f"stdout: {jr.stdout}\nstderr: {jr.stderr}"
            )
        return float(m.group(0).replace("D", "E").replace("d", "e")), None


def run_pair(pair_name, rhf_case, uhf_case, env):
    t0 = time.time()
    rhf_val, rhf_err = _run_case(rhf_case, env)
    uhf_val, uhf_err = _run_case(uhf_case, env)
    elapsed = time.time() - t0

    if rhf_err or uhf_err:
        detail = ""
        if rhf_err:
            detail += f"RHF ({rhf_case}): {rhf_err}\n"
        if uhf_err:
            detail += f"UHF ({uhf_case}): {uhf_err}\n"
        return pair_name, "ERROR", elapsed, detail.strip()

    diff = abs(rhf_val - uhf_val)
    status = "PASS" if diff <= DEFAULT_TOLERANCE else "FAIL"
    detail = (f"RHF {rhf_val:.10f}  UHF {uhf_val:.10f}  "
              f"diff {diff:.2e} (tol {DEFAULT_TOLERANCE:.1e})")
    return pair_name, status, elapsed, detail


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=4,
                         help="pairs to run concurrently (default 4)")
    parser.add_argument("--pair", action="append", dest="pairs",
                         help="run only these pair names (repeatable); "
                              "default is all pairs")
    args = parser.parse_args()

    if shutil.which("xringcl") is None:
        sys.exit("xringcl not found on PATH -- module load intel/2025.1.0 "
                  "openmpi/5.0.7 first, and make sure ~/Develop/bin is on PATH")

    env = os.environ.copy()
    ld = env.get("LD_LIBRARY_PATH", "")
    if _EXTRA_LD_PATH not in ld:
        env["LD_LIBRARY_PATH"] = _EXTRA_LD_PATH + (":" + ld if ld else "")

    all_pair_names = [p[0] for p in PAIRS]
    pair_names = args.pairs if args.pairs else all_pair_names
    unknown = set(pair_names) - set(all_pair_names)
    if unknown:
        sys.exit(f"unknown pair(s): {', '.join(sorted(unknown))}")
    pairs_by_name = {p[0]: p for p in PAIRS}

    print(f"Running {len(pair_names)} pair(s), --jobs {args.jobs}\n")
    print(f"{'PAIR':<20} {'TIME':>8}  {'STATUS':<6}  DETAIL")
    print("-" * 90)

    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futures = {
            ex.submit(run_pair, name, pairs_by_name[name][1],
                      pairs_by_name[name][2], env): name
            for name in pair_names
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
        print(f"{name:<20} {elapsed:>6.2f}s  {status:<6}  {detail_first_line}")
        if status in ("FAIL", "ERROR") and "\n" in detail:
            for line in detail.splitlines()[1:]:
                print(f"    {line}")

    print("-" * 90)
    print(f"{len(pair_names)} pair(s), {n_pass} passed, {n_fail} failed, "
          f"{n_error} errored, {total_time:.1f}s total CPU-time-of-runs")

    sys.exit(0 if (n_fail == 0 and n_error == 0) else 1)


if __name__ == "__main__":
    main()
