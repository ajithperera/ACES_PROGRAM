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
import functools
import os
import resource
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
DEFAULT_TIMEOUT = 600  # seconds, per-case xaces3 wall-clock budget

# The budgets below are sized for a HiPerGator compute node. A slower or
# busier machine (a workstation running several cases at once, say) needs the
# same runs to be allowed more wall clock before they count as a failure --
# a timeout says nothing about correctness. Scale them all with --timeout-scale
# rather than editing the table per machine.
TIMEOUT_SCALE = float(os.environ.get("ACES3_TIMEOUT_SCALE", "1"))

# Per-case timeout overrides. The gten_* cases are full CCSD g-tensor
# property runs (SCF+transform+CCSD+lambda+response, np=4) that measured
# ~454s/~358s wall clock end-to-end on a dedicated reference run -- close
# enough to the 600s default that this suite's own (--oversubscribe,
# tempfile-scratch-on-/tmp) conditions blow past it. gten_no confirmed this
# directly: it finished at 599.20s under the old default, a hair from a
# spurious timeout-FAIL. Give real headroom rather than nudging the number
# until it barely passes.
TIMEOUT_OVERRIDES = {
    "gten_h2op": 1200,
    "gten_cn": 1200,
    "gten_no": 1200,
    "ccsd_f12": 1200,
    "ccsdt_f12": 1200,
    "dea_eomcc_f12": 2700,
}

# Per-case tolerance overrides, for methods that may not reproduce to the
# same tightness every run (e.g. loosely-converging EOM/STEOM iterative
# roots). Empty for now -- every case reproduced to ~1e-10 or better when
# this suite was built, so the flat default is used unless proven otherwise.
TOLERANCE_OVERRIDES = {}

RMS_ERROR_RE = re.compile(r"RMS error\s+([-\d.eEdD+]+)")

# g-tensor table verification (gten_* cases). SCFENEG alone doesn't check
# the actual g-tensor property calculation these cases exist to exercise --
# see reference_gtensor.txt in each such case dir, extracted from a known-
# good run's summary.out (unit 66, written by
# src/sia/sip/aces_instructions/print_rel_info.F). A case dir carrying a
# reference_gtensor.txt gets its own summary.out parsed and diffed
# elementwise against it, on top of the normal xtest_compare check.
GTENSOR_TOLERANCE = 1e-3  # ppm, RMS over all tensor elements across tables
GTENSOR_TOLERANCE_OVERRIDES = {}

GTENSOR_HEADER_RE = re.compile(r"g-tensor:\s*(.+?)\s*\(ppm\)")
GTENSOR_COLHDR_RE = re.compile(r"^\s*X\s+Y\s+Z\s*$")
GTENSOR_ROW_RE = re.compile(
    r"^\s*([XYZ])\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s*$")

# General property-block verification (dipole/quadrupole/polar/... cases).
# TOTENERG alone only proves the CCSD ground state reproduced -- it says
# nothing about the property calculation these cases exist to exercise (the
# 2026-07-31 eom_ee_so case is the cautionary example: TOTENERG and the
# singlet EOM roots matched the archive exactly while the triplet roots
# silently failed to converge to all-zero, which only showed up by actually
# diffing the printed root table). A case dir carrying a
# reference_properties.txt (a verbatim snippet copied from a known-good
# summary.out, header line first) gets that exact number of lines located in
# the actual summary.out by matching the header line, then every number in
# both blocks is extracted in order and compared elementwise.
PROPERTY_TOLERANCE = 1e-4  # a.u., RMS over all numbers in the reference block
PROPERTY_TOLERANCE_OVERRIDES = {}

NUMBER_RE = re.compile(r"[-+]?\d+\.\d+(?:[DdEe][-+]?\d+)?")


def parse_numbers(text):
    return [float(tok.replace("D", "E").replace("d", "e"))
            for tok in NUMBER_RE.findall(text)]


def compare_property_block(name, ref_text, actual_text):
    ref_lines = ref_text.splitlines()
    if not ref_lines:
        return False, 0.0, "reference_properties.txt is empty"
    header = ref_lines[0].strip()
    actual_lines = actual_text.splitlines()
    start = next((i for i, l in enumerate(actual_lines)
                  if l.strip() == header), None)
    if start is None:
        return False, 0.0, (
            f"reference_properties.txt header line not found in actual "
            f"summary.out: {header!r}"
        )
    actual_block = "\n".join(actual_lines[start:start + len(ref_lines)])
    ref_nums = parse_numbers(ref_text)
    actual_nums = parse_numbers(actual_block)
    if len(ref_nums) != len(actual_nums) or not ref_nums:
        return False, 0.0, (
            f"property block number count mismatch: expected {len(ref_nums)}, "
            f"got {len(actual_nums)}"
        )
    diffs = [a - r for a, r in zip(actual_nums, ref_nums)]
    rms = (sum(d * d for d in diffs) / len(diffs)) ** 0.5
    tol = PROPERTY_TOLERANCE_OVERRIDES.get(name, PROPERTY_TOLERANCE)
    return rms <= tol, rms, f"property RMS error {rms:.3e} (tol {tol:.1e})"


def parse_gtensor_tables(text):
    """Extract {table name: 3x3 matrix (rows/cols in X,Y,Z order)} from a
    summary.out-shaped block of text. Table headers and their column-header/
    data rows can have unrelated lines (timing prints, blank lines) between
    them, so this scans forward from each header looking for the next X/Y/Z
    column-header line, then the next 3 row lines after that."""
    tables = {}
    lines = text.splitlines()
    i = 0
    pending_name = None
    while i < len(lines):
        m = GTENSOR_HEADER_RE.search(lines[i])
        if m:
            pending_name = re.sub(r"\s+", " ", m.group(1)).strip()
            i += 1
            continue
        if pending_name and GTENSOR_COLHDR_RE.match(lines[i]):
            rows = {}
            j = i + 1
            while j < len(lines) and len(rows) < 3:
                rm = GTENSOR_ROW_RE.match(lines[j])
                if rm:
                    rows[rm.group(1)] = [float(rm.group(k)) for k in (2, 3, 4)]
                j += 1
            if len(rows) == 3:
                tables[pending_name] = [rows["X"], rows["Y"], rows["Z"]]
            pending_name = None
            i = j
            continue
        i += 1
    return tables


def compare_gtensor(name, ref_text, actual_text):
    ref_tables = parse_gtensor_tables(ref_text)
    actual_tables = parse_gtensor_tables(actual_text)
    missing = [k for k in ref_tables if k not in actual_tables]
    if missing or not ref_tables:
        return False, 0.0, (
            f"g-tensor table(s) missing from actual summary.out: "
            f"{', '.join(missing) if missing else '(no reference tables parsed)'}"
        )
    diffs = []
    for tname, ref_mat in ref_tables.items():
        act_mat = actual_tables[tname]
        for r in range(3):
            for c in range(3):
                diffs.append(act_mat[r][c] - ref_mat[r][c])
    rms = (sum(d * d for d in diffs) / len(diffs)) ** 0.5
    tol = GTENSOR_TOLERANCE_OVERRIDES.get(name, GTENSOR_TOLERANCE)
    return rms <= tol, rms, f"g-tensor RMS error {rms:.3e} (tol {tol:.1e})"


def raise_stack_limit():
    """Give the ranks the stack the SOC superinstructions need.

    steom_so_diag.F (and the other aces_instructions_extnd routines it is
    modelled on) declares its working matrices as automatic arrays --
    AAA/BBB/CMATR/VL/VR are all msize x msize, and msize is 4 * nocc * nvirt --
    so case steom_so alone wants a few hundred MB of *stack* per rank. Under a
    login shell's default 8MB that is an immediate SIGSEGV inside
    steom_so_diag at line 32, which looks like a code bug and is not one; a
    batch environment like HiPerGator's simply has no such limit, which is why
    the baseline never saw it. Children of this process inherit the limit, so
    raising it here covers every mpirun the suite launches.
    """
    try:
        soft, hard = resource.getrlimit(resource.RLIMIT_STACK)
        if soft != hard:
            resource.setrlimit(resource.RLIMIT_STACK, (hard, hard))
    except (ValueError, OSError):
        pass  # not fatal: only the SOC cases care


@functools.lru_cache(maxsize=1)
def mpirun_extra_args():
    """--oversubscribe is an OpenMPI flag: it is required there whenever the
    suite runs more ranks than the box has cores, but Intel MPI's Hydra
    launcher rejects it outright (it oversubscribes by default), so hardcoding
    it makes every case fail with an argument error on an Intel MPI build.
    Ask mpirun which one it is."""
    try:
        out = subprocess.run(["mpirun", "--version"], capture_output=True,
                             text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return []
    banner = (out.stdout or "") + (out.stderr or "")
    if "Intel(R) MPI" in banner or "Intel MPI" in banner:
        return []
    return ["--oversubscribe"]



# xaces3 says why it died in a diagnostic line on stdout, then calls
# dump_block_ids()/array_block_summary(), which print hundreds of lines of
# block tables before mpi_abort finally exits. Keeping only the stdout tail
# therefore keeps only the dump and throws the reason away -- that is how the
# gten_* "will not fit on any stack" abort stayed invisible for a month --
# and stderr, which carries the Fortran traceback, was dropped entirely.
# Pull out the diagnostic lines by name, then the tail, then stderr.
ABORT_LINE_RE = re.compile(
    r"error|cannot|will not fit|abort|assertion|insufficient|exceeds",
    re.IGNORECASE)


def summarize_failure(stdout, stderr, ncontext=2, nmatch=12, ntail=8):
    out = []
    lines = (stdout or "").splitlines()
    hits = [i for i, l in enumerate(lines) if ABORT_LINE_RE.search(l)]
    if hits:
        out.append("-- diagnostics --")
        shown = set()
        for i in hits[:nmatch]:
            for j in range(i, min(i + 1 + ncontext, len(lines))):
                if j not in shown:
                    shown.add(j)
                    out.append(lines[j])
    if lines:
        out.append("-- last %d stdout lines --" % ntail)
        out.extend(lines[-ntail:])
    err = (stderr or "").splitlines()
    if err:
        out.append("-- last %d stderr lines --" % ntail)
        out.extend(err[-ntail:])
    return "\n".join(out)

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

        timeout = int(TIMEOUT_OVERRIDES.get(name, DEFAULT_TIMEOUT) * TIMEOUT_SCALE)
        t0 = time.time()
        try:
            run = subprocess.run(
                ["mpirun", *mpirun_extra_args(), "-np", str(np), XACES3],
                cwd=scratch, env=env, capture_output=True, text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired:
            elapsed = time.time() - t0
            return name, False, elapsed, f"xaces3 timed out after {timeout}s"
        elapsed = time.time() - t0

        if run.returncode != 0:
            return name, False, elapsed, (
                f"xaces3 exited {run.returncode}\n"
                + summarize_failure(run.stdout, run.stderr)
            )

        env["ACES_EXE_PATH"] = FAST_DIR
        cmp_run = subprocess.run(
            ["mpirun", *mpirun_extra_args(), "-np", "1", XTEST_COMPARE, name],
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

        gtensor_ref = os.path.join(case_dir, "reference_gtensor.txt")
        if os.path.isfile(gtensor_ref):
            summary_path = os.path.join(scratch, "summary.out")
            if not os.path.isfile(summary_path):
                passed = False
                detail += "\nreference_gtensor.txt present but no summary.out produced"
            else:
                with open(gtensor_ref) as f:
                    ref_text = f.read()
                with open(summary_path) as f:
                    actual_text = f.read()
                g_ok, g_rms, g_detail = compare_gtensor(name, ref_text, actual_text)
                passed = passed and g_ok
                detail += "\n" + g_detail

        property_ref = os.path.join(case_dir, "reference_properties.txt")
        if os.path.isfile(property_ref):
            summary_path = os.path.join(scratch, "summary.out")
            if not os.path.isfile(summary_path):
                passed = False
                detail += "\nreference_properties.txt present but no summary.out produced"
            else:
                with open(property_ref) as f:
                    ref_text = f.read()
                with open(summary_path) as f:
                    actual_text = f.read()
                p_ok, p_rms, p_detail = compare_property_block(name, ref_text, actual_text)
                passed = passed and p_ok
                detail += "\n" + p_detail

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
    parser.add_argument("--timeout-scale", type=float, default=None,
                        help="multiply every per-case timeout by this factor "
                             "(also settable as ACES3_TIMEOUT_SCALE); use >1 on "
                             "a machine slower than a HiPerGator node")
    parser.add_argument("--case", action="append", dest="cases",
                         help="run only these case names (repeatable); "
                              "default is all cases in tests/fast/cases/")
    args = parser.parse_args()

    raise_stack_limit()

    if args.timeout_scale is not None:
        global TIMEOUT_SCALE
        TIMEOUT_SCALE = args.timeout_scale

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
