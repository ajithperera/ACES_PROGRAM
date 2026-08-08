#!/usr/bin/env python3
"""Drive pyaces_run_case.py over every SCF/CC/EOM/PCCD case classified as
pyaces-compatible by classify_suite.py. Fresh subprocess per case, fresh
isolated working dir per case (ZMAT copied in, GENBAS symlinked).
"""
import csv
import os
import re
import shutil
import subprocess
import sys
import time

ACES_DIAG_RE = re.compile(r"@[A-Z][A-Z0-9_]*:")

BASE = "/blue/ufhpc/perera/acespy_build"
CLASSIFICATION = os.path.join(BASE, "suite_classification.tsv")
SUITE_CASES = os.path.expanduser("~/Develop/ACESII_build/tests/fast/cases")
GENBAS = os.path.expanduser("~/Develop/ACESII_build/tests/fast/GENBAS")
RUN_DIR = "/blue/ufhpc/perera/pyaces_suite_test"
RUNNER = os.path.join(BASE, "pyaces_run_case.py")
TOLERANCE = 1e-6

CASE_FILTER = set(sys.argv[1:]) or None


def main():
    with open(CLASSIFICATION) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    candidates = [r for r in rows if r["category"] in ("SCF", "CC", "EOM", "PCCD")]
    if CASE_FILTER:
        candidates = [r for r in candidates if r["name"] in CASE_FILTER]

    os.makedirs(RUN_DIR, exist_ok=True)
    results = []
    for i, row in enumerate(candidates, 1):
        name = row["name"]
        entrypoint = row["entrypoint"]
        iuhf = row["iuhf"]
        record = row["record"]
        ref_value = float(row["ref_value"].replace("D", "E").replace("d", "e"))

        case_dir = os.path.join(RUN_DIR, name)
        shutil.rmtree(case_dir, ignore_errors=True)
        os.makedirs(case_dir)
        shutil.copy(os.path.join(SUITE_CASES, name, "ZMAT"), os.path.join(case_dir, "ZMAT"))
        os.symlink(GENBAS, os.path.join(case_dir, "GENBAS"))

        t0 = time.time()
        try:
            proc = subprocess.run(
                ["python3", "-u", RUNNER, case_dir, entrypoint, iuhf, record],
                cwd=case_dir, capture_output=True, text=True, timeout=300,
            )
            timed_out = False
        except subprocess.TimeoutExpired as e:
            proc = None
            timed_out = True
        elapsed = time.time() - t0

        log_path = os.path.join(case_dir, "run.log")
        with open(log_path, "w") as f:
            if timed_out:
                f.write("TIMEOUT after 300s\n")
            else:
                f.write("=== stdout ===\n" + proc.stdout + "\n=== stderr ===\n" + proc.stderr)

        # CRASH/FAIL bucketing: subprocess.run reports a NEGATIVE returncode
        # only when the child was killed by a signal (real SIGSEGV/SIGABRT/etc
        # -- an actual memory-corruption-class crash). A POSITIVE nonzero
        # returncode means the process terminated itself deliberately -- most
        # commonly a Fortran-side aces_exit()/STOP (clean assertion failure,
        # non-convergence) or an ifort runtime "severe" I/O error (abort()
        # inside forrtl, e.g. missing scratch file) -- neither is a crash bug
        # in the collision-hunting sense, both were previously indistinguishable
        # from a real SIGSEGV in this harness's "CRASH" bucket.
        # Named markers for cases where the "@ROUTINENAME:" banner itself
        # didn't make it out (aces_exit's own C exit() sometimes drops its
        # last buffered Fortran write, e.g. case 117a: halts right after
        # "@LOCATE: Label ... not found", the aces_exit banner that should
        # follow never appears in the captured log at all).
        CLEAN_FAIL_MARKERS = (
            "Assertion failed", "did not converge", "@VCC:",
            "forrtl: severe", "System error:", "PYACES_CASE_RESULT:",
        )
        status = "ERROR"
        detail = ""
        if timed_out:
            status = "TIMEOUT"
        elif proc.returncode != 0:
            combined_log = proc.stdout + proc.stderr
            if proc.returncode < 0:
                status = "CRASH"
                detail = f"killed by signal {-proc.returncode}"
            elif any(marker in combined_log for marker in CLEAN_FAIL_MARKERS) or ACES_DIAG_RE.search(combined_log):
                status = "FAIL_CLEAN"
                detail = f"exit_code={proc.returncode}"
            else:
                status = "CRASH"
                detail = f"exit_code={proc.returncode} (unrecognized, no clean-exit marker found)"
            for line in proc.stdout.splitlines():
                if "PYACES_CASE_RESULT" in line:
                    detail += " | " + line
                    break
        else:
            got = None
            for line in proc.stdout.splitlines():
                if line.startswith("PYACES_CASE_RESULT: SUCCESS"):
                    got = float(line.split()[-1])
            if got is None:
                status = "NO_RESULT"
            else:
                diff = abs(got - ref_value)
                if diff <= TOLERANCE:
                    status = "PASS"
                else:
                    status = "MISMATCH"
                detail = f"got={got:.10f} ref={ref_value:.10f} diff={diff:.2e}"

        results.append((name, row["category"], entrypoint, status, detail, elapsed))
        print(f"[{i}/{len(candidates)}] {name:20s} {entrypoint:10s} {status:10s} {detail} ({elapsed:.1f}s)")
        sys.stdout.flush()

    print()
    counts = {}
    for r in results:
        counts[r[3]] = counts.get(r[3], 0) + 1
    print("=== Summary ===")
    for k, v in sorted(counts.items(), key=lambda x: -x[1]):
        print(f"  {k}: {v}")

    summary_path = os.path.join(RUN_DIR, "summary.tsv")
    with open(summary_path, "w") as f:
        f.write("name\tcategory\tentrypoint\tstatus\tdetail\telapsed\n")
        for r in results:
            f.write("\t".join(str(x) for x in r) + "\n")
    print(f"\nWrote {summary_path}")


if __name__ == "__main__":
    main()
