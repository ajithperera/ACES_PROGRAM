#!/usr/bin/env python3
"""OO-pCCD slow regression tier (OO_METHOD=RFO, cc-pVDZ/cc-pVTZ).

Each case in tests/pccd_slow/cases/<name>/ has a ZMAT, optional data files
(oorot.in, c0.dat, pccd_*.in), all copied into the work directory, and an
EXPECT file:

    energy <E>          TOTENERG of the converged minimum (tol 1e-6)
    grads <N>           gradients the reference run needed; more than
                        1.5*N+3 is reported as SLOW (not a failure)
    fail <text>         the run must stop with a nonzero exit and <text>
                        in its output (instead of energy/grads)
    stationary <E> <l>  with fail: the last checked stationary point
                        (energy, lowest Hessian eigenvalue) it reports
    long                run only with --long
    timeout <s>         wall-clock limit (default 7200 s)

A converged case must also print "converged to a minimum" (the RFO
Hessian stability check passed). Jobs run through xaces2 (PROGRAM=PCCD
hands off to xrunpccd). The output of debug builds is huge (hundreds of
GB for the TZ cases), so it is streamed through a filter that drops the
matrix dumps and only run.log is kept, in the work directory of failed
cases (successful ones are deleted unless --keep).

Usage (real compute: submit through sbatch, see submit.sbatch):
    module load intel/2025.1.0 openmpi/5.0.7
    python3 tests/pccd_slow/run_slow.py [--jobs N] [--long]
            [--case NAME ...] [--workdir DIR] [--keep]
"""
import argparse
import concurrent.futures
import importlib.util
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

SLOW_DIR = os.path.dirname(os.path.abspath(__file__))
CASES_DIR = os.path.join(SLOW_DIR, "cases")
FAST_DIR = os.path.join(os.path.dirname(SLOW_DIR), "fast")
GENBAS = os.path.join(SLOW_DIR, os.pardir, os.pardir, os.pardir, "basis", "GENBAS")

# JOBARC record comparison and environment setup shared with the fast suite
# (no __pycache__ left in tests/fast).
sys.dont_write_bytecode = True
_spec = importlib.util.spec_from_file_location(
    "fast_run_tests", os.path.join(FAST_DIR, "run_tests.py"))
fast = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fast)

DEFAULT_TIMEOUT = 7200
STAB_TOL_E = 1e-6
STAB_TOL_L = 1e-4

# Matrix/vector dumps of the debug builds (same filter as the validation
# scripts of 10/2026).
DUMP = re.compile(r"^[ 0-9.EeDd+-]*$|^ +ROW |^ +COLUMN |^\[ *[0-9]")
MINIMUM = re.compile(r"converged to a minimum in +(\d+) gradients")
STATIONARY = re.compile(r"last stationary point checked: E = *(\S+), "
                        r"lowest Hessian eigenvalue *(\S+)")


def load_expect(name):
    exp = {"timeout": DEFAULT_TIMEOUT, "long": False}
    with open(os.path.join(CASES_DIR, name, "EXPECT")) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            key, _, val = line.partition(" ")
            val = val.strip()
            if key == "energy":
                exp["energy"] = float(val)
            elif key == "grads":
                exp["grads"] = int(val)
            elif key == "fail":
                exp["fail"] = val
            elif key == "stationary":
                e, l = val.split()
                exp["stationary"] = (float(e), float(l))
            elif key == "long":
                exp["long"] = True
            elif key == "timeout":
                exp["timeout"] = int(val)
            else:
                raise ValueError(f"{name}/EXPECT: unknown key {key!r}")
    if ("energy" in exp) == ("fail" in exp):
        raise ValueError(f"{name}/EXPECT: needs exactly one of energy, fail")
    return exp


def run_filtered(cwd, env, timeout):
    """Run xaces2 in cwd, keep the filtered output in run.log."""
    log = os.path.join(cwd, "run.log")
    t0 = time.time()
    with open(log, "w", errors="replace") as out:
        proc = subprocess.Popen(["xaces2"], cwd=cwd, env=env,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT)
        timed_out = False
        for raw in proc.stdout:
            line = raw.decode("ascii", errors="replace")
            if not DUMP.match(line.rstrip("\n")):
                out.write(line)
            if time.time() - t0 > timeout:
                proc.kill()
                timed_out = True
                break
        proc.wait()
    return proc.returncode, timed_out, time.time() - t0, log


def run_case(name, exp, env, workroot, keep):
    work = os.path.join(workroot, name)
    shutil.rmtree(work, ignore_errors=True)
    shutil.copytree(os.path.join(CASES_DIR, name), work,
                    ignore=shutil.ignore_patterns("EXPECT"))
    os.symlink(GENBAS, os.path.join(work, "GENBAS"))

    rc, timed_out, elapsed, log = run_filtered(work, env, exp["timeout"])
    with open(log, errors="replace") as f:
        text = f.read()
    tail = "\n".join(text.splitlines()[-12:])

    if timed_out:
        status, detail = "ERROR", f"timed out ({exp['timeout']} s)\n{tail}"
    elif "fail" in exp:
        status, detail = check_fail(exp, rc, text, tail)
    else:
        status, detail = check_minimum(name, exp, rc, text, tail, work,
                                       env, elapsed)

    if status in ("PASS", "SLOW") and not keep:
        shutil.rmtree(work, ignore_errors=True)
    else:
        for junk in ("MOINTS", "MOABCD", "IIII", "GAMLAM", "FOCKCD",
                     "MOL_ERDOED", "IIJJ", "IJIJ", "IJKL"):
            try:
                os.remove(os.path.join(work, junk))
            except OSError:
                pass
        detail += f"\n(work directory kept: {work})"
    return name, status, elapsed, detail


def check_minimum(name, exp, rc, text, tail, work, env, elapsed):
    if rc != 0:
        return "ERROR", f"xaces2 exited {rc}\n{tail}"
    m = MINIMUM.findall(text)
    if not m:
        return "FAIL", f"no verified minimum in the output\n{tail}"
    grads = int(m[-1])
    _, status, _, detail = fast._compare_record(
        work, env, name, "TOTENERG", exp["energy"], elapsed)
    detail += f"  grads {grads}"
    if "grads" in exp:
        detail += f" (ref {exp['grads']})"
        if status == "PASS" and grads > 1.5 * exp["grads"] + 3:
            status = "SLOW"
    return status, detail


def check_fail(exp, rc, text, tail):
    if rc == 0:
        return "FAIL", f"expected a failure, xaces2 exited 0\n{tail}"
    if exp["fail"] not in text:
        return "FAIL", f"exit {rc} without {exp['fail']!r}\n{tail}"
    detail = f"stopped as expected ({exp['fail']})"
    if "stationary" in exp:
        m = STATIONARY.findall(text)
        if not m:
            return "FAIL", f"no last stationary point reported\n{tail}"
        e = float(m[-1][0])
        l = float(m[-1][1].replace("D", "E"))
        re_, rl = exp["stationary"]
        if abs(e - re_) > STAB_TOL_E or abs(l - rl) > STAB_TOL_L:
            return "FAIL", (f"last stationary point E {e:.10f} l {l:.4e}, "
                            f"ref E {re_:.10f} l {rl:.4e}")
        detail += f"; stationary point E {e:.10f} l {l:.4e}"
    return "PASS", detail


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--jobs", type=int, default=4,
                        help="cases to run concurrently (default 4)")
    parser.add_argument("--case", action="append", dest="cases",
                        help="run only these cases (repeatable)")
    parser.add_argument("--long", action="store_true",
                        help="also run the cases marked long")
    parser.add_argument("--workdir",
                        help="where the work directories go (default: a "
                             "new directory under $TMPDIR)")
    parser.add_argument("--keep", action="store_true",
                        help="keep the work directories of passing cases")
    args = parser.parse_args()

    if shutil.which("xaces2") is None:
        sys.exit("xaces2 not found on PATH -- module load intel/2025.1.0 "
                 "openmpi/5.0.7 first, and put ~/Develop/bin on PATH")

    env = os.environ.copy()
    ld = env.get("LD_LIBRARY_PATH", "")
    if fast._EXTRA_LD_PATH not in ld:
        env["LD_LIBRARY_PATH"] = fast._EXTRA_LD_PATH + (":" + ld if ld else "")
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")

    all_cases = sorted(d for d in os.listdir(CASES_DIR)
                       if os.path.isfile(os.path.join(CASES_DIR, d, "EXPECT")))
    expects = {c: load_expect(c) for c in all_cases}
    if args.cases:
        unknown = set(args.cases) - set(all_cases)
        if unknown:
            sys.exit(f"unknown case(s): {', '.join(sorted(unknown))}")
        cases = args.cases
    else:
        cases = [c for c in all_cases if args.long or not expects[c]["long"]]

    if args.workdir:
        workroot = os.path.abspath(args.workdir)
        os.makedirs(workroot, exist_ok=True)
    else:
        workroot = tempfile.mkdtemp(prefix="pccd_slow_")

    print(f"Running {len(cases)} case(s), --jobs {args.jobs}, "
          f"work in {workroot}\n")
    print(f"{'CASE':<24} {'TIME':>9}  {'STATUS':<6}  DETAIL")
    print("-" * 80, flush=True)

    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futures = [ex.submit(run_case, c, expects[c], env, workroot, args.keep)
                   for c in cases]
        for fut in concurrent.futures.as_completed(futures):
            name, status, elapsed, detail = fut.result()
            results.append((name, status))
            lines = detail.splitlines()
            print(f"{name:<24} {elapsed:>8.1f}s  {status:<6}  {lines[0]}")
            if status not in ("PASS",):
                for line in lines[1:]:
                    print(f"    {line}")
            sys.stdout.flush()

    counts = {s: sum(1 for _, st in results if st == s)
              for s in ("PASS", "SLOW", "FAIL", "ERROR")}
    print("-" * 80)
    print(f"{len(cases)} case(s), {counts['PASS']} passed, "
          f"{counts['SLOW']} slow, {counts['FAIL']} failed, "
          f"{counts['ERROR']} errored")
    sys.exit(0 if counts["FAIL"] == 0 and counts["ERROR"] == 0 else 1)


if __name__ == "__main__":
    main()
