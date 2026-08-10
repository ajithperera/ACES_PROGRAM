#!/usr/bin/env python3
"""Run ONE ACES II test-suite case through pyaces, in its own process (matches
the proven smoke-test recipe -- aces2py's global COMMON-block state has only
ever been exercised for one job per process; running multiple different
ZMATs in one long-lived process is untested and not worth risking here).

Usage: python3 -u pyaces_run_case.py <workdir> <entrypoint> <iuhf> <record>

entrypoint is one of Runscf / Runcc / Runpccd / Runee.
Runee is self-contained (calls Joda_()/aces_init() itself, same as Init()
does) -- skip Init() for it, call Runee directly. All others: Init() first.

FNO (Frozen Natural Orbitals, FNO_KEEP=... in the ZMAT) is special-cased for
Runcc: the classic driver's `Runfno` (aces2/runfno.f) runs a whole
preliminary Vtran/Intproc/Vcc/Fno sub-job (with overridden IFLAGS) just to
pick which virtuals to freeze, before the REAL production Vtran/Vcc runs on
the reduced space. Doing that preliminary pass via aces2py's OWN in-process
Vtran/Vcc would mean calling them twice in one process -- confirmed unsafe
this session (project_pyaces.md memory: silently converges to a wrong
answer, not just crashes). So the preliminary pass is run here as genuine
OS subprocesses against the already-built standalone xvtran/xintprc/xvcc/
xfno executables (exactly matching how the classic driver itself works --
job_control.F shells out to these same executables for EVERY case, FNO or
not). aces2py's own Vtran/Vcc are still only ever invoked ONCE per process,
via the final Runcc(iuhf) call below, same as every other case.
"""
import os
import re
import subprocess
import sys
import traceback

WORKDIR, ENTRYPOINT, IUHF, RECORD = sys.argv[1:5]
IUHF = int(IUHF)

REPO_DIR = os.path.dirname(os.path.abspath(__file__))
STANDALONE_BIN = os.path.join(os.path.dirname(REPO_DIR), "aces2", "bin")

os.chdir(WORKDIR)
sys.path.insert(0, os.path.join(REPO_DIR, "Main", "build"))

try:
    import aces2py
except Exception:
    print("PYACES_CASE_RESULT: IMPORT_FAILED")
    traceback.print_exc()
    sys.exit(1)

import numpy as np


def _zmat_has_fno():
    try:
        with open("ZMAT") as f:
            text = f.read()
    except OSError:
        return False
    return re.search(r"FNO_KEEP\s*=", text, re.IGNORECASE) is not None


def run_fno_preliminary_pass():
    """Python port of aces2/runfno.f's orchestration, using real subprocess
    calls instead of in-process Fortran calls for the risky double-pass
    part. See module docstring above for why."""
    print("PYACES_CASE: FNO detected, running preliminary NO-generation pass")

    # NOTE: unlike runfno.f (called from job_control.F at a point where the
    # session is NOT yet open), aces2py's JOBARC session is already open
    # here (Runscf() just ran and deliberately never closes it). Calling
    # Acesjainit() on an already-open session is a hard aces_exit, not a
    # no-op (confirmed empirically) -- do NOT call it here.

    iflags = np.zeros(100, dtype=np.int64)
    aces2py.IGetrecpy("IFLAGS  ", iflags)
    icalc = int(iflags[1])          # IFLAGS(2)
    iabcd = int(iflags[92])         # IFLAGS(93)
    ivtran = int(iflags[82])        # IFLAGS(83)
    igamma_abcd = int(iflags[99])   # IFLAGS(100)
    iref = int(iflags[10])          # IFLAGS(11), reference type (0=RHF)

    iflags[1] = 1
    iflags[92] = 0
    iflags[82] = 0
    iflags[99] = 0
    aces2py.Iputrecpy("IFLAGS  ", iflags)

    ndropa = np.zeros(1, dtype=np.int64)
    aces2py.IGetrecpy("NUMDROPA", ndropa)
    aces2py.Iputrecpy("FNOFREEZ", ndropa)
    zero1 = np.zeros(1, dtype=np.int64)
    aces2py.Iputrecpy("NUMDROPA", zero1)
    if iref != 0:
        aces2py.Iputrecpy("NUMDROPB", zero1)

    # Close the JOBARC session before shelling out -- matches runfno.f's own
    # aces_ja_fin-before-runit() discipline, so the subprocess sees a fully
    # flushed, consistent on-disk JOBARC.
    aces2py.Acesjafin()

    for exe in ("xvtran", "xintprc", "xvcc", "xfno"):
        logname = f"fno_prelim_{exe}.log"
        r = subprocess.run([os.path.join(STANDALONE_BIN, exe)],
                            capture_output=True, text=True)
        with open(logname, "w") as f:
            f.write(r.stdout)
            f.write(r.stderr)
        print(f"PYACES_CASE: FNO prelim {exe} rc={r.returncode} (see {logname})")
        if r.returncode != 0:
            raise RuntimeError(
                f"FNO preliminary pass: {exe} failed, rc={r.returncode}, "
                f"see {logname}")

    aces2py.Acesjainit()
    iflags = np.zeros(100, dtype=np.int64)
    aces2py.IGetrecpy("IFLAGS  ", iflags)
    iflags[1] = icalc
    iflags[92] = iabcd
    iflags[82] = ivtran
    iflags[99] = igamma_abcd
    aces2py.Iputrecpy("IFLAGS  ", iflags)
    aces2py.Acesjafin()

    if not os.path.exists("IIII"):
        print("PYACES_CASE: FNO prelim pass consumed IIII, regenerating")
        for exe in ("xvmol", "xvmol2ja"):
            logname = f"fno_prelim_{exe}.log"
            r = subprocess.run([os.path.join(STANDALONE_BIN, exe)],
                                capture_output=True, text=True)
            with open(logname, "w") as f:
                f.write(r.stdout)
                f.write(r.stderr)
            print(f"PYACES_CASE: FNO prelim {exe} rc={r.returncode} (see {logname})")
            if r.returncode != 0:
                raise RuntimeError(
                    f"FNO preliminary pass: {exe} regeneration failed, "
                    f"rc={r.returncode}, see {logname}")

    # Defensive reopen before handing back to Runcc() -- don't rely on any
    # lazy-reopen behavior in the low-level Getrec/Putrec primitives.
    aces2py.Acesjainit()
    print("PYACES_CASE: FNO preliminary pass complete")


try:
    if ENTRYPOINT != "Runee":
        print("PYACES_CASE: calling Init()")
        aces2py.Init()

    if ENTRYPOINT == "Runcc" and _zmat_has_fno():
        print(f"PYACES_CASE: calling Runscf({IUHF}) [FNO pre-stage]")
        aces2py.Runscf(IUHF)
        print("PYACES_CASE: Runscf returned OK")
        run_fno_preliminary_pass()

    fn = getattr(aces2py, ENTRYPOINT)
    print(f"PYACES_CASE: calling {ENTRYPOINT}({IUHF})")
    fn(IUHF)
    print(f"PYACES_CASE: {ENTRYPOINT} returned OK")
except Exception:
    print(f"PYACES_CASE_RESULT: {ENTRYPOINT.upper()}_FAILED")
    traceback.print_exc()
    sys.exit(1)

try:
    dest = np.zeros(1, dtype=np.float64)
    aces2py.Dgetrecpy(RECORD, dest)
    print(f"PYACES_CASE_VALUE: {dest[0]!r}")
    if dest[0] == 0.0:
        print("PYACES_CASE_RESULT: RECORD_ZERO_OR_MISSING")
        sys.exit(1)
    print(f"PYACES_CASE_RESULT: SUCCESS {dest[0]:.10f}")
except Exception:
    print("PYACES_CASE_RESULT: DGETRECPY_FAILED")
    traceback.print_exc()
    sys.exit(1)
