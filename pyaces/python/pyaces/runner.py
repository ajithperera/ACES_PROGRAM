"""Run one ACES II calculation in a child process.

aces_exit() ends the whole process, so an ACES II error inside the caller's
interpreter cannot be caught. run() starts pyaces_run_case.py -- the same
worker the test suite uses, including its FNO handling -- in a fresh process
and turns its result markers into a return value or a CalculationError.
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass

from ._ext import ext_dir

_PYACES_DIR = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
_WORKER = os.path.join(_PYACES_DIR, 'pyaces_run_case.py')

# method name -> aces2py entry point
ENTRY_POINTS = {
    'scf': 'Runscf', 'cc': 'Runcc', 'ee': 'Runee', 'pccd': 'Runpccd',
    'dirmp2': 'Rundirmp2', 'ks': 'Runintgrt', 'hbar': 'Runhbar',
    'props': 'Runprops',
}


class CalculationError(RuntimeError):
    """The child process did not produce the record. status is the worker's
    PYACES_CASE_RESULT tag (e.g. RUNCC_FAILED) or CRASH/TIMEOUT."""

    def __init__(self, message, status, returncode, workdir, log_tail):
        super().__init__(message)
        self.status = status
        self.returncode = returncode
        self.workdir = workdir
        self.log_tail = log_tail


@dataclass
class Result:
    value: float
    record: str
    workdir: str
    log: str


def iuhf_from_zmat(text):
    """0 for an RHF reference, 1 otherwise (REF/REFERENCE=UHF or ROHF, or QRHF)."""
    kw = text.upper().replace(' ', '')
    m = re.search(r'\bREF(?:ERENCE)?=(\w+)', kw)
    ref = m.group(1) if m else 'RHF'
    qrhf = re.search(r'\bQRHF_G', kw) is not None and 'MAKERHF' not in kw
    return 0 if ref == 'RHF' and not qrhf else 1


def run(method, zmat, genbas, iuhf=None, record='TOTENERG', workdir=None,
        timeout=None, keep=False):
    """Run one calculation and return Result(value, record, workdir, log).

    method  'scf', 'cc', 'ee', 'pccd', 'dirmp2', 'ks', 'hbar', 'props', or an
            aces2py entry-point name such as 'Runcc'.
    zmat    ZMAT text, or a path to a ZMAT file.
    genbas  path to the GENBAS basis library (symlinked into the work dir).
    iuhf    0 (RHF) or 1 (UHF/ROHF); guessed from ZMAT when None.
    workdir directory to run in (created, must not hold another job); a
            temporary one is used when None and removed afterwards unless
            keep=True or the run fails.
    """
    entry = ENTRY_POINTS.get(method, method)
    if entry not in ENTRY_POINTS.values():
        raise ValueError(f'unknown method {method!r}; use one of {sorted(ENTRY_POINTS)}')
    if os.path.isfile(zmat):
        with open(zmat) as f:
            zmat = f.read()
    if iuhf is None:
        iuhf = iuhf_from_zmat(zmat)

    temporary = workdir is None
    workdir = tempfile.mkdtemp(prefix='pyaces-') if temporary else os.path.abspath(workdir)
    os.makedirs(workdir, exist_ok=True)
    with open(os.path.join(workdir, 'ZMAT'), 'w') as f:
        f.write(zmat)
    link = os.path.join(workdir, 'GENBAS')
    if not os.path.exists(link):
        os.symlink(os.path.abspath(genbas), link)

    env = dict(os.environ, PYACES_EXT_DIR=ext_dir())
    log_path = os.path.join(workdir, 'run.log')
    try:
        with open(log_path, 'w') as log:
            proc = subprocess.run(
                [sys.executable, '-u', _WORKER, workdir, entry, str(iuhf), record],
                stdout=log, stderr=subprocess.STDOUT, env=env, timeout=timeout)
        rc = proc.returncode
    except subprocess.TimeoutExpired:
        rc = None
    with open(log_path, errors='replace') as f:
        text = f.read()

    m = re.search(r'^PYACES_CASE_RESULT: SUCCESS (\S+)', text, re.M)
    v = re.search(r'^PYACES_CASE_VALUE: (?:np\.float64\()?([-+0-9.eE]+)', text, re.M)
    if rc == 0 and m and v:
        result = Result(float(v.group(1)), record, workdir, text)
        if temporary and not keep:
            shutil.rmtree(workdir, ignore_errors=True)
            result.workdir = None
        return result

    tag = re.search(r'^PYACES_CASE_RESULT: (\S+)', text, re.M)
    status = 'TIMEOUT' if rc is None else (tag.group(1) if tag else 'CRASH')
    tail = '\n'.join(text.splitlines()[-25:])
    raise CalculationError(
        f'{entry}(iuhf={iuhf}) failed: {status} (exit {rc}); log in {log_path}',
        status, rc, workdir, tail)
