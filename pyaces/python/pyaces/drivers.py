"""In-process calculation drivers: lower-case names for aces2py's Run* entry points.

iuhf is 0 for RHF and 1 for UHF/ROHF. Call init() first, except before
run_ee(), which reads ZMAT itself. Do not call the same
stage twice in one process (Vtran/Vcc twice has converged to wrong answers),
and do not call run_hbar() after run_cc() -- it would redo Vtran.
"""

from ._ext import ext


def init():
    """Read ZMAT (Joda) and set up core memory, flags and lists."""
    ext().Init()


def run_ints(iuhf):
    ext().Runints(iuhf)


def run_scf(iuhf):
    ext().Runscf(iuhf)


def run_cc(iuhf):
    ext().Runcc(iuhf)


def run_hbar(iuhf):
    ext().Runhbar(iuhf)


def run_props(iuhf):
    ext().Runprops(iuhf)


def run_pccd(iuhf):
    ext().Runpccd(iuhf)


def run_ee(iuhf):
    ext().Runee(iuhf)


def run_dirmp2(iuhf):
    ext().Rundirmp2(iuhf)


def run_intgrt(iuhf):
    ext().Runintgrt(iuhf)


def close_jobarc():
    """Flush and close JOBARC (e.g. before an external program reads it)."""
    ext().Acesjafin()


def open_jobarc():
    """Reopen JOBARC after close_jobarc(); reopening an open session is fatal."""
    ext().Acesjainit()
