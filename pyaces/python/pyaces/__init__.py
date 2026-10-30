"""Python-side API for ACES II through the aces2py extension module.

Two ways to drive a calculation:

* ``pyaces.run(...)`` runs one calculation in a child process and returns the
  requested JOBARC record. An ACES II error (which calls exit() and would take
  the interpreter with it) or a crash becomes a ``CalculationError`` here
  instead. Use this unless you need to look inside the session.

* The in-process functions (``init``, ``run_scf``, ``run_cc``, ...,
  ``get_record``, ``put_record``, ``get_list``) call aces2py directly in this
  interpreter. They are thin, lower-case wrappers over the raw f2py functions
  (``aces2py.Init``, ``aces2py.Runcc``, ``aces2py.Dgetrecpy``, ...), which stay
  available unchanged as ``pyaces.ext()``. The ACES II rules still apply: run
  one calculation per process, work inside the directory holding ZMAT and
  GENBAS, and expect any ACES II error to terminate the process.
"""

from ._ext import ext
from .records import (record_length, get_record, put_record,
                      get_string_record, put_string_record)
from .lists import (list_shape, get_list, put_list, get_all_irreps,
                    put_all_irreps)
from .integrals import get_1e_ints
from .drivers import (init, run_ints, run_scf, run_cc, run_hbar, run_props,
                      run_pccd, run_ee, run_dirmp2, run_intgrt,
                      close_jobarc, open_jobarc)
from .runner import run, Result, CalculationError, ENTRY_POINTS, iuhf_from_zmat
from . import moio

__all__ = [
    'ext', 'run', 'Result', 'CalculationError', 'ENTRY_POINTS', 'iuhf_from_zmat',
    'record_length', 'get_record', 'put_record', 'get_string_record',
    'put_string_record', 'list_shape', 'get_list', 'put_list',
    'get_all_irreps', 'put_all_irreps', 'get_1e_ints', 'init', 'run_ints',
    'run_scf', 'run_cc', 'run_hbar',
    'run_props', 'run_pccd', 'run_ee', 'run_dirmp2', 'run_intgrt',
    'close_jobarc', 'open_jobarc', 'moio',
]
