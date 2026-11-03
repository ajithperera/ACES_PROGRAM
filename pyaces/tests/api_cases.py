#!/usr/bin/env python3
"""API tests for the in-process pyaces functions, one calculation per test.

run_api_tests.py runs each test in its own process (an ACES II error ends the
process, and one calculation per process is the rule anyway) inside a fresh
directory holding the test's ZMAT and a GENBAS link:

    python3 api_cases.py <test name>

prints API_TEST_RESULT: PASS, or a traceback and API_TEST_RESULT: FAIL.

The reference energies come from the classic xaces2 on the same inputs (oneAPI
2024.0, 2026-10-30), not from pyaces.
"""

import os
import sys
import traceback

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'python'))
import pyaces  # noqa: E402

_H2O = """O
H 1 R
H 1 R 2 A

R=0.96
A=104.5

"""

ZMATS = {
    'rhf_scf': 'H2O RHF SCF, STO-3G\n' + _H2O +
               '*ACES2(BASIS=STO-3G,CALC=SCF\nREF=RHF,MEM_SIZE=1gb,scf_conv=10)\n\n',
    'uhf_scf': 'H2O+ UHF SCF, STO-3G\n' + _H2O +
               '*ACES2(BASIS=STO-3G,CALC=SCF\nCHARGE=1,MULTIPLICITY=2\n'
               'REF=UHF,MEM_SIZE=1gb,scf_conv=10)\n\n',
    'ccsd': 'H2O CCSD, STO-3G\n' + _H2O +
            '*ACES2(BASIS=STO-3G,CALC=CCSD\nREF=RHF,MEM_SIZE=1gb,cc_conv=10,scf_conv=10)\n\n',
    'ccsd_c1': 'H2O CCSD, STO-3G, no symmetry\n' + _H2O +
               '*ACES2(BASIS=STO-3G,CALC=CCSD,SYMMETRY=OFF\n'
               'REF=RHF,MEM_SIZE=1gb,cc_conv=10,scf_conv=10)\n\n',
}

# xaces2 on the ZMATS above.
E_SCF_RHF = -74.963319004820
E_SCF_UHF = -74.6567025644
E_CCSD = -75.013036960911
TOL = 1e-8


def check(cond, what):
    if not cond:
        raise AssertionError(what)


def close(got, ref, what, tol=TOL):
    check(abs(got - ref) <= tol, f'{what}: got {got!r}, expected {ref!r} (tol {tol})')


def scalar(name, dtype=np.float64):
    return pyaces.get_record(name, dtype, 1)[0]


# --- records ---------------------------------------------------------------

def test_records():
    """Record length, double/integer/string round trips, missing records."""
    pyaces.init()
    ext = pyaces.ext()

    # Getreclenpy returns the length (before f97a9759 it returned nothing).
    # Joda's records: three atoms, nine coordinates.
    n = ext.Getreclenpy('COORD   ')
    check(isinstance(n, (int, np.integer)), f'Getreclenpy returned {n!r}, not an int')
    check(n == 9, f'COORD length {n}')
    check(pyaces.record_length('COORD') == 9, 'record_length(COORD)')
    check(int(scalar('NATOMS', np.int64)) == 3, 'NATOMS for H2O is 3')

    check(pyaces.record_length('NOSUCHRC') is None, 'missing record has a length')
    try:
        pyaces.get_record('NOSUCHRC')
        raise AssertionError('get_record of a missing record did not raise')
    except KeyError:
        pass

    d = np.array([1.5, -2.25, 3.0e-12, np.pi])
    pyaces.put_record('PYTESTD1', d)
    check(pyaces.record_length('PYTESTD1') == d.size, 'double record length')
    check(np.array_equal(pyaces.get_record('PYTESTD1'), d), 'double record round trip')

    i = np.array([7, -3, 2**40], dtype=np.int64)
    pyaces.put_record('PYTESTI1', i)
    check(np.array_equal(pyaces.get_record('PYTESTI1', np.int64), i), 'integer record round trip')

    # Cgetrecpy returns the string (before f97a9759 Dest was never handed back).
    s = b'pyaces API test 0123456789'
    pyaces.put_string_record('PYTESTC1', s)
    raw = ext.Cgetrecpy('PYTESTC1', len(s))
    check(raw is not None, 'Cgetrecpy returned None')
    check(bytes(raw)[:len(s)] == s, f'Cgetrecpy returned {bytes(raw)[:len(s)]!r}')
    check(pyaces.get_string_record('PYTESTC1', len(s)) == s, 'string record round trip')
    # Counted in 8-byte words, as aces2's putcrec does (it was 7 words, then
    # 4-byte ones, and only the first 8 characters were stored).
    check(pyaces.record_length('PYTESTC1') == 4, f"PYTESTC1 is {pyaces.record_length('PYTESTC1')} words")

    # Records survive closing and reopening JOBARC.
    pyaces.close_jobarc()
    pyaces.open_jobarc()
    check(np.array_equal(pyaces.get_record('PYTESTD1'), d), 'double record after reopen')


# --- drivers and the iuhf argument -------------------------------------------

def test_scf_rhf():
    pyaces.init()
    pyaces.run_scf(0)
    close(scalar('TOTENERG'), E_SCF_RHF, 'RHF SCF energy')
    check(pyaces.record_length('HBAR') is None, 'HBAR record before Runhbar')
    # A character record written by aces2 (vscf's orbital irreps, 8 characters
    # per orbital). pyaces used to ask for twice its length and end the process.
    n = pyaces.record_length('EVCSYMAF')
    check(n == 7, f'EVCSYMAF is {n} words')
    labels = pyaces.get_string_record('EVCSYMAF', 8 * n).split()
    check(labels == [b'A1', b'A1', b'B2', b'A1', b'B1', b'A1', b'B2'], f'EVCSYMAF = {labels}')


def test_scf_uhf():
    """iuhf=1 must reach the Fortran side as 1 (it was a 4-byte int before f97a9759)."""
    pyaces.init()
    pyaces.run_scf(1)
    close(scalar('TOTENERG'), E_SCF_UHF, 'UHF SCF energy')
    check(pyaces.record_length('SCFEVCB0') is not None, 'no beta orbitals: UHF was not run')


# --- Runhbar ---------------------------------------------------------------
# Runhbar takes one of three branches. "SCF done" (SCFEVCA0 exists) and "from
# scratch" are tested below. The "CC done" branch needs a CC_DONE record and
# no SCFEVCA0, and no ACES II program writes CC_DONE, so it cannot be reached
# through the API. Before f97a9759 both branches called the COMMON /HBAR/
# instead of the subroutine and crashed.

def _hbar_checks():
    close(scalar('TOTENERG'), E_CCSD, 'CCSD energy')
    # hbar writes the HBAR record when it has built the elements.
    check(pyaces.record_length('HBAR') is not None, 'no HBAR record: hbar did not finish')


def test_hbar_from_scratch():
    pyaces.init()
    pyaces.run_hbar(0)
    _hbar_checks()


def test_hbar_after_scf():
    pyaces.init()
    pyaces.run_scf(0)
    close(scalar('TOTENERG'), E_SCF_RHF, 'SCF energy before Runhbar')
    pyaces.run_hbar(0)
    _hbar_checks()


# --- MO lists --------------------------------------------------------------

def test_lists():
    """list_shape, get_list column by column vs get_all_irreps, errors, put_list."""
    pyaces.init()
    pyaces.run_cc(0)
    close(scalar('TOTENERG'), E_CCSD, 'CCSD energy')
    nirrep = int(scalar('NIRREP', np.int64))
    check(nirrep == 4, f'C2v expected, NIRREP={nirrep}')

    lst = 46   # T2(Ij,Ab)
    parts = []
    for irr in range(1, nirrep + 1):
        shape = pyaces.list_shape(lst, irr)
        if shape is None:
            continue
        nr, nc = shape
        a = pyaces.get_list(lst, irr, 1, nc)
        check(a.size == nr * nc, 'get_list size')
        if nc > 1:
            check(np.array_equal(pyaces.get_list(lst, irr, 2, 1), a[nr:2 * nr]),
                  f'column 2 of irrep {irr}')
        parts.append(a)
    whole = np.concatenate(parts)
    check(np.array_equal(whole, pyaces.get_all_irreps(lst, 1, whole.size)),
          'irrep-by-irrep get_list differs from get_all_irreps')
    check(np.abs(whole).max() > 1e-3, 'T2 amplitudes are all zero')

    check(pyaces.list_shape(lst, 9) is None, 'irrep 9 exists')
    nr, nc = pyaces.list_shape(lst, 1)
    for args in ((lst, 1, nc, 2), (lst, 1, 0, 1)):
        try:
            pyaces.get_list(*args)
            raise AssertionError(f'get_list{args} did not raise')
        except ValueError:
            pass
    try:
        pyaces.get_list(lst, 1, 1, 1, size=nr - 1)
        raise AssertionError('get_list with a short size did not raise')
    except ValueError:
        pass
    try:
        pyaces.get_list(lst, 9, 1, 1)
        raise AssertionError('get_list of a missing list did not raise')
    except KeyError:
        pass
    buf = np.full(nr - 1, 7.0)
    pyaces.ext().Getlstpy(1, 1, 0, 1, lst, buf, buf.size)
    check((buf == 7.0).all(), 'Getlstpy wrote into a too-short array')

    a = pyaces.get_list(lst, 1, 1, 1)
    pyaces.put_list(2 * a, lst, 1, 1, 1)
    check(np.array_equal(pyaces.get_list(lst, 1, 1, 1), 2 * a), 'put_list did not write')
    pyaces.put_list(a, lst, 1, 1, 1)
    check(np.array_equal(pyaces.get_list(lst, 1, 1, 1), a), 'put_list did not restore')


# --- one-electron integrals ------------------------------------------------

def _one_electron(expect_nirrep):
    pyaces.init()
    pyaces.run_scf(0)
    nirrep = int(scalar('NIRREP', np.int64))
    check(nirrep == expect_nirrep, f'NIRREP={nirrep}')
    n = int(scalar('NBASTOT', np.int64))
    S = pyaces.get_1e_ints('overlap')
    T = pyaces.get_1e_ints('kinetic')
    H = pyaces.get_1e_ints('oneh')
    for M, name in ((S, 'S'), (T, 'T'), (H, 'H')):
        check(M.shape == (n, n), f'{name} shape')
        check(np.array_equal(M, M.T), f'{name} not symmetric')
    nbfirr = pyaces.get_record('NUMBASIR', np.int64, nirrep)
    off = np.concatenate(([0], np.cumsum(nbfirr)))
    inside = np.zeros((n, n), bool)
    for k in range(nirrep):
        inside[off[k]:off[k + 1], off[k]:off[k + 1]] = True
    check(not S[~inside].any() and not H[~inside].any(), 'nonzero outside the irrep blocks')
    check(np.all(np.linalg.eigvalsh(T) > 0), 'kinetic matrix not positive definite')

    # The SCF eigenvectors are in the same basis: C^T S C = 1.
    C = pyaces.get_record('SCFEVCA0').reshape(n, n, order='F')
    check(np.abs(C.T @ S @ C - np.eye(n)).max() < 1e-12, 'C^T S C is not 1')

    # RHF: E(SCF) = sum over occupied (h_ii + e_i) + E_nuc, against xaces2's E(SCF).
    e = pyaces.get_record('SCFEVLA0', length=n)
    occ = np.argsort(e)[:5]
    hmo = np.einsum('pi,pq,qi->i', C[:, occ], H, C[:, occ])
    close(hmo.sum() + e[occ].sum() + scalar('NUCREP'), E_SCF_RHF,
          'E(SCF) rebuilt from Get1ehpy oneh')


def test_one_electron_c2v():
    """Block-diagonal SO matrices; STO-3G H2O has an empty irrep (A2)."""
    _one_electron(4)


def test_one_electron_c1():
    _one_electron(1)


# --- table -----------------------------------------------------------------

# name -> (ZMAT, function)
TESTS = {
    'records': ('rhf_scf', test_records),
    'scf_rhf': ('rhf_scf', test_scf_rhf),
    'scf_uhf': ('uhf_scf', test_scf_uhf),
    'hbar_from_scratch': ('ccsd', test_hbar_from_scratch),
    'hbar_after_scf': ('ccsd', test_hbar_after_scf),
    'lists': ('ccsd', test_lists),
    'one_electron_c2v': ('ccsd', test_one_electron_c2v),
    'one_electron_c1': ('ccsd_c1', test_one_electron_c1),
}


if __name__ == '__main__':
    name = sys.argv[1]
    try:
        TESTS[name][1]()
    except Exception:
        traceback.print_exc()
        sys.stdout.flush()
        print('API_TEST_RESULT: FAIL', flush=True)
        sys.exit(1)
    print('API_TEST_RESULT: PASS', flush=True)
