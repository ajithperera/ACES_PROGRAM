"""MO-list (MOIO) access for the current in-process session.

List numbers are in pyaces.moio. Arguments are named for what Getlstpy
actually does with them: aces2py.Getlstpy's Left_index is the irrep and
Right_index the list number, and its Irrep/Ndim arguments are ignored.
"""

import numpy as np

from ._ext import ext


def get_list(lst, irrep, first_column, n_columns, size):
    """Read n_columns columns of MO list lst, irrep irrep, into a float64 array of length size."""
    out = np.zeros(size, dtype=np.float64)
    ext().Getlstpy(first_column, n_columns, 0, irrep, lst, out, size)
    return out


def put_list(values, lst, irrep, first_column, n_columns):
    """Write n_columns columns of MO list lst, irrep irrep."""
    a = np.ascontiguousarray(values, dtype=np.float64).ravel()
    ext().Putlstpy(a, first_column, n_columns, 0, irrep, lst)


def get_all_irreps(lst, irrep, size):
    """Read every irrep block of list lst with total symmetry irrep."""
    out = np.zeros(size, dtype=np.float64)
    ext().Getallpy(out, irrep, lst)
    return out


def put_all_irreps(values, lst, irrep):
    """Write every irrep block of list lst with total symmetry irrep."""
    a = np.ascontiguousarray(values, dtype=np.float64).ravel()
    ext().Putallpy(a, irrep, lst)
