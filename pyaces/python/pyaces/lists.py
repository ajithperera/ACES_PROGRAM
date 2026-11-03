"""MO-list (MOIO) access for the current in-process session.

List numbers are in pyaces.moio. A list is addressed as MOIO(irrep, list):
aces2py.Getlstpy/Putlstpy take the irrep of the columns as Left_index and the
list number as Right_index; their Irrep argument is not used. get_list and
put_list check the request against list_shape() and raise ValueError, where
aces2py would only print a message and do nothing.
"""

import numpy as np

from ._ext import ext


def list_shape(lst, irrep):
    """(rows, columns) of MO list lst in irrep irrep, or None if it does not exist."""
    nrows, ncols = ext().Getlstdimpy(irrep, lst)
    if nrows < 1 or ncols < 1:
        return None
    return int(nrows), int(ncols)


def _check(lst, irrep, first_column, n_columns, size):
    shape = list_shape(lst, irrep)
    if shape is None:
        raise KeyError(f'MO list ({irrep}, {lst}) does not exist')
    nrows, ncols = shape
    if first_column < 1 or n_columns < 1 or first_column + n_columns - 1 > ncols:
        raise ValueError(f'columns {first_column}..{first_column + n_columns - 1} '
                         f'requested; list ({irrep}, {lst}) has {ncols}')
    if size is not None and size < nrows * n_columns:
        raise ValueError(f'size {size} is too small for {n_columns} columns '
                         f'of {nrows} rows')
    return nrows


def get_list(lst, irrep, first_column, n_columns, size=None):
    """Read n_columns columns of MO list lst, irrep irrep, into a float64 array.

    size defaults to exactly the columns read (rows * n_columns).
    """
    nrows = _check(lst, irrep, first_column, n_columns, size)
    if size is None:
        size = nrows * n_columns
    out = np.zeros(size, dtype=np.float64)
    ext().Getlstpy(first_column, n_columns, 0, irrep, lst, out, size)
    return out


def put_list(values, lst, irrep, first_column, n_columns):
    """Write n_columns columns of MO list lst, irrep irrep."""
    a = np.ascontiguousarray(values, dtype=np.float64).ravel()
    _check(lst, irrep, first_column, n_columns, a.size)
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
