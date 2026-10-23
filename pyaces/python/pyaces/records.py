"""JOBARC record access for the current in-process session."""

import numpy as np

from ._ext import ext, recname


def record_length(name):
    """Length of a record in integer words, or None if it does not exist."""
    n = int(ext().Getreclenpy(recname(name)))
    return None if n < 0 else n


def get_record(name, dtype=np.float64, length=None):
    """Read a record and return it as a new 1-D array.

    dtype is np.float64 or np.int64 (one 8-byte word per element). The length
    defaults to the record's own length. A missing record raises KeyError
    here; asking aces2py for it directly would end the process.
    """
    dtype = np.dtype(dtype)
    if length is None:
        length = record_length(name)
        if length is None:
            raise KeyError(f'JOBARC has no record {name!r}')
    out = np.zeros(length, dtype=dtype)
    if dtype == np.float64:
        ext().Dgetrecpy(recname(name), out)
    elif dtype == np.int64:
        ext().IGetrecpy(recname(name), out)
    else:
        raise TypeError(f'records are float64 or int64, not {dtype}')
    return out


def put_record(name, values):
    """Write a float64 or int64 record (the dtype of values decides which)."""
    a = np.ascontiguousarray(values)
    if a.dtype.kind == 'f':
        ext().Dputrecpy(recname(name), a.astype(np.float64, copy=False).ravel())
    elif a.dtype.kind in 'iu':
        ext().Iputrecpy(recname(name), a.astype(np.int64, copy=False).ravel())
    else:
        raise TypeError(f'records are float64 or int64, not {a.dtype}')


def get_string_record(name, nbytes):
    """Read nbytes (at most 1024) of a character record as bytes."""
    return bytes(ext().Cgetrecpy(recname(name), nbytes))[:nbytes]


def put_string_record(name, data):
    """Write a character record; data is bytes or str."""
    if isinstance(data, str):
        data = data.encode()
    buf = np.array(data, dtype=f'S{len(data)}')
    ext().Cputrecpy(recname(name), len(data), buf)
