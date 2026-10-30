"""AO/SO integrals of the current in-process session, from the IIII file."""

import numpy as np

from ._ext import ext, recname

ONE_ELECTRON = ('overlap', 'oneh', 'kinetic')


def get_1e_ints(name):
    """One-electron integrals as a dense (nbas, nbas) array.

    name is 'overlap', 'oneh' (core Hamiltonian) or 'kinetic'. With symmetry
    the basis is the SO basis with the functions in irrep order (the order of
    NUMBASIR), so the matrix is block diagonal; without symmetry it is the AO
    basis. Run init() and run_ints() first.
    """
    if name not in ONE_ELECTRON:
        raise ValueError(f'name is one of {ONE_ELECTRON}, not {name!r}')
    if ext().Getreclenpy(recname('NBASTOT')) < 1:
        raise RuntimeError('NBASTOT is not in JOBARC; run init() first')
    nbas = np.zeros(1, dtype=np.int64)
    ext().IGetrecpy(recname('NBASTOT'), nbas)
    n = int(nbas[0])
    out = np.zeros(n * n, dtype=np.float64)
    ext().Get1ehpy(name, out, n)
    return out.reshape(n, n, order='F')
