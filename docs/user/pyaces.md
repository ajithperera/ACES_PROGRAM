# pyaces: ACES II from Python

pyaces builds the ACES II programs into one Python extension module,
`aces2py`, with F2PY. A Python script can then run the steps of a
calculation in-process and read and write `JOBARC` records and lists as
NumPy arrays: for example, to inspect intermediate quantities, or to pass
integrals and amplitudes to other Python codes.

## Building and importing

The superbuild builds pyaces after ACES II (see
{doc}`../getting_started/index`). The module ends up in
`pyaces/Main/build/` as `aces2py.cpython-*.so`. Load the same modules used
to build it, and put that directory on the Python path:

```python
import sys
sys.path.insert(0, "/path/to/ACES_PROGRAM/pyaces/Main/build")
import aces2py
```

## A first calculation

As with `xaces2`, the job runs in a directory that holds `ZMAT` and
`GENBAS`. This script runs RHF-CCSD for water in STO-3G:

```text
Water
O
H 1 R
H 1 R 2 A

R=0.96
A=104.5

*ACES2
BASIS=STO-3G,CALC=CCSD

```

```python
import sys
sys.path.insert(0, "/path/to/ACES_PROGRAM/pyaces/Main/build")
import numpy as np
import aces2py

aces2py.Init()          # read ZMAT (the xjoda step) and open JOBARC
aces2py.Runscf(0)       # SCF; the argument is 0 for RHF, 1 for UHF
aces2py.Runcc(0)        # transformation, integral sort, CC
e = np.zeros(1)
aces2py.Dgetrecpy("TOTENERG", e)
print("total energy", e[0])
```

It prints −75.0130370, the same as `xaces2` on the same input to within
10<sup>-7</sup> hartree.

```{important}
Run **one calculation per Python process**. The ACES II programs keep their
state in Fortran common blocks, and running a second job, or the same step
twice, in one process has been seen to converge to wrong answers. Start a
new process (for example with `subprocess`) for each job.
```

## Functions

Running steps (the argument `iuhf` is 0 for RHF, 1 for UHF/ROHF):

| Function | Does |
|---|---|
| `Init()` | Parse `ZMAT`, set up `JOBARC`, initialize ACES; call first |
| `Runints(iuhf)` | AO integrals |
| `Runscf(iuhf)` | Integrals and SCF |
| `Runcc(iuhf)` | Integral transformation and sorting, then the correlated method in `CALC` |
| `Runpccd(iuhf)` | Single-point pCCD: integrals, SCF, transformation and pCCD (orbital-optimized pCCD, `PROGRAM=PCCD`, is not available) |
| `Runhbar(iuhf)`, `Buildhbar()` | Similarity-transformed Hamiltonian |
| `Runprops(iuhf)` | Properties |
| `Runee(iuhf)` | EOM-CC excitation energies; self-contained, so call it **instead of** `Init()` |
| `Acesjainit()`, `Acesjafin()` | Open and close `JOBARC` (for example around running a stand-alone program with `subprocess`) |

Records and lists (array arguments are NumPy arrays that the function fills
or reads; record names are 8 characters, padded with blanks):

| Function | Does |
|---|---|
| `Dgetrecpy(name, dest)`, `Dputrecpy(name, src)` | Read or write a double-precision `JOBARC` record |
| `IGetrecpy(name, dest)`, `Iputrecpy(name, src)` | The same for integer records (`dtype=np.int64`) |
| `Cgetrecpy(name, length, dest)`, `Cputrecpy(...)` | Character records |
| `Getreclenpy(name, length)` | Length of a record |
| `Getlstpy(first, ncols, irrep, left, right, dest, ndim)`, `Putlstpy(...)` | Columns of a list: `left` is the irrep of the columns, `right` the list number, `irrep` is not used (see {doc}`../programmer/io`) |
| `Getlstdimpy(left, right)` | Rows and columns of a list (0, 0 if it does not exist) |
| `Getallpy(dest, irrep, right)`, `Putallpy(...)` | A whole list |
| `Get1ehpy(name, dest, n)`, `Get2ehpy(name, dest, n)` | One- and two-electron integrals as dense arrays (with symmetry, the one-electron matrix is in the SO basis, block diagonal by irrep) |

For example, the job flags:

```python
iflags = np.zeros(100, dtype=np.int64)
aces2py.IGetrecpy("IFLAGS  ", iflags)
```

## Test suite

`pyaces/pyaces_suite_runner.py` runs the ACES II test inputs through
`aces2py` (one process per case, through `pyaces_run_case.py`) and compares
the energies with the ACES II references; see {doc}`../programmer/testing`.
`pyaces_run_case.py` is also the best worked example of driving each entry
point, including a frozen-natural-orbital job that mixes `aces2py` calls
with stand-alone ACES II programs.
