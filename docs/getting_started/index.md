# Building, testing and a first run

## Requirements

- Intel oneAPI compilers (`ifx`, `icx`, `icpx`) and Intel MKL
- An MPI library with Intel compiler wrappers (OpenMPI 5 is used on HiPerGator); ACES III needs it
- Python 3 with NumPy, for pyaces and the test runners
- CMake 3.18 or newer and GNU make

On HiPerGator the default module set is loaded for you by every build and test step:

```bash
module load intel/2025.1.0 openmpi/5.0.7 python/3.12
```

## Building

The top-level CMake build is a *superbuild*: it does not compile code itself
but drives each part's own build system with one shared toolchain. ACES II
builds with its GNU makefiles, ACES III with autotools, and pyaces with
`pyaces/build.sh`.

```bash
git clone https://github.com/ajithperera/ACES_PROGRAM.git
cd ACES_PROGRAM
cmake -S . -B build
cmake --build build
```

Options, given as `-D<option>=<value>` to the first `cmake` command:

| Option | Default | Meaning |
|---|---|---|
| `ACES_MODULES` | `intel/2025.1.0 openmpi/5.0.7 python/3.12` | Modules loaded before each step. Use `""` to build with the compilers and libraries already in your environment. |
| `ACES_NP` | `4` | MPI ranks and parallel cases used by the test suites |
| `ACES_BUILD_ACES2` | `ON` | Build ACES II |
| `ACES_BUILD_ACES3` | `ON` | Build ACES III |
| `ACES_BUILD_PYACES` | `ON` | Build pyaces (needs ACES II) |

On machines other than HiPerGator, the ACES II makefiles fall back to a
generic Linux/Intel configuration when `ifx` is on `PATH`; set
`ACES_HOST=<name>` to choose a specific host entry.

Where the results go:

- ACES II executables (`xaces2`, `xjoda`, `xvscf`, ...) and libraries:
  `aces2/bin` and `aces2/lib`. Developers who want them elsewhere can create
  an untracked `aces2/Makefiles/local.mk` containing, for example,
  `DIR_INSTBIN := $(HOME)/Develop/bin`.
- ACES III: `aces3/bin` (`xaces3`, and the compiled SIAL programs in `aces3/bin/sio`)
- pyaces: the extension module in `pyaces/`

On HiPerGator, build through Slurm rather than on a login or interactive node.

## Running the test suites

```bash
ctest --test-dir build
```

This runs the ACES II fast suite, the ACES III fast suite, the pyaces suite,
and a `rollup` step that prints one summary table. Each suite reports
every case as PASS, FAIL, ERROR (and, for pyaces, MISMATCH or TIMEOUT)
against stored reference values.

```{important}
Run the ACES III suite inside a **single-node** job (`sbatch --nodes=1`).
Each case keeps its scratch files in a node-local directory, so a job
spread over several nodes fails with timeouts and "Block not found"
errors. The runner refuses to start when `SLURM_JOB_NUM_NODES` is greater
than 1. The full ACES III suite takes about 40 minutes with `ACES_NP=4`.
```

## A first ACES II calculation

ACES II reads its input from a file named `ZMAT` and its basis sets from
`GENBAS` in the current directory. Put `aces2/bin` on your `PATH`, make an
empty working directory, and copy `basis/GENBAS` into it (the basis set
library that the test suites use; `basis/ECPDATA` holds the ECPs):

```text
Water CCSD(T)/DZP single point
H
O 1 R
H 2 R 1 A

R=0.943056
A=105.96879

*ACES2
BASIS=DZP
CALC=CCSD(T),SPHERICAL=OFF
```

```bash
xaces2 > out
```

`xaces2` is a driver: it runs the member executables (`xjoda`, `xvmol`,
`xvscf`, `xvcc`, ...) in turn, and they pass data to each other through
files in the working directory such as `JOBARC` and `MOINTS`. Run each job
in its own directory. {doc}`../user/aces2_input` describes the ZMAT file, and
{ref}`keyword-reference` lists every keyword.
