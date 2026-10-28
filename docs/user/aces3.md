# Running ACES III

ACES III is the parallel member of the system. Its methods are written in
SIAL (Super Instruction Assembly Language), compiled to `.sio` files, and
executed by the SIP runtime (Super Instruction Processor) under MPI. Data
are held as distributed blocks, so memory and work scale across processes.
ACES III shares the ZMAT input, `GENBAS`, and `JOBARC` with ACES II.

The test suite covers:

- SCF, MBPT(2), CCSD and CCSD(T) energies, RHF and UHF, with and without
  frozen core
- MBPT(2) and CCSD gradients and geometry optimization; CCSD frequencies
- EOM-CCSD excitation energies, ionization potentials and electron
  affinities; double ionization and double attachment (DIP/DEA-EOM-CC);
  STEOM
- dipole and quadrupole moments, static, dynamic and complex
  polarizabilities, and CCSD g-tensors
- explicitly correlated (F12) variants of several of these

## Running a job

Put `ZMAT` and `GENBAS` in an empty directory, set `ACES_EXE_PATH` to the
directory that holds `xaces3` and its `sio/` subdirectory, and start
`xaces3` under MPI:

```bash
module load intel/2025.1.0 openmpi/5.0.7
export ACES_EXE_PATH=/path/to/ACES_PROGRAM/aces3/bin
mpirun -np 4 $ACES_EXE_PATH/xaces3 > out
```

`xaces3` reads the input with the ACES II `xjoda`, which it runs as
`xjoda aces3` in the work directory at the start of the job and after each
gradient of an optimization or finite-difference run. `xjoda` must be on
`PATH` (for a superbuild, `ACES_PROGRAM/aces2/bin`), or name it with
`ACES_XJODA`:

```bash
export ACES_XJODA=/path/to/ACES_PROGRAM/aces2/bin/xjoda
```

ACES III runs in C1 only. In this mode `xjoda` sets `SYMMETRY=OFF`: the
calculation is done in C1, but the molecule is put in the standard
orientation of its full point group, which is found and stored. With
`SYMMETRY=NONE` or `NOREORI=ON` in the input, the input orientation is kept
instead (and no point group is found). A-, G- and D-tensor jobs
(`PROPS=A-TEN`, `G-TEN`, `D-TEN`) always keep the input orientation and
origin, because the gauge origin of these properties is the coordinate
origin. `xjoda` also uses the ACES III
defaults `CC_MAXCYC=0` (100 iterations) and `EOM_MAXCYC=30` unless they are
given.

### Orbital symmetry labels

At the end of every job `xaces3` runs `xa2proc symadapt`, which gives each
C1 SCF orbital the irreducible representation of the full point group and
appends the table to `summary.out` (for water, 1a<sub>1</sub>,
2a<sub>1</sub>, 1b<sub>2</sub>, 3a<sub>1</sub>, 1b<sub>1</sub>, the same
labels and orbital energies as an ACES II run in C<sub>2v</sub>). When an
EOM or CIS program has written the excited state densities to `JOBARC`, the
natural orbitals of each state are labeled too. `xa2proc` is found like
`xjoda`: `$ACES_XA2PROC`, else next to `$ACES_XJODA`, else on `PATH`. A
`SYMMETRY=NONE` job is not labeled. The labels of fully symmetric orbitals
of single atoms are not assigned yet (they print as `XXXX`).

Keep all processes of a job on **one node** unless the scratch directory
(`LOCAL_PATH`, below) is on a shared file system. Each process keeps its
block files in the run directory, and the processes must see each other's.
Under Slurm, use `#SBATCH --nodes=1`.

## Input

The `*ACES2` namelist works as in ACES II and selects the method. ACES III
adds a `*SIP` section after it with its own parameters, one per line:

```{literalinclude} ../../aces3/tests/fast/cases/ccsd_rhf/ZMAT
:language: text
```

(Test `ccsd_rhf`: RHF CCSD for N<sub>2</sub>H<sub>4</sub>, reference total
energy −111.611012182.)

### Choosing the SIAL programs

Each `SIAL_PROGRAM` line names one `.sio` program; they run in the order
given and pass data through `JOBARC` and the block files. If the `*SIP`
section gives no `SIAL_PROGRAM`, `xaces3` picks a *jobflow*, a standard
sequence, from `CALC`, `REF`, `DROPMO`, the geometry and frequency
keywords, `EXCITE` and the DIP/DEA keywords. It writes the parameters it
used, including the chosen programs, to `ZMAT.AUTO`.

The standard jobflows (from `aces3/src/aces/aces_sial/default_jobflows`):

```{include} aces3_jobflows.md
```

Naming the programs explicitly, as all test-suite inputs do, is the most
reliable choice. The inputs in `aces3/tests/fast/cases/` are working
examples for every capability listed above; for example, CCSD(T) with a UHF
reference (`ccsdt_uhf`):

```{literalinclude} ../../aces3/tests/fast/cases/ccsdt_uhf/ZMAT
:language: text
```

and EOM-CCSD ionization potentials (`ip_eomcc`):

```{literalinclude} ../../aces3/tests/fast/cases/ip_eomcc/ZMAT
:language: text
```

### `*SIP` parameters

`MAXMEM`
: Memory per process in MB (default 900). Set it for any real job; the
  test inputs use 2000–2800.

`SIAL_PROGRAM`
: A `.sio` program to run; repeat for a sequence. See above.

`COMPANY`, `IOCOMPANY`
: Four numbers each, describing the worker processes and the I/O-server
  processes; only the third, the number of processes, is used. By default
  three quarters of the processes are workers and one quarter I/O servers.

`SIP_MX_SEGSIZE`, `SIP_MX_OCC_SEGSIZE`, `SIP_MX_VIRT_SEGSIZE`
: Block (segment) sizes for the AO, occupied and virtual index ranges.
  When not given they are chosen by `AUTO_SEG_ALGORITHM`:
  `SEGMENT_OPTIMIZED` (the default) uses few, large segments;
  `MEMORY_OPTIMIZED` keeps blocks small.

`SIP_SUB_SEGSIZE`, `SIP_SUB_OCC_SEGSIZE`, `SIP_SUB_VIRT_SEGSIZE`, `SIP_SUB_AO_SEGSIZE`
: Sub-segment sizes used by programs that subdivide blocks (default 5).

`LOCAL_PATH`
: Directory for the I/O servers' block files (default: the run directory).

`INTEGRAL_PACKAGE`
: `ERD` (the default and the only one in use).

`FAST_ERD_MEMCALC`
: `YES` estimates the ERD integral memory instead of computing it exactly;
  useful for large Hessian jobs.

`NWORKTHREAD`
: Number of communication buffers per worker (default 20).

`DEBUG_PRINT`, `TIMERS`
: `YES` prints debugging information or detailed timings for each SIAL
  program. Both can make the output very large.

`ITRIP`, `SUB`, `IHESS`, `JHESS`
: Split a large CCSD(T) job (`ITRIP` *first last*, a range of occupied
  orbitals; `SUB`, the occupied range held in memory) or a Hessian
  (`IHESS`/`JHESS` *first last*, ranges of atoms) into pieces that can run
  as separate jobs. The partial results must then be added up by hand.

`REGULARIZER`, `MOLFILE`, `TRACE`, `AAT_DATABASE`
: Specialized options read by particular programs.

`MASTER_IS_WORKER`
: Obsolete; the master is always a worker.

## Restarts

ACES III checkpoints between and within some SIAL programs. If the file
`aces3.restart` exists in the run directory when `xaces3` starts, it
resumes from the last checkpoint.
