# How ACES II is organized

## Member executables and the driver

ACES II is not one program but a set of *member executables*, each built
from one source directory: `joda/` gives `xjoda`, `vscf/` gives `xvscf`,
and so on. They never call each other. They run one after another in the
same directory and pass everything through files there, chiefly `JOBARC`
(named records, see {doc}`io`) and the list files `MOINTS` and friends.

`xaces2` (`aces2/aces2/`) is the driver. It runs `xjoda`, reads the flags
that `xjoda` stored in `JOBARC`, and then runs the programs the job needs in
order (`aces2/aces2/job_control.F`). A CCSD single-point energy, for
example, runs:

| Step | Program | Does |
|---|---|---|
| 1 | `xjoda` | read ZMAT, find symmetry, set up `JOBARC` |
| 2 | `xvmol` | AO integrals |
| 3 | `xvmol2ja` | AO/SO transformation and ordering records |
| 4 | `xvscf` | SCF |
| 5 | `xvtran` | AO → MO integral transformation |
| 6 | `xintprc` | sort MO integrals into lists; MBPT(2) energy |
| 7 | `xvcc` | CC amplitude equations and energy |

With `INTEGRALS=ERD` (see the keyword reference), steps 2 and 3 are
`xerdints` alone, and gradients use `xerdgrad` in place of `xvdint`.

Gradients add `xlambda`, `xdens`, `xanti`, `xbcktrn` and `xvdint`;
excited states add `xlambda` and `xvee` (or `xvea`, `xvip`); `PROGRAM=PCCD`
hands the whole job to `xrunpccd`. Because the sequence is just a list of
programs, any step can be run by hand, which is often the quickest way to
debug one module: run `xaces2` once, then rerun only the program in question
in the same directory.

`xjoda` also does everything *between* single points: geometry steps,
finite-difference displacements, and restarts. In an optimization `xaces2`
loops, calling `xjoda` after each gradient until `xjoda` reports
convergence.

## The main programs

| Program | Role |
|---|---|
| `xjoda` | Input parsing (the `*ACES2` namelist, see `joda/asv_ctl.c`), symmetry, geometry optimization, finite differences, restarts |
| `xvmol` | One- and two-electron AO integrals (VMOL, by J. Almlöf and P. R. Taylor, with ECPs added) |
| `xvmol2ja` | Transformation matrices between computational and ZMAT order, Cartesian and spherical functions |
| `xerdints` | `INTEGRALS=ERD`: one- and two-electron SO integrals from Flocke's OED/ERD, in VMOL's file format, with the records `xvmol2ja` writes; builds the SO basis itself |
| `xvprops` | One-electron property integrals |
| `xnddo`, `xmopac` | Semi-empirical starting density; starting Hessian |
| `xvscf` | Hartree–Fock and Kohn–Sham SCF (`xintgrt` for the numerical integration, `xintpack` for OEP) |
| `xdirmp2`, `xscfgrd` | Direct MBPT(2) and SCF gradients with GAMESS integrals |
| `xvtran` | Four-index AO → MO transformation |
| `xintprc` | Sorts MO integrals into the OOOO, OOOV, OOVV, OVVV and VVVV lists; MBPT(2) energy; SCF stability analysis |
| `xvcc`, `xvcc5t`, `xvcc5q` | CC amplitude equations and non-iterative corrections; also finite-order MBPT |
| `xccsdtq`, `xccsdtqf`, `xlcct` | Higher-order CC (CCSDTQ and quadruples/pentuples corrections) |
| `xfno` | Frozen natural orbitals |
| `xlambda` | CC Λ equations |
| `xvee`, `xvea`, `xvip` | EOM-CC excitation energies, electron affinities, ionization potentials; excited-state densities |
| `xvcceh` | EOM-CCSD polarizabilities and NMR spin–spin coupling constants |
| `xdens` | Correlated one- and two-particle densities in the MO basis |
| `xprops` | First-order properties, Mulliken analysis, relativistic corrections |
| `xanti`, `xbcktrn` | Sort the two-particle density and transform it to the AO basis |
| `xvdint`, `xvksdint` | Integral derivatives and gradients (ABACUS, by T. Helgaker, P. Jørgensen, H. J. Aa. Jensen and P. R. Taylor); GIAO integrals |
| `xerdgrad` | `INTEGRALS=ERD`: gradients from OED/ERD derivative integrals over the symmetry-unique shell pairs and quartets, using the densities `xvdint` reads (SCF vectors, `RELDENS`, `IINTERM`, `NHFDENS`, the back-transformed Γ on `MOINTS`) |
| `xcphf` | Coupled-perturbed Hartree–Fock |
| `xnmr` | Correlated NMR shielding (paramagnetic part) |
| `xtdhf`, `xrpa` | Time-dependent HF and RPA excitation energies |
| `xrunpccd`, `xpccd`, `xoorot` | pCCD, OO-pCCD, DOCI |
| `xringcc`, `xringcl` | Ring and direct-ring CCD |
| `xa2proc` | Utilities, see below |

## Libraries

| Directory | Contents |
|---|---|
| `acescore` | The I/O layer: `aces_init`/`aces_fin`, `JOBARC` records (`GETREC`, `PUTREC`), the list files (`GETLST`, `PUTLST`, `aces_list_touch`), and their cache |
| `libr`, `libr2`, `librt3` | The main ACES library: routines that need ACES common blocks or parameters |
| `tools` | Utilities that do not depend on anything in ACES |
| `sbtools`, `ecplib`, `eispack`, `linpack`, `lb` | Support code: namelist reader, ECP integrals, eigensolvers, and reference BLAS/LAPACK (do not edit these; they are netlib copies) |
| `crust` | A wrapper around `xgemm` that guards against zero leading dimensions |
| `include` | Shared headers: `*.h` files hold only preprocessor defines (usable from C and Fortran), `*.com` files hold common blocks |

## Starting a program: memory and common blocks

Every member program starts with `aces_init` (or the older `crapsi`, which
calls it). It opens `JOBARC`, reads the job flags into `IFLAGS`/`IFLAGS2`,
and, when asked to allocate memory, allocates the work heap (`MEMORY_SIZE`)
and starts the list I/O system and its cache. The heap is addressed through
the integer array `ICORE`: `aces_init` returns the index `ICORENDX` of its
first element and its length `ICOREDIM`, and programs carve their arrays
out of it by offsets. Programs end with `aces_fin`.

## The utility program `xa2proc`

`xa2proc help` lists its modules. The most useful when developing:

`jarec` / `jareq` *type* *NAME* *dim*
: Print a `JOBARC` record. Types: `i`, `d` (vectors), `ai`, `ad` (arrays,
  dimension `RxC`). Example: `xa2proc jareq d TOTENERG 1`.

`jasum`, `iosum`
: Summaries of the `JOBARC` records and of the list files.

`zerorec` *NAME* ...
: Overwrite records with zeros.

`rmfiles`
: Delete the list files and reset their bookkeeping so the next program
  can start the I/O system cleanly.

`clrdirty`
: Clear the "dirty" flag on the file set. In an optimization or
  frequency job with `RESTART=ON`, the file set is flagged dirty while a
  step runs; `xaces2` clears it before calling `xjoda`. If you run the
  programs by hand, clear it yourself, or `xjoda` will treat the step as a
  crash, restore the previous files, and repeat the same geometry forever.

`mem` *amount*
: Change the memory size stored in an existing `JOBARC` (for example
  `xa2proc mem 2GB`) until the next `xjoda` run.

`molden`, `hyperchem`, `xyz`
: Write the geometry, orbitals and frequencies for viewers, or print the
  current geometry.

`symadapt`
: ACES III: label the C1 SCF orbitals (and excited state natural
  orbitals) of an `xaces3` run with the irreps of the full point group,
  appending the table to `summary.out`. `xaces3` runs it at the end of each
  job. The code is in `a2proc/a3_symadapt_*.F`; routines that have other
  versions elsewhere in ACES II carry the prefix `A3SA_`.

`parfd`, `extrap`, `factor`, `statthermo`, `test`
: Finite-difference bookkeeping, energy/gradient extrapolation, and helpers
  used by scripts and the old test harness.
