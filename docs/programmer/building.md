# How the build works

{doc}`../getting_started/index` covers building from a fresh clone. This
chapter explains what happens underneath and how to rebuild one piece
while developing.

## The superbuild

`CMakeLists.txt` at the top of the repository compiles nothing itself. Each
part keeps its own build system, and CMake runs it as an `ExternalProject`
through a script in `cmake/`:

| Target | Script | Underlying build |
|---|---|---|
| `aces2` | `cmake/build_aces2.sh` | `gmake all install` with the ACES II makefiles |
| `aces3` | `cmake/build_aces3.sh` | `configure` (autotools) and `gmake` in four dependency levels |
| `pyaces` | `cmake/build_pyaces.sh` | `pyaces/build.sh` (F2PY), after `aces2` |

Every script first sources `cmake/env.sh`, whose `aces_load_modules` does
`module purge; module load $ACES_MODULES`, or nothing when `ACES_MODULES`
is empty. The test targets (`cmake/run_suite_*.sh`, `cmake/rollup.sh`) work
the same way. The older `configure`, `build.sh` and `run_tests.sh` at the
top level are the first, shell-only version of the same orchestration and
still work.

## ACES II makefiles

Every source directory's `GNUmakefile` is a symbolic link to
`aces2/Makefiles/GNUmakefile`, which includes `GNUmakefile.src` (also a
link). The top-level `aces2/makefile` (`GNUmakefile.tl`) recurses through
the directories, libraries first.

**Install locations.** `GNUmakefile.src` installs executables into the
checkout's own `aces2/bin` and libraries into `aces2/lib`, creating them if
needed. A developer who wants executables elsewhere puts an override in the
untracked `aces2/Makefiles/local.mk` (copy `local.mk.example`):

```make
DIR_INSTBIN := $(HOME)/Develop/bin
```

Use `:=` here, not `?=`. `GNUmakefile` sets `DIR_INSTBIN` before it includes
`GNUmakefile.src`, so a conditional assignment would never take effect.
The superbuild passes `DIR_INSTBIN` and `DIR_INSTLIB` on the `gmake`
command line instead.

**Rebuilding one program.** In the program's directory:

```bash
module load intel/2025.1.0 openmpi/5.0.7
cd aces2/vscf
gmake && gmake install
```

If you change a library (`libr`, `acescore`, ...), install it and then relink
every program that uses it. Always run the installed binary: running a stale
copy from another directory is the most common reason a fixed bug seems to
come back. `md5sum $(which xvscf) vscf/xvscf` tells you whether they match.

**Hosts.** `GNUmakefile` holds a table of known machines with their
compilers and flags. On a machine not in the table it falls back to a
generic Linux x86-64 Intel configuration when `ifx` is on `PATH`;
`gmake ACES_HOST=<name>` picks a table entry by name. The Intel build uses
`ifx`/`icx` with `-i8` (8-byte default integers), `-O2`, and MKL.

**Debug output.** `GNUmakefile.src` defines `_DEBUG_LVL0`, which turns on
the debug prints guarded by `#ifdef _DEBUG_LVL0` throughout the code. This
is deliberate during development; production builds turn all debug output
off. Expect large output files from big jobs while it is on.

## ACES III

`cmake/build_aces3.sh` runs `aces3/configure` through a site file
(`aces3/configure_line.ufhpc` by default; set `ACES3_CONFIGURE_LINE` to
use another, such as `configure_line.ttu`). It then builds in four levels,
because later directories need the products of earlier ones:

1. the SIP runtime and the integral, ECP and other libraries, and the SIAL
   compiler's support code;
2. the SIAL compiler, which produces `aces3/bin/sial`;
3. `aces_sial`, which compiles every `.sial` program to `.sio` in
   `aces3/bin/sio`;
4. `src/main`, the `xaces3` executable.

After changing a `.sial` file, rebuilding level 3 is enough; after changing
Fortran or C++ in the runtime, rebuild the library and relink `xaces3`.

ACES III has no input parser of its own: `xaces3` runs the ACES II `xjoda`
(see {doc}`../user/aces3`), so ACES II must be built
first. Keywords that only ACES III uses are in the same table,
`aces2/joda/asv_ctl.c`, and `xaces3` reads them from `IFLAGS2` by slot
number (`aces3/src/main/init_data.F`); move both together.

## pyaces

`pyaces/build.sh` builds the `aces2py` extension module with F2PY,
linking against the libraries in `aces2/lib`. The script carries comments
on each problem it works around; the main ones are:

- Only the source files that define routines declared in the `.pyf`
  interface are given to F2PY. Passing every source makes F2PY wrap every
  subroutine it finds.
- Routines whose source names are not all lower case need explicit symbol
  handling, because the Intel compiler lower-cases external names while
  F2PY assumes the source spelling.
- Several ACES II directories that hold their own copies of library
  routines (`a2mix`, `aux`, `rpa`, `scrnc`) are left out of the link. Their
  copies would otherwise be linked ahead of the library versions and change
  results.

`PYACES_MODULES=""` builds with the current environment instead of loading
modules.
