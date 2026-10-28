# Coding conventions and traps

ACES II is fixed-form Fortran 77 with C preprocessing, implicit typing, and
no interfaces between routines. The compiler catches very little of what
goes wrong, so these rules matter.

## Column 72

Fixed-form source ends at column 72. The compiler ignores everything after
it **without a warning**. If a trailing comma falls in column 73, the next
continuation line joins the previous argument, implicit typing turns the
result into a new variable name, and the call compiles with one argument
missing. Check line lengths after every edit:

```bash
awk 'length($0)>72 && !/^[Cc!*#]/ {print FILENAME": "FNR}' file.F
```

Break long lines with a continuation character in column 6.

## One include-using routine per file

The common-block include files (`machsp.com`, `sym.com`, `syminf.com`,
`sympop.com`, `flags.com`, ...) are protected by `#ifndef` guards. The guard
applies to the whole preprocessed **file**, not to each subroutine. A second
subroutine in the same `.F` file that includes the same file gets nothing:
with implicit typing its `NIRREP`, `IINTFP`, `POP`, ... silently become
uninitialized local variables, and array-valued names become external
functions. Put each subroutine that includes a `.com` file in its own `.F`
file, as the existing code does.

## Implicit typing

The code relies on `Implicit Double Precision (A-H,O-Z)`, with names
starting I–N integer. The build's `-i8` flag makes those integers 8 bytes.
Follow that convention. Do not declare `Integer*8` (or any other explicit
width); either rely on implicit typing or use `Implicit None` and declare
every variable in that routine. Remember that a name like `Dissiz` is real
unless declared otherwise.

## Record lengths

Lengths passed to `GETREC`, `PUTREC` and the list routines are in integer
words. Write `n*IINTFP` for a double-precision array; never a literal `2`.

## Debug output

`#ifdef _DEBUG_LVL0` blocks are part of the development build (see
{doc}`building`). Keep them; production builds switch them off.

## Symmetrize response matrices

Property matrices built from left and right response vectors (for example
polarizabilities or W matrices) are not exactly symmetric numerically. The
established ACES II practice, as in `vcceh/secprp.F`, is to finish with an
explicit symmetrization, `call Symmet2(A, n)` (A(p,q) = (A(p,q)+A(q,p))/2).
A small asymmetry in a new property code usually means that step is
missing.

## When a result is wrong: try a lower optimization level first

The Intel compiler at `-O2` has miscompiled several old ACES routines that
depend on argument aliasing or memory layout the Fortran standard does not
guarantee:

- `joda/mkerdoed.F`: a run of scalar copies between common blocks became an
  aligned vector load on unaligned data (crash).
- `oed/oed__*_prepare_ctr.f`: an overlapping reverse-order array copy
  became `memcpy` (wrong derivative integrals).
- `libr/modmatmul.F`: called with the same array for two arguments by
  `UNITRY` (wrong results).

For an unexplained wrong answer or crash, recompile the suspect files at
`-O1` or `-O0` **before** tracing the algorithm by hand. If the result
changes, keep the file at the lower level (per-file flags in its
`GNUmakefile`) and record why.

## Checking a result

Bug-free code matters more than speed. Before trusting a new method,
check it every way available: against symmetry (results with `SYMMETRY=ON`
and `OFF` must agree), RHF against UHF for closed shells, alternative code
paths, analytic against finite-difference derivatives, and an independent
program where one exists. When there is no independent reference, write a
small, deliberately simple stand-alone code to compute the same quantity.
