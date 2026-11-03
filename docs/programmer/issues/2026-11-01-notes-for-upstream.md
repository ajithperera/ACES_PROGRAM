# Notes for upstream: untracked `.f` files and the ring-CC ω

The push of 2026-11-01 (`6d3d08c9..84f41c7c` on `ajithperera/ACES_PROGRAM`)
brings two things that need someone on the upstream side. First, the
cpp-generated `.f` files of aces2 are no longer in git, so a `git pull`
deletes them, and local edits made to a `.f` must be moved to the `.F`.
Second, the EOM ring-CC doubles correction uses ω of the first root of
an irrep for every root; whether that is intended needs the authors of
the method.

- **Status:** open (needs action or an answer from upstream)
- **Area:** aces2 build (`aces2/.gitignore`, every module directory);
  aces2 `vee` (`form_dbls_correctns_2rcc.F`)
- **Found with:** commits `6809d831`, `24d517c9`, `6ec048c4` (`.f`);
  `cfbea037` and {doc}`2026-10-29-ring-ne-ag-correction` (ω)

## 1. `.f` files are no longer tracked

### What changed

The makefile of every aces2 module runs cpp on each `.F` and compiles the
`.f` it writes (`${F_f} : %.f : %.F` in the module `GNUmakefile`, with
`process_F = ${CPP} ${cppflags} -D__fortran -D__fortran77`). These
generated `.f` files used to be tracked in git, next to their `.F`.

- `6809d831` and `24d517c9` stopped tracking every `.f` that has a `.F`
  next to it (424 + 2860 files).
- `24d517c9` renamed the 3688 real `.f` sources (no `.F` next to them)
  to `.F`. cpp gives back the same text for all of them except a few
  comment lines; the exceptions are listed under *Generated `.f` files*
  in {doc}`../building`.
- `6ec048c4` renamed `fblocal/fblocal.f` to `.F` and changed its second
  `#ifdef _debug_LVL0` to `_DEBUG_LVL0`, so both debug prints still
  appear as before.
- `aces2/.gitignore` ignores `/*/*.f`, so the `.f` files that the build
  writes do not show up in `git status`.

Only the module directories themselves changed. `aces2/configure/`
(no `GNUmakefile`), subdirectories such as `old/`, `save/`, `Keep/` and
`Test/`, and all of aces3 still track their `.f` files.

### Why

The makefile regenerates a `.f` only when its `.F` is newer, and it does
not follow `#include`. A checkout or a `git push` to a server writes the
tracked `.f` with the same time stamp as the `.F`, so after a header such
as `dftfunc.com` changes, the build compiles `.f` files that still hold
the old header. On the ccsd server this mixed old and new common-block
layouts, and every Kohn–Sham Fock matrix came out NaN.

### What to do after pulling

- **No local changes to `.f` files.** Nothing to do. The pull deletes the
  tracked `.f` files, and the next `make` writes them again from the
  `.F`.
- **Uncommitted edits to a `.f`.** `git pull` stops with "Your local
  changes to the following files would be overwritten". Run
  `git diff -- '*.f'` and save the output, then move the edits:
  - For a `.f` that has a `.F` next to it, make the same edit in the
    `.F`. An edit made only to the `.f` would be lost at the next build
    anyway, because the build writes the `.f` again.
  - For a `.f` that was renamed to `.F`, `git stash`, `git pull`,
    `git stash pop` usually carries the edit over to the `.F` by rename
    detection. Check the result with `git diff`.
- **Local commits that touch a `.f`.** The merge with upstream gives
  modify/delete conflicts for `.f` files that have a `.F` next to them.
  Move the change into the `.F` and resolve with `git rm` on the `.f`.
  Renamed files usually merge on their own.
- **New code.** Add new sources as `.F`, and edit the `.F`, never the
  generated `.f`. A `.F` goes through `cpp -traditional -P`, so watch the
  things cpp changes: the words `unix` and `linux` (both defined as 1),
  `/* … */` (removed), and a line ending in `\` (joined to the next
  line). See {doc}`../building`.

## 2. ω in the EOM ring-CC doubles correction

### The question

`vee/form_dbls_correctns_2rcc.F`, inside the loop `Do My_root = 1,
Nroot(Irrepx)`, calls

    Call Form_d2_4rcc(Work,Rcc_eigs,Len_ov,Length,Iuhf,Irrepx,
   +                  Listd0,1)

so the denominator ω − D2 of every root uses ω of the first root of the
irrep. The CIS version (`form_dbls_correctns_2cis.F`) reads each root's
own energy into `Root` and passes that, so there every root gets its own
ω.

### What changes with the root's own ω

With `My_root` in place of the last `1`, the `Ag` and `B1g` corrections
of the Ne D states agree exactly (triplet −4.276406, singlet −4.340012
eV in both), as they must for components of one atomic state. With ω of
the first root they differ, because the first root of `Ag` (20.730 eV)
and the first root of `B1g` (20.944 eV) are different states.

The change also moves the corrections of all the other roots that are
not first in their irrep, by up to 0.5 eV, and most of the ring tests
then fail against the current references. The module docstring of `aces2/tests/jcp2020_eom_rccd`
records that the current numbers match Table I of the paper. The code
was therefore left as it is.

### What is needed

An answer from the authors of the method: is ω of the first root
intended (for example, to match how Table I was computed), or should
each root use its own ω? If it should be the root's own ω, the change
is one argument, and the suite references then need to be generated
again.

Not related to this question: the Ne `Ag` failures themselves were
caused by left and right eigenvectors that were not biorthonormal and
are fixed (`vee/biorth_rcc_vecs.F`, `cfbea037`; ring 49/49 with ifx
2024.0 and ifx 2023.2).
