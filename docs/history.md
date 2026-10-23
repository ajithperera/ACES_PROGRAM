# History and citation

## ACES II

ACES II was begun in 1990 at the Quantum Theory Project, University of
Florida. The first version was written by J. F. Stanton, J. Gauss,
J. D. Watts, W. J. Lauderdale and R. J. Bartlett. It is a set of member
programs that communicate through files, built around many-body
perturbation theory and coupled-cluster methods, with full use of
D<sub>2h</sub> symmetry and analytic gradients and second derivatives.

The citation given with release 2.5.0 (2006) is:

> ACES II is a program product of the Quantum Theory Project, University of
> Florida. Authors: J. F. Stanton, J. Gauss, S. A. Perera, J. D. Watts,
> A. D. Yau, M. Nooijen, N. Oliphant, P. G. Szalay, W. J. Lauderdale,
> S. R. Gwaltney, S. Beck, A. Balková, D. E. Bernholdt, K. K. Baeck,
> P. Rozyczko, H. Sekino, C. Huber, J. Pittner, W. Cencek, D. Taylor, and
> R. J. Bartlett. Integral packages included are VMOL (J. Almlöf and
> P. R. Taylor); VPROPS (P. Taylor); ABACUS (T. Helgaker, H. J. Aa. Jensen,
> P. Jørgensen, J. Olsen, and P. R. Taylor); HONDO/GAMESS (M. W. Schmidt,
> K. K. Baldridge, J. A. Boatz, S. T. Elbert, M. S. Gordon, J. J. Jensen,
> S. Koseki, N. Matsunaga, K. A. Nguyen, S. Su, T. L. Windus, M. Dupuis,
> J. A. Montgomery).

Basis sets from the PNNL Basis Set Exchange (formerly the EMSL Basis Set
Library) should also be cited as that library asks.

Contributions named in the 2006 manual:

- Two-determinant CCSD energies: P. G. Szalay and A. Balková; analytic
  derivatives: P. G. Szalay
- Dropped orbitals in analytic derivatives (RHF, UHF, ROHF): K.-K. Baeck
- EOM-CCSD dynamic polarizabilities: J. F. Stanton, S. A. Perera and
  M. Nooijen
- EOM-CCSD NMR spin–spin coupling constants: S. A. Perera and M. Nooijen
- Partitioned EOM-CCSD excitation energies: S. R. Gwaltney and M. Nooijen
- EOM-CCSD excited-state gradients: J. F. Stanton and J. Gauss

```{note}
Methods added since 2006 (among them pCCD and orbital-optimized pCCD,
DOCI, ring and direct-ring CCD and their EOM extensions, and RPA) are not
yet credited here.
```

## ACES III

ACES III, the parallel program built on the Super Instruction Architecture,
was developed at the Quantum Theory Project in the 2000s (the source
copyrights run from 2003 to 2010). Its user
and SIAL programmer guides list V. Lotrich, M. Ponton, N. Flocke,
A. Perera, E. Deumens and R. J. Bartlett as authors.

## This repository

The three codes were brought together in one repository, `ACES_PROGRAM`,
in 2026, with the history of each kept through `git subtree` from its
former repository (`ajithperera/ACES-II`, `ajithperera/acesiii`,
`ajithperera/acespy`). The work added a common CMake build, automated
regression suites for all three parts, and this manual.

## Earlier manuals

This manual replaces:

- *ACES II Release 2.5.0 User Manual*, draft, Quantum Theory Project,
  29 January 2006
- *ACES II Programmer's Manual* (1990s)
- *ACES III User Guide*, *SIAL Programmer Guide* and *SIP Design*, in
  `aces3/doc/`

Their content has been carried over where it still describes the code, and
checked against the code where possible.
