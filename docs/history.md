# History and citation

ACES (originally *Advanced Concepts and Esoteric Studies*, now *Advanced
Concepts in Electronic Structure*) is the program system of
**Prof. Rodney J. Bartlett** and his group at the Quantum Theory Project,
University of Florida. Professor Bartlett has led its development from the
start. ACES was built to put many-body perturbation theory (MBPT),
coupled-cluster (CC) and equation-of-motion coupled-cluster (EOM-CC)
methods into practical use. Many of the CC and MBPT methods now found in
other quantum chemistry programs first appeared in ACES.

This page summarizes the history given in the program's review article:

> A. Perera, R. J. Bartlett, B. A. Sanders, V. F. Lotrich and J. N. Byrd,
> "Advanced concepts in electronic structure (ACES) software programs",
> *J. Chem. Phys.* **152**, 184105 (2020),
> [doi:10.1063/5.0002581](https://doi.org/10.1063/5.0002581).

Cite that article for ACES as a whole, along with the papers for the
methods you use.

## Origins and ACES I

Work on ACES began in the late 1970s at Battelle Memorial Institute in
Columbus, Ohio, to turn diagrammatic many-body theory into general programs.
That work helped move the field toward size-extensive methods and away
from the variational configuration-interaction approach.

ACES I, in the 1980s, brought the first general-purpose implementations of:

- CCSD, for closed- and open-shell molecules, together with CCSDT-1,
  CCSD[T] and MBPT up to sixth order
- analytic gradients for non-variational many-body methods, through their
  response (relaxed) density matrices
- full CCSDT and CCSDTQ

## ACES II

ACES II was begun in 1990, rewritten for the vector supercomputers of the
time. The first version was written by J. F. Stanton, J. Gauss,
J. D. Watts, W. J. Lauderdale and R. J. Bartlett. It is a set of member
programs that communicate through files. It was the first CC program to use
molecular point-group symmetry (D<sub>2h</sub> and its subgroups) throughout
the correlated calculation, and it has analytic first and second derivatives.

Methods that first appeared in ACES II include:

- full CCSDT for open-shell molecules
- CCSD(T) for non-Hartree–Fock references (ROHF, quasi-RHF, Brueckner and
  Kohn–Sham orbitals)
- EOM-CCSD for excited states (1993), then triples corrections to it,
  analytic EOM-CC gradients, and IP- and EA-EOM-CC for ionized and
  electron-attached states
- similarity-transformed EOM-CC (STEOM-CC, mid-1990s), with the MRCC module
  written to support it, multireference EOM-CC, and extended STEOM for
  doubly excited states, generated with an automatic program generator

Some of the ACES II developers later took part of the code into the
separate program CFOUR.

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

## ACES III and Aces4

ACES III, begun around 2000, carries CC theory to massively parallel
computers. It is written in SIAL, the language of the Super Instruction
Architecture, and has its own integral library. It has been public since
2007 (the source copyrights run from 2003 to 2010). Its user and SIAL
programmer guides list V. Lotrich, M. Ponton, N. Flocke, A. Perera,
E. Deumens and R. J. Bartlett as authors; B. A. Sanders designed much of
the Super Instruction Architecture. The parallel EOM-CCSD gradients are
unique to ACES III, as are its parallel DIP- and DEA-EOM-CCSD, IP/EA-STEOM,
F12 variants, ESR tensors and spin–orbit coupling with EOM and STEOM.
Its reference paper is V. Lotrich, N. Flocke, M. Ponton, A. D. Yau,
A. Perera, E. Deumens and R. J. Bartlett, *J. Chem. Phys.* **128**, 194104
(2008).

Aces4 is a later rewrite of the parallel runtime by V. Lotrich, J. N. Byrd,
N. Jindal, A. Perera, B. A. Sanders, R. J. Bartlett and others. It was one
of the first electronic structure programs to work with block-sparse arrays.
Aces4 is not part of this repository.

ACES II was funded by the U.S. Air Force Office of Scientific Research,
ACES III by the DoD High Performance Computing Modernization Program, and
Aces4 by the Air Force Research Laboratory and ENSCO, Inc. The U.S. Army
Research Office also supported the work.

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
