# The ACES program system

ACES is a quantum chemistry program system from the Quantum Theory Project at
the University of Florida, built around many-body perturbation theory and
coupled-cluster methods. This repository holds its three parts:

ACES II (`aces2/`)
: The serial Fortran code: SCF, MBPT, CC and EOM-CC energies, analytic
  gradients and second derivatives, properties, pCCD, and RPA/ring-CC
  methods. It is also the correctness reference for ACES III.

ACES III (`aces3/`)
: The parallel code. Methods are written in SIAL, the Super Instruction
  Assembly Language, and run by the SIP runtime under MPI.

pyaces (`pyaces/`)
: A Python interface (F2PY) to the ACES II libraries.

```{toctree}
:maxdepth: 2
:caption: Getting started

getting_started/index
```

```{toctree}
:maxdepth: 2
:caption: User guide

user/aces2_input
user/methods
user/aces3
user/pyaces
```

```{toctree}
:maxdepth: 1
:caption: Reference

reference/keywords
```

```{toctree}
:maxdepth: 2
:caption: Programmer manual

programmer/index
```

```{toctree}
:maxdepth: 1
:caption: About

history
```
