# Methods and calculation types

A handful of keywords decide what an ACES II job does. Each row links to
the chapter that describes it:

| What | Keyword | Examples | Chapter |
|---|---|---|---|
| Reference wavefunction | `REFERENCE`, `SCF_CONV`, `HFSTABILITY` | `REF=UHF`, `HFSTABILITY=FOLLOW` | {doc}`scf` |
| Kohn–Sham DFT | `SCF_TYPE=KS`; `kspot`, `func` in `*INTGRT` | `func=cam-qtp-01`, `func=b3lyp` | {doc}`dft` |
| Correlation method | `CALCLEVEL` | `MBPT(2)`, `CCSD`, `CCSD(T)` | {doc}`correlated` |
| Excited, ionized, attached states | `EXCITE`, `IP_CALC`, `EA_CALC`, `DIP_CALC`, `DEA_CALC` | `EXCITE=EOMEE` | {doc}`excited_states` |
| Pair coupled cluster | `CALCLEVEL`, `PCCD_REF`, `OO_METHOD` | `CALC=PCCD`, `PCCD_REF=OPTIMIZED` | {doc}`pccd` |
| Ring coupled cluster and RPA | `CALCLEVEL`, `EXCITE`, `EOMREF` | `CALC=RCCD`, `EXCITE=DRPA` | {doc}`ring_cc` |
| Properties | `PROPS`, `ESTATE_PROP` | `PROPS=FIRST_ORDER` | {doc}`properties` |
| Geometry | `*` on Z-matrix coordinates, `GEOM_OPT`, `OPT_METHOD` | `OPT_METHOD=EVFTS` | {doc}`geometry_optimization` |
| Frequencies | `VIBRATION` | `VIB=EXACT`, `VIB=FINDIF` | {doc}`frequencies` |

Every example in these chapters is an input from the ACES II test suite,
included directly from `aces2/tests/`, with the reference energy the suite
checks against. They are small enough to run in seconds to minutes.

```{toctree}
:maxdepth: 1

scf
dft
correlated
excited_states
pccd
ring_cc
properties
geometry_optimization
frequencies
```
