# Properties

`PROPS=FIRST_ORDER` (or `PROPS=1`) computes the dipole and quadrupole
moments, electric field gradients, spin densities and the Cowan–Griffin
relativistic correction from the relaxed density of the chosen method
(test `026`, CCSD, reference −76.2725454575):

```{literalinclude} ../../aces2/tests/fast/cases/026/ZMAT
:language: text
```

Other values of `PROPS`:

- `SECOND_ORDER`: static polarizabilities (CPHF, with `CALC=SCF`)
- `NMR`: NMR chemical shifts with GIAOs (SCF and MBPT(2), with
  `SPHERICAL=OFF`), as in test `037`
- `J_FC`, `J_SD`, `J_SO`, `JSC_ALL`: contributions to NMR spin–spin
  coupling constants (`CALC=CCSD`; the Fermi-contact and spin-orbit terms
  need `REF=UHF`; test `046`)
- `EOM_NLO`: static and frequency-dependent EOM-CCSD polarizabilities
  (tests `s_075b`, `s_075c`)
