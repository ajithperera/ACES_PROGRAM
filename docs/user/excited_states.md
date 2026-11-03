# Excited, ionized and electron-attached states

Excited states of the same charge come from `EXCITE`. `EOMEE` is EOM-CCSD on
a CCSD ground state; `TDA` (CIS), `RPA`, `CIS(D)`, `P-EOMEE`, `STEOM`, and
EOM methods built on ring-CCD ground states ({doc}`ring_cc`) are also available. The number
of states in each irrep is given with `ESTATE_SYM`, one number per irrep,
separated by `/`. `ESTATE_PROP=EXPECTATION` also computes transition moments,
oscillator strengths, and one-electron properties of the excited states.

EOM-CCSD for water, one state in each of the four irreps (test `054a`,
reference ground-state energy −76.2725089278):

```{literalinclude} ../../aces2/tests/fast/cases/054a/ZMAT
:language: text
```

`IP_CALC` and `EA_CALC` compute ionization potentials and electron
affinities of a closed-shell parent state, for example with `IP_EOMCC` and
`EA_EOMCC`. `IP_SYM` and `EA_SYM` give the number of states in each irrep.
`DIP_CALC` and `DEA_CALC` give double ionization and double attachment
(TDA, EOM-CC or STEOM). The electron affinity of the sodium atom with
EA-EOM-CCSD on a UHF reference (test `042`, reference −161.973978215):

```{literalinclude} ../../aces2/tests/fast/cases/042/ZMAT
:language: text
```
