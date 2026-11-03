# Correlated energies

`CALCLEVEL` (shortened to `CALC`) selects the method. The main families:

Many-body perturbation theory
: `MBPT(2)`, `MBPT(3)`, `SDQ-MBPT(4)`, `MBPT(4)`; `OO-MP2` with orbital optimization

Coupled cluster
: `CCD`, `CCSD`, `CCSD(T)`, `CCSD[T]`, `ACCSD`, `ACCSD(T)`, `CC2`, `CC3`,
  the iterative triples approximations `CCSDT-1` to `CCSDT-4`, `CCSDT`,
  quadruples corrections and approximations up to `CCSDTQ`, and the linear
  methods `LCCD`, `LCCSD`, `LCCSDT`

Configuration interaction and QCI
: `CID`, `CISD`, `QCISD`, `QCISD(T)`, `QCISD(TQ)`

Pair coupled cluster and ring methods
: `PCCD` and its variants, `OO-CCD`, `OO-LCCD` ({doc}`pccd`); `rCCD`, `drCCD` ({doc}`ring_cc`)

The {ref}`keyword-reference` lists every allowed value.

A CCSD energy for water (test `109`, reference −76.2697698466):

```{literalinclude} ../../aces2/tests/fast/cases/109/ZMAT
:language: text
```

`INCORE` keeps integral lists in memory to cut disk I/O, and
`ABCDTYPE=AOBASIS` computes the four-virtual-index contribution from AO
integrals, which avoids storing the largest MO integral list.

## Frozen orbitals and truncated virtual spaces

`DROPMO` removes orbitals from the correlation treatment, counted from the
lowest. `DROPMO=1` freezes the oxygen 1s, and `DROPMO=1-2/10` freezes
orbitals 1, 2 and 10. `TRUNC_ORBS=FNO` replaces the virtual space with
frozen natural orbitals, keeping the percentage given by `FNO_KEEP`
(test `124a.fno`, reference −76.2496011642):

```{literalinclude} ../../aces2/tests/fast/cases/124a.fno/ZMAT
:language: text
```

## Brueckner orbitals

`BRUECKNER=ON` repeats the CC calculation, rotating the orbitals until the
T1 amplitudes vanish to within `BRUCK_CONV`. Test `022` runs CCSD with
`BRUECKNER=ON`, which gives Brueckner CCD (B-CCD).
