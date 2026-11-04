# Optimized effective potential

The optimized effective potential (OEP) method finds the local,
multiplicative potential whose orbitals minimize a given orbital-dependent
energy. ACES II has two OEP methods, both for closed-shell molecules:

`SCF_TYPE=OEPX`
: Exchange-only OEP. The energy is the Hartree–Fock energy expression,
  and the exchange operator in the Fock matrix is replaced by a local
  exchange potential. This is exact-exchange Kohn–Sham theory (Ivanov,
  Hirata, Bartlett, PRL 83, 5455 (1999)).

`SCF_TYPE=OEP2`
: Second-order OEP. The energy adds the second-order many-body
  perturbation theory (MBPT(2)) correlation energy, computed with the
  Kohn–Sham orbitals, and the potential gets a matching local correlation
  part. This is the *ab initio* DFT of the Bartlett group (Bartlett,
  Grabowski, Hirata, Ivanov, JCP 122, 034104 (2005); Schweigert, Lotrich,
  Bartlett, JCP 125, 104108 (2006)).

The code is in `aces2/oeptools` and in the `oep_*` routines of
`aces2/vscf`.

## Running an OEP calculation

The minimum is one keyword in `*ACES2`:

```text
SCF_TYPE=OEPX
```

or `SCF_TYPE=OEP2`. Everything else is optional:

- an `*OEP(...)` block with the OEP options (see [](#oep-block));
- an `*INTGRT` namelist for the numerical grid, as in {doc}`dft`.
  `kspot` and `func` are not used.

When `SCF_TYPE` is `OEPX` or `OEP2`, `xjoda` changes two keywords
without asking:

- `SYMMETRY=OFF`. The calculation runs in C1.
- `SPHERICAL=OFF`. The basis uses Cartesian d and f functions, so
  cc-pVDZ water has 25 functions instead of 24.

Compare OEP energies with Hartree–Fock or MP2 energies computed with
`SPHERICAL=OFF`, or the basis sets differ.

## What the run does

`xaces2` runs three extra programs before `xvscf`:

1. `xintpack` computes the three-index overlap integrals
   ∫ χ<sub>μ</sub> χ<sub>ν</sub> χ<sub>P</sub> on the numerical grid.
   The auxiliary basis {χ<sub>P</sub>}, in which the local potential is
   expanded, is the orbital basis itself. Linearly dependent
   combinations are removed, and the number left is written to `NUMAUX`.
2. `xintpack2vmol` transforms the integrals to the VMOL ordering and
   writes them to the `OEPINT` file.
3. `xvscf` does the SCF iterations.

The grid enters only through these integrals. With the default grid
(`potradpts=50`, `potgrid=4`) the OEPX energy of the water example below
changes by 5×10<sup>-8</sup> hartree on going to `potradpts=100`,
`potgrid=6`, `enegrid=6`.

`xvscf` starts with ordinary Hartree–Fock iterations. When the density
change falls below 10<sup>-OEP1ON</sup>, it switches to the exchange OEP.
For `OEP2`, the correlation potential is added once the density change
falls below 10<sup>-OEP2ON</sup>. The iteration table marks the OEP
iterations:

```text
      10           -76.0269920394              0.1090127648D-02
      10 x         -76.0269946990              0.1090127648D-02
      11 x         -76.0269948684              0.3657760547D-01
```

`x` marks an exchange-only OEP iteration and `xc` an OEP2 iteration. For
`OEP2`, convergence is not declared until an iteration in which the
correlation potential was updated satisfies `SCF_CONV`.

After convergence `xvscf` prints:

- `E (OEP) = ...`, the total energy. It is also stored as `TOTENERG` in
  `JOBARC`.
- The HOMO condition. The local exchange potential is defined up to a
  constant. In a complete basis the HOMO expectation values of the local
  and nonlocal exchange are equal; in a finite basis they are not, and
  the difference is printed. All orbital eigenvalues are then shifted by
  that difference.
- The energy split into nuclear repulsion, one-electron, Coulomb,
  exchange and correlation parts, and into zero-, first- and
  second-order parts.

(oep-block)=
## The `*OEP` block

The options go in a block of the form

```text
*OEP(KEYWORD=VALUE,KEYWORD=VALUE)
```

anywhere in `ZMAT` after `*ACES2`. Pairs are separated by commas or new
lines, and the block ends with `)`:

```text
*OEP(REFERENCE=SEMI,
     PRINTASV=TRUE)
```

The parser is simpler than the one for `*ACES2`:

- Keywords are case sensitive. Write them, and string values, in upper
  case. An unknown or misspelled keyword is ignored without a message.
- Keyword names cannot be abbreviated.
- Do not put `#` comments in the block.

Use `PRINTASV=TRUE` to have the values actually used printed in the
output.

| Keyword | Type | Meaning | Default |
|---|---|---|---|
| `OEP1ON` | integer | Switch from Hartree–Fock to OEP iterations when the density change is below 10<sup>-OEP1ON</sup> | 1 |
| `OEP2ON` | integer | Add the second-order potential when the density change is below 10<sup>-OEP2ON</sup> (`OEP2` only) | 2 |
| `OEP2OFF` | integer | Stop updating the second-order potential when the density change is below 10<sup>-OEP2OFF</sup> | 20 |
| `OEP2ALWAYS` | integer | Update the second-order potential every iteration when the density change is below 10<sup>-OEP2ALWAYS</sup> | 6 |
| `OEP2UPDATE` | integer | Otherwise update it every `OEP2UPDATE` iterations | 1 |
| `REFERENCE` | string | Zeroth-order Hamiltonian for the second-order energy and potential: `OEP` (the OEP eigenvalues, also accepted as `KS`), `DIAG` (diagonal Fock matrix elements, also `DF`), `SEMI` (diagonal Fock elements in semicanonical orbitals, also `SC`) | `OEP` |
| `HYBRID` | integer | Percentage of nonlocal (Hartree–Fock) exchange in the exchange operator; the rest is the local OEP exchange. `HYBRID=100` gives Hartree–Fock | 0 |
| `ASYMPTOTIC` | string | Correction to the long-range behavior of the exchange potential: `NONE`, `FA` (scaled Coulomb, Fermi–Amaldi), `SLATER` (Slater potential, evaluated on the grid) | `NONE` |
| `SVDTHR` | integer | Singular values below 10<sup>-SVDTHR</sup> are dropped when the response function is inverted | 8 |
| `DFF`, `DFW`, `DWW`, `FDF`, `WDW` | boolean | Include the named diagram in the second-order potential (`OEP2` only) | `TRUE` |
| `PRINTASV` | boolean | Print a table of all OEP options | `FALSE` |
| `OEPR12`, `R12ON` | boolean, integer | MBPT(2)-F12 contribution to the second-order potential. Untested | `FALSE`, 2 |

`REFERENCE` has a large effect on `OEP2` (see the table below). With
the default, the OEP eigenvalues, the second-order correlation energy of
water is 60 % larger than the MP2 value. The semicanonical choice,
`REFERENCE=SEMI`, is close to MP2.

## Results for water

Water in cc-pVDZ with Cartesian functions, at the geometry of the
examples below. `HF` and `MP2` are ordinary runs with `SPHERICAL=OFF`;
E<sub>c</sub> is the correlation part of the energy.

| Method | Input | Energy (hartree) | E<sub>c</sub> (hartree) |
|---|---|---|---|
| Hartree–Fock | `SCF_TYPE=HF` | −76.0269950587 | – |
| OEPX | `SCF_TYPE=OEPX` | −76.0256741392 | – |
| OEPX | `*OEP(ASYMPTOTIC=FA)` | −76.0257846774 | – |
| OEPX | `*OEP(ASYMPTOTIC=SLATER)` | −76.0254770668 | – |
| 50 % hybrid | `*OEP(HYBRID=50)` | −76.0266873270 | – |
| 100 % hybrid | `*OEP(HYBRID=100)` | −76.0269950587 | – |
| MP2 | `CALC=MBPT(2)` | −76.2347027494 | −0.2077 |
| OEP2 | `SCF_TYPE=OEP2` | −76.3400419983 | −0.3296 |
| OEP2 | `*OEP(REFERENCE=DIAG)` | −76.2404440111 | −0.2181 |
| OEP2 | `*OEP(REFERENCE=SEMI)` | −76.2352132622 | −0.2124 |

The OEPX energy lies 1.3 millihartree above Hartree–Fock, as it must:
both minimize the same energy expression, and the OEP restricts the
orbitals to those of a local potential. `HYBRID=100` reproduces the
Hartree–Fock energy to every printed digit. The integral program does
not matter: `INTEGRALS=ERD` gives the same OEPX energy as VMOL.

## Examples

All inputs below are ACES II test cases. The energies are the references
the suite checks.

Exchange-only OEP for water (test `oepx_h2o`, −76.0256741392). No
`*OEP` block and no `*INTGRT` namelist, so every option has its default:

```{literalinclude} ../../aces2/tests/fast/cases/oepx_h2o/ZMAT
:language: text
```

OEP2 with the default reference Hamiltonian (test `oep2_h2o`,
−76.3400419983):

```{literalinclude} ../../aces2/tests/fast/cases/oep2_h2o/ZMAT
:language: text
```

OEP2 with the semicanonical reference Hamiltonian (test `oep2sc_h2o`,
−76.2352132622):

```{literalinclude} ../../aces2/tests/fast/cases/oep2sc_h2o/ZMAT
:language: text
```

## Limits

- **Closed shell only.** `REFERENCE=UHF` stops in `xvscf` with
  `The OEP method for an unrestriced reference is not implemented yet`.
- **No symmetry and no spherical functions**, as described above.
- **Energies only.** Gradients and correlated calculations that start
  from OEP orbitals have not been tested.
- **Run in an empty directory.** `xintpack` creates its files with
  `status='new'` and stops if `AUXTRANSF` or the other files from an
  earlier run are still there.
