# EOM ring-CC: Ne `Ag` doubles correction differs between compilers

Four Ne cases of the `jcp2020_eom_rccd` suite fail with ifort 2024.0 and
with ifx 2023.2: the doubles correction of some `Ag` states differs from
the reference by 0.015–1.4 eV, and the two compilers give different wrong
values. All other numbers match. The failing states are exactly the
`Ag` components of degenerate atomic states, which points to an ill-posed
quantity rather than a compiler bug.

- **Status:** open (known upstream as ifort-only; ifx 2023.2 fails too)
- **Area:** aces2 `vee` (EOM ring-CC doubles correction), test suite
  `aces2/tests/jcp2020_eom_rccd`
- **Found with:** local workstation, ifort 2024.0 + Intel MPI 2021.11
  (`./configure --setvars=/opt/intel/oneapi/setvars.sh --fortran=ifort`);
  ccsd server, ifx 2023.2 + Intel MPI 2021.10 (2026-10-29, job 8098);
  upstream's tag `baseline-20261028` lists it for ifort 2023.2: "Known,
  ifort-only (oneAPI 2023.2): jcp2020 Ne Ag EOM correction column"

## Symptom

`./run_tests.sh ring` gives 41/49 with both compilers. The four Ne cases
fail both directly and through `xaces2` (`--xaces2`); everything else in
the ring suites passes. The same rows fail with both compilers, but the
values differ:

| Case | Row | Base value (eV) | ifort 2024.0 | ifx 2023.2 | Reference |
|---|---|---|---|---|---|
| `ne_eomsrccd` | 0, 1 | 21.026280 | −4.283900 | −4.316053 | −4.298957 |
| `ne_eomsfrccd` | 0, 1 | 19.943569 | −4.254007 | −4.251804 | −4.172985 |
| `ne_eomsdrccd` | 0, 1 | 25.590032 | −4.240939 | −3.850035 | −4.500696 |
| `ne_eomsdrccd` | 4 | 72.687813 | −0.735584 | −0.910249 | −0.881774 |
| `ne_eomsfdrccd` | 0, 1 | 25.030008 | −3.178928 | −4.623036 | −4.606569 |
| `ne_eomsfdrccd` | 4 | 71.840892 | −0.924047 | −0.891052 | −0.682103 |

The base values (first column) match to all printed digits. The
"Exc. Eng." column differs by the same amount as the correction.

## How to reproduce

```
./configure --setvars=/opt/intel/oneapi/setvars.sh --fortran=ifort
./build.sh aces2
cd aces2 && python3 tests/jcp2020_eom_rccd/run_tests.py --case ne_eomsrccd
```

The pipeline is `xjoda xvmol xvmol2ja xvprops xvscf xvtran xintprc
xringcc xringcl xvee`; the table is the "summary of the EOM(S)-RCCD
results" at the end of the `xvee` output (written by
`vee/process_uhf_ccd_roots.F`). Each case takes about 2 s.

## What is known about the cause

Ne is an atom, run in D2h. A D state has five components: two in `Ag` and
one each in `B1g`, `B2g` and `B3g`. In every failing row the base value
also appears in `B1g`, `B2g` and `B3g`; the `Ag` rows that pass have no
such partner. So the failures are the `Ag` components of D states, which
are degenerate within `Ag`.

The reference itself shows that the correction is not the same for all
components of one state. For `ne_eomsrccd`, base value 21.026280 eV:

| Irrep | Reference | ifort 2024.0 | ifx 2023.2 |
|---|---|---|---|
| `Ag` (2 rows) | −4.298957 | −4.283900 | −4.316053 |
| `B1g`, `B2g`, `B3g` | −4.335160 | −4.335160 | −4.335160 |

If the correction were a property of the state, all five components would
give the same number. They do not, even in the reference. When two roots
have the same energy, any rotation of their eigenvectors is an equally
valid solution, and a correction computed from each eigenvector
separately depends on which rotation the eigensolver happened to return.
That choice depends on rounding, so it can change with the compiler while
the eigenvalues stay the same. The `B1g`/`B2g`/`B3g` components are alone
in their irreps and match on every compiler.

This is a working explanation, not yet checked in the code. It fits all
six failing rows and none of the passing ones, and it explains why two
compilers give two different wrong values for the same rows. It also
means "ifort-only" is not the right description: any compiler can land on
another rotation, and ifx 2025.1 matching the reference may be luck.

## Possible fixes

1. In the suite, compare only the base value for `Ag` rows whose base
   value also occurs in another irrep (or for the Ne cases only), and say
   why in the reference table.
2. In `vee`, make the correction invariant inside a degenerate block, for
   example by averaging it over the roots of one block, or by fixing the
   rotation of degenerate eigenvectors before the correction is formed.
   This changes published-style output and needs the method authors.

## Open questions

- Is the local ifx 2024.0 build also affected? ifx 2023.2 on the ccsd
  server is (see the table above); upstream reports ifx 2025.1 on
  HiPerGator as 49/49.
- Why do both `Ag` rows of one state print the same correction? With two
  rotated eigenvectors one would expect two different values. The two
  rows may be the same root found twice.
