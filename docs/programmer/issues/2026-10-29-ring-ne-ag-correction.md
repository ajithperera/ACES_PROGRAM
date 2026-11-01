# EOM ring-CC: Ne `Ag` doubles correction differs between compilers

Four Ne cases of the `jcp2020_eom_rccd` suite fail with ifort 2024.0 and
with ifx 2023.2: the doubles correction of some `Ag` states differs from
the reference by 0.015–1.4 eV, and the two compilers give different wrong
values. All other numbers match. The failing states are exactly the
`Ag` components of degenerate atomic states. The cause is that `vee`
used the left and right eigenvectors of the EOM matrix as `DGEEV`
returns them, without making them biorthonormal.

- **Status:** fixed (`vee/biorth_rcc_vecs.F`, 2026-11-01); a second,
  separate question about the denominator is open (see *Open questions*)
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

## Cause

`vee/eom_rcc_diags.F` diagonalizes the EOM matrix of one irrep with
`DGEEV`, sorts the roots and passes the left (L) and right (R)
eigenvectors to `form_dbls_correctns_2rcc.F`. The doubles correction of
root i is formed in `dbls_correctns_2rpa.F` as

    E(D)_i = <L_i| W |R_i>

and is only right when <L_i|R_j> = δ_ij. `DGEEV` does not give that: it
scales every L and every R to unit length on its own.

- For a non-degenerate root, L_i is orthogonal to the other R_j anyway,
  but <L_i|R_i> is not 1. The correction is multiplied by it. The matrix
  is nearly symmetric, so this is small: <L|R> = 0.9993–1.0000 in the
  suite, up to 1.6e-3 eV in N2.
- For a block of degenerate roots, `DGEEV` picks a basis of the left
  block and one of the right block independently. <L_i|R_j> is then not
  diagonal, and <L_i|R_i> depends on rounding. In the local Ne run
  (`ne_eomsrccd`, ifx 2024.0) it was 0.988058 for the two `Ag` roots of
  the triplet D state and 0.991022 for the two of the singlet D state,
  against 1.000000 for every root of the other irreps.

The `Ag` correction printed was the invariant value times <L|R>:
−4.283900 / 0.991022 = −4.322711, which is what the fixed code gives on
any choice of vectors. With a rotation of the same angle on both roots,
<L_1|R_1> = <L_2|R_2>; this is why the two `Ag` rows printed the same
value. The degenerate pairs of `B1g`, `B2g` and `B3g` at 21.055947 eV
(spin label `??????`, not compared by the suite) had the same problem:
their corrections ranged from −3.75 to −4.21 eV and became −4.319 to
−4.321 eV with the fix.

## Fix

`vee/biorth_rcc_vecs.F`, called right after `Sort_rcc_eigs`, goes
through the sorted roots in blocks of equal eigenvalues (to 1e-8 a.u.)
and replaces the left vectors of each block by

    L <- L S^(-T),   S = L^T R,

so that <L_i|R_j> = δ_ij. For a block of one root this divides L_i by
<L_i|R_i>. The right vectors, which `process_*_ccd_roots` prints, are
not changed.

## Verification

- Local ifx 2024.0: the four Ne cases now give the same correction for
  both `Ag` rows of a state, and it no longer depends on which vectors
  `DGEEV` returns.
- The suite references were updated for the rows that changed by more
  than 1e-5 eV (85 rows; the tolerance is 1e-4 eV). Apart from the Ne
  `Ag` rows (0.02–0.24 eV), the largest change is 1.6e-3 eV (N2 `B1u`).
  The base values did not change.
- `./run_tests.sh ring` 49/49, and aces2 and pyaces unchanged.

## Open questions

- **ω in the denominator.** `form_dbls_correctns_2rcc.F` calls
  `Form_d2_4rcc(..., Listd0, 1)`, so the denominator ω − D2 of every root
  uses ω of the first root of the irrep. The CIS caller
  (`form_dbls_correctns_2cis.F`) passes each root's own energy. With the
  root's own ω, the `Ag` and `B1g` corrections of the Ne D states agree
  exactly (triplet −4.276406, singlet −4.340012 eV in both), while with
  ω of the first root they differ, because the first root of `Ag`
  (20.730 eV) and of `B1g` (20.944 eV) are different states. But the
  change moves the corrections of all higher roots by up to 0.5 eV, and
  the module docstring of the suite records that the current numbers
  match Table I of the paper. The code was therefore left as it is; this
  needs the authors of the method.
- Not yet checked with ifx 2023.2 and ifort 2024.0 (the ccsd run will
  show it).
