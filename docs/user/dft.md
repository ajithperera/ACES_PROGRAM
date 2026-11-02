# Density functional theory

ACES II runs Kohn–Sham (KS) DFT with local, gradient-corrected (GGA),
global hybrid and range-separated hybrid functionals. This includes
CAM-B3LYP and the QTP family from the Quantum Theory Project (CAM-QTP-00,
CAM-QTP-01, CAM-QTP-02, LC-QTP and QTP17). This chapter describes the
input. {doc}`../programmer/range_separated_dft` describes how the
range-separated and QTP hybrids are implemented.

## Running a KS calculation

Kohn–Sham DFT needs two things in `ZMAT`:

- `SCF_TYPE=KS` in `*ACES2`;
- a `*INTGRT` namelist that names the functional.

`REFERENCE` works the same way as for Hartree–Fock: `RHF` (the default)
gives restricted KS, and `UHF` gives unrestricted KS. `SCF_CONV`,
`SCF_MAXCYC` and the extrapolation keywords control the KS iterations as
they do for HF.

The `*INTGRT` namelist comes after `*ACES2` and ends with `*end`. It has
one keyword per line, written as in the test inputs:

```text
*INTGRT
 ksmem=20000000
 kspot=cam-qtp-01
 func=cam-qtp-01
 enegrid=6
 potgrid=6
 potradpts=100
*end
```

Many test inputs also contain a `*SCF` block with `ks` and `ksmem`. No
program reads it: without it the energies are the same to every printed
digit. Give `ksmem` in `*INTGRT`.

## The `*INTGRT` namelist

| Keyword | Meaning | Default |
|---|---|---|
| `kspot` | Functional used for the KS potential during the SCF iterations (see [](#dft-writing-functionals)) | `hf` (no DFT potential) |
| `func` | Functional used for the final energy. INTGRT evaluates it with the converged density | `none` |
| `rs_omega` | Range-separation parameter ω in bohr⁻¹, used by `hf_lr` and `becke_sr` | the ω of the hybrid named in `func`, else in `kspot`, else 0 |
| `ksmem` | Size of the memory block for the numerical integration, allocated in addition to `MEMORY_SIZE` | 1500000 |
| `potradpts` | Radial points around each atom | 50 |
| `enegrid` | Lebedev angular grid at each radial point, for the energy | 4 |
| `potgrid` | Lebedev angular grid at each radial point, for the potential | 4 |
| `radtyp` | Radial quadrature: `Handy`, `Gauss-Legendre`, `Becke`, `TA`, `multiEX` | `Handy` |
| `radscal` | Radial scaling: `none`, `Slater`, `Bragg`, `TA`, `multiEX` | `Slater` |
| `parttyp` | Partitioning of space between the atoms: `rigid`, `fuzzy` | `fuzzy` |
| `partpoly` | Sizes of the partition polyhedra: `equal`, `bsrad` (Bragg–Slater radii), `dynamic` | `bsrad` |
| `fuzzyiter` | Iterations in fuzzy partitioning | 4 |
| `radlimit` | Extent of the radial integration | 3.0 |
| `cutoff` | Densities below 10<sup>-cutoff</sup> are taken as zero | 12 |
| `print_scf`, `print_int` | When to print the SCF summary and the integration parameters: `never`, `dft`, `ks`, `always` | `dft` |
| `print_mos`, `print_occ`, `print_atom`, `print_size` | Same choices, for the MOs, occupations, per-atom integration data and size information | `never` |

The default grid is coarse. The test cases use `potradpts=100`,
`enegrid=6` and `potgrid=6`. With these values the energies differ from
PySCF's large grids by about 10<sup>-7</sup> hartree.

Name the functional in both `kspot` and `func`. `kspot` decides which
potential the orbitals are optimized for, and `func` decides which energy
is reported. Without `func`, the reported energy comes from the `kspot`
functional, as in tests `103.ks` to `105.ks`.

(dft-writing-functionals)=
## Writing a functional

`kspot` and `func` take the same values, in one of three forms:

`name`
: A single term, or one of the named hybrids below. A named hybrid must
  appear on its own.

`exchange,correlation`
: An exchange term and a correlation term, each with coefficient 1, for
  example `kspot=lda,vwn` (SVWN) or `kspot=becke,lyp` (BLYP).

`term,coef,term,coef[,term,coef]...`
: Any linear combination of at least two terms. Each term appears only
  once, for example `kspot=hf_lr,1.0,becke_sr,1.0,lyp,1.0` (LC-BLYP).

### Exchange terms

| Term | Functional |
|---|---|
| `lda` | Slater (LDA) exchange |
| `becke` | Becke 88 |
| `pbe_ex` | Perdew–Burke–Ernzerhof |
| `pw91_ex` | Perdew–Wang 91 |
| `hf` | Exact (Hartree–Fock) exchange |
| `hf_lr` | Long-range exact exchange, from erf(ω r₁₂)/r₁₂ integrals |
| `becke_sr` | Short-range Becke exchange, attenuated as in Iikura, Tsuneda, Yanai and Hirao, JCP 115, 3540 (2001) |

### Correlation terms

| Term | Functional |
|---|---|
| `vwn` | Vosko–Wilk–Nusair, the VWN5 fit (as in Gaussian's B3LYP) |
| `vwn_rpa` | Vosko–Wilk–Nusair fitted to the RPA (libxc `LDA_C_VWN_RPA`) |
| `lyp` | Lee–Yang–Parr |
| `pbe_cor` | Perdew–Burke–Ernzerhof |
| `pw91_cor` | Perdew–Wang 91 |
| `wl` | Wilson–Levy (energy only, no potential) |
| `wi` | Wilson–Ivanov (energy only, no potential, open shell only) |

### Named hybrids

Range-separated hybrids split the electron repulsion as

$$
\frac{1}{r_{12}} =
\underbrace{\frac{\alpha + \beta\,\mathrm{erf}(\omega r_{12})}{r_{12}}}_{\text{exact exchange}}
+ \underbrace{\frac{1 - \alpha - \beta\,\mathrm{erf}(\omega r_{12})}{r_{12}}}_{\text{DFT exchange}},
$$

so that exact exchange goes from α at short range to α + β at long range.

| Name | α | β | ω (bohr⁻¹) | Exchange | Correlation | Reference |
|---|---|---|---|---|---|---|
| `b3lyp` | – | – | – | 0.20 hf + 0.08 lda + 0.72 becke | 0.19 vwn + 0.81 lyp | Becke, JCP 98, 5648 (1993) |
| `cam-b3lyp` | 0.19 | 0.46 | 0.33 | 0.19 hf + 0.46 hf_lr + 0.35 becke + 0.46 becke_sr | 0.19 vwn + 0.81 lyp | Yanai, Tew, Handy, CPL 393, 51 (2004) |
| `cam-qtp-00` | 0.54 | 0.37 | 0.29 | 0.54 hf + 0.37 hf_lr + 0.09 becke + 0.37 becke_sr | 0.20 vwn + 0.80 lyp | Verma, Bartlett, JCP 140, 18A534 (2014) |
| `cam-qtp-01` | 0.23 | 0.77 | 0.31 | 0.23 hf + 0.77 hf_lr + 0.77 becke_sr | 0.20 vwn + 0.80 lyp | Jin, Bartlett, JCP 145, 034107 (2016) |
| `cam-qtp-02` | 0.28 | 0.72 | 0.335 | 0.28 hf + 0.72 hf_lr + 0.72 becke_sr | lyp | Haiduke, Bartlett, JCP 148, 184106 (2018) |
| `lc-qtp` | 0 | 1 | 0.475 | hf_lr + becke_sr | lyp | Haiduke, Bartlett, JCP 148, 184106 (2018) |
| `qtp17` | – | – | – | 0.62 hf + 0.38 lda | 0.20 vwn_rpa + 0.80 lyp | Jin, Bartlett, JCP 149, 064111 (2018) |

The named hybrids set their own ω. Give `rs_omega` only to override it,
or when you build a range-separated functional from `hf_lr` and
`becke_sr` yourself.

## Range-separated hybrids

The long-range exact exchange needs integrals over
erf(ω r₁₂)/r₁₂. When ω ≠ 0, the integral program runs its two-electron
pass a second time and writes these integrals to their own files
(`IIIILR`, `IJIJLR`, `IIJJLR`, `IJKLLR`) next to the ordinary ones. Both
integral programs can do this:

- VMOL, the default;
- `xerdints`, selected with `INTEGRALS=ERD` in `*ACES2`.

The two give the same energies to 5×10<sup>-11</sup> hartree. Point-group
symmetry is used as usual.

Range-separated hybrids have these limits:

- **Energies only.** There are no analytic gradients for `hf_lr` and
  `becke_sr`. Gradients and geometry optimizations stop with
  `no analytic gradients for range-separated hybrids`.
- **VMOL or ERD integrals only.** SEWARD and direct SCF do not write the
  long-range files, and the run stops with an error.
- **One ω per run.** `func` and `kspot` may not name hybrids with
  different ω. For example, `func=cam-b3lyp` with `kspot=lc-qtp` is an
  error.
- **ω is required.** `hf_lr` or `becke_sr` with ω = 0 is an error. A
  negative `rs_omega` is also an error.

## Examples

All inputs below are ACES II test cases. The energies are the references
the suite checks.

An LDA (SVWN) calculation of water in a minimal basis with the default
grid (test `103.ks`, −74.7326254367):

```{literalinclude} ../../aces2/tests/fast/cases/103.ks/ZMAT
:language: text
```

QTP17 for water, a global hybrid (test `qtp17_h2o_rks`,
−76.151936309928):

```{literalinclude} ../../aces2/tests/fast/cases/qtp17_h2o_rks/ZMAT
:language: text
```

CAM-QTP-01 for triplet CH₂ with an unrestricted reference (test
`camqtp01_ch2_uks`, −39.090393038052):

```{literalinclude} ../../aces2/tests/fast/cases/camqtp01_ch2_uks/ZMAT
:language: text
```

LC-QTP for water with the ERD integrals (test `lcqtp_h2o_rks_erd`,
−76.266733480240). The only change from the VMOL input is
`INTEGRALS=ERD`:

```{literalinclude} ../../aces2/tests/fast/cases/lcqtp_h2o_rks_erd/ZMAT
:language: text
```

LC-BLYP built from its terms, with an explicit ω (test `lcblyp_h2o_rks`,
−76.265013867952):

```{literalinclude} ../../aces2/tests/fast/cases/lcblyp_h2o_rks/ZMAT
:language: text
```

## Test cases

The DFT cases in `aces2/tests/fast/cases` check the total energy
(`TOTENERG`). The hybrid cases use cc-pVDZ, the grid above, H₂O (closed
shell, RKS) and triplet CH₂ (UKS). The three LDA cases use the default
grid. An `_erd` case is
the same input with `INTEGRALS=ERD`.

| Functional | H₂O RKS | CH₂ UKS | ERD cases |
|---|---|---|---|
| SVWN, LDA | `103.ks` (H₂O, STO-3G), `104.ks` (C₂H₆, STO-3G), `105.ks` (NH₃, cc-pVTZ) | – | – |
| B3LYP | `b3lyp_h2o_rks` | `b3lyp_ch2_uks` | – |
| CAM-B3LYP | `camb3lyp_h2o_rks` | `camb3lyp_ch2_uks` | both |
| LC-BLYP (terms) | `lcblyp_h2o_rks` | – | – |
| CAM-QTP-00 | `camqtp00_h2o_rks` | `camqtp00_ch2_uks` | H₂O |
| CAM-QTP-01 | `camqtp01_h2o_rks` | `camqtp01_ch2_uks` | H₂O |
| CAM-QTP-02 | `camqtp02_h2o_rks` | `camqtp02_ch2_uks` | H₂O |
| LC-QTP | `lcqtp_h2o_rks` | `lcqtp_ch2_uks` | H₂O |
| QTP17 | `qtp17_h2o_rks` | `qtp17_ch2_uks` | both |

The hybrid energies agree with PySCF to about 2×10<sup>-7</sup>
hartree. The remaining difference comes from the different grids.
