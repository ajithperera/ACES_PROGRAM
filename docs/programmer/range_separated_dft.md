# Range-separated and QTP hybrid functionals

This chapter describes how ACES II computes range-separated hybrid
functionals (CAM-B3LYP, CAM-QTP-00/01/02, LC-QTP) and the global hybrid
QTP17. It covers the working equations, the program flow with either
integral program (VMOL, the default, or ERD with `INTEGRALS=ERD`), every
routine that was added or changed, the input rules, the limits, and how
the code was checked. {doc}`../user/methods` shows how to run these methods.

## Summary

### What ACES did not have

- No erf(ω r₁₂)/r₁₂ integrals: exact exchange could only be mixed in by a
  constant factor (`coef_pot_nonlocal`).
- No erfc-attenuated DFT exchange (the ITYH `becke_sr`).
- Only the VWN5 fit of VWN correlation (`vwn`), not the RPA form
  (`vwn_rpa`).
- B3LYP was the only hybrid, and its coefficients were written into
  several routines.
- ERD (`aces2/erd`, used by `xerdints` with `INTEGRALS=ERD`) computed only
  $1/r_{12}$ integrals; the erf, erfc and attenuated operators appear
  nowhere in it. The `erd_F12` copy in ACES III has F12, K12 and Yukawa
  operators, but not erf, and ACES II does not use it.

### What was added

- **VMOL**
  - `rsattn.F` turns the ordinary integrals into long-range ones by
    changing ρ and the prefactor.
  - `rsfiles.F` opens the long-range integral files `IIIILR`, `IJIJLR`,
    `IIJJLR` and `IJKLLR`.
- **ERD and `xerdints`**
  - `erd/erd__rs_attenuate.F` turns the Rys quadrature of $1/r_{12}$ into
    that of erf(ω r₁₂)/r₁₂.
  - `erd/erd__rsomega.inc` holds ω for ERD.
  - `erdints/erd_driver.F` writes the same long-range files as VMOL in a
    second pass.
- **ω from the input:** `intgrt/rsomega.F`, shared by VMOL and `xerdints`.
- **VSCF:** `lrxfock.F` adds the long-range exchange to the KS Fock matrix.
- **INTGRT**
  - `lrexch.F` computes the long-range exchange energy.
  - `ityh_attn.F`, `becke_sr_spin.F` and `func/pot_exch_becke_sr.F`
    compute the short-range Becke exchange.
  - `func/pot_corr_vwnrpa.F` computes VWN_RPA.
  - `hybparm.F` is the single table of the coefficients and ω of every
    hybrid.
- **New functionals:** `cam-b3lyp`, `cam-qtp-00`, `cam-qtp-01`,
  `cam-qtp-02`, `lc-qtp` and `qtp17`.
- **New terms and keyword:** `hf_lr`, `becke_sr`, `vwn_rpa` and `rs_omega`.

### What was changed

- **`vmol/molecu.F`, `conlor.F`, `conlox.F`:** when ω≠0, the integral
  pass (TWOEL) runs a second time to write the long-range integrals.
- **`libr/mkrhff.F`, `mkuhff.F`:** call `LRXFOCK`; stop with an error
  when combined with SEWARD or direct SCF.
- **`intgrt/initintgrt.F`:** reads the hybrids from the `hybparm` table,
  determines ω and checks the input.
- **`integ*.F`, `numint*.F`, `anlytc.F`**
  - compute the potential and energy of the new terms;
  - get the hybrid energies from a generic loop over the `hybparm` table
    instead of the formula written for B3LYP.
- **Include files**
  - `dftfunc.com`: more functional indices.
  - `intgrtflags.com`: `coef_pot_lrx` and `int_rs_omega`.
  - `energ.com`: `ehflx`.
- **`integ*.F`, `vksdint/integAG.F`:** zero every functional slot of the
  grid potential before it is summed. Unused slots used to hold stack
  garbage; with ifx 2023.2 that was NaN, and 0×NaN made every KS Fock
  matrix NaN.
- **`vksdint`:** stops with an error when gradients are requested for a
  range-separated hybrid.
- **ERD:** `erd__e0f0_pcgto_block.f` and `erd__ssss_pcgto_block.f`
  attenuate when ω≠0; `erd__gener_eri_batch.f` and
  `erd__memory_eri_batch.f` then send all s/p quartets except (ss|ss) to
  the general Rys path.
- **`GNUmakefile`:** `xvmol` links `intgrt`; `xerdints` links `intgrt` and
  `sbtools`.
- **Tests and docs:** thirteen regression tests, `methods.md` and this
  chapter.

### Remaining limits

- Energies only; there are no gradients yet.
- VMOL or ERD integrals only; SEWARD and direct SCF are not supported.
- One ω per run.

The sections below give the details.

## Theory

### Splitting the exchange operator

A range-separated hybrid uses different amounts of exact exchange at
short and long range. The Coulomb-attenuating method (CAM) of Yanai, Tew
and Handy splits the exchange interaction as

$$
\frac{1}{r_{12}} =
\underbrace{\frac{\alpha + \beta\,\mathrm{erf}(\omega r_{12})}{r_{12}}}_{\text{exact exchange}}
+
\underbrace{\frac{1 - \alpha - \beta\,\mathrm{erf}(\omega r_{12})}{r_{12}}}_{\text{DFT exchange}} .
$$

At $r_{12}\to 0$ the exact-exchange fraction is $\alpha$, and at
$r_{12}\to\infty$ it is $\alpha+\beta$. The DFT part is written with
$\mathrm{erf} = 1 - \mathrm{erfc}$:

$$
\frac{1-\alpha-\beta\,\mathrm{erf}(\omega r_{12})}{r_{12}}
= (1-\alpha-\beta)\frac{1}{r_{12}} + \beta\,\frac{\mathrm{erfc}(\omega r_{12})}{r_{12}} .
$$

ACES therefore builds every range-separated hybrid from four exchange
terms, each of which has its own name in `*INTGRT`:

| Term | Name | Operator | Where it is computed |
|---|---|---|---|
| $\alpha$ | `hf` | $1/r_{12}$, exact exchange | the usual SCF integrals |
| $\beta$ | `hf_lr` | $\mathrm{erf}(\omega r_{12})/r_{12}$, exact exchange | new long-range integrals from VMOL or ERD |
| $1-\alpha-\beta$ | `becke` | $1/r_{12}$, Becke 88 | numerical grid |
| $\beta$ | `becke_sr` | $\mathrm{erfc}(\omega r_{12})/r_{12}$, Becke 88 | numerical grid, attenuated |

Correlation is not split.

### Long-range two-electron integrals

For primitive Gaussians with exponents $p$ (bra pair) and $q$ (ket pair),
the ordinary Coulomb integral depends on the reduced exponent
$\rho = pq/(p+q)$ through the prefactor and the Boys function
$F_n(\rho\,|PQ|^2)$. The integral over $\mathrm{erf}(\omega r_{12})/r_{12}$
has exactly the same form after two replacements:

$$
\rho \;\to\; \rho' = \rho\,\frac{\omega^2}{\rho+\omega^2},
\qquad
\text{prefactor} \;\to\; \text{prefactor}\times\sqrt{\frac{\omega^2}{\rho+\omega^2}} .
$$

$\rho'$ is used everywhere $\rho$ appears: in the Boys function argument and
in the $(-2\rho)^n$ factors of its derivatives, which build the higher
angular momenta. Nothing else in the integral code changes. As
$\omega\to\infty$, $\rho'\to\rho$ and the ordinary integrals are recovered.

ERD uses Rys quadrature instead. There the Coulomb integral is an integral
over $t\in[0,1]$ of $e^{-Tt^2}$ times a polynomial in $t^2$, with
$T=\rho|PQ|^2$. For $\mathrm{erf}(\omega r_{12})/r_{12}$ the same integral
runs only to $t=\sqrt f$, $f=\omega^2/(\rho+\omega^2)$; the substitution
$t^2 = f s^2$ brings it back to $[0,1]$:

$$
\int_0^{\sqrt f} e^{-Tt^2}\,G(t^2)\,dt
= \sqrt f\int_0^1 e^{-fT s^2}\,G(f s^2)\,ds .
$$

So the Rys roots and weights are taken at $fT$, every root $t^2$ is
multiplied by $f$, and the prefactor by $\sqrt f$. The formulas that build
the 2D integrals from the roots stay as they are. This is the same
attenuation as in VMOL, written for a different integral scheme.

### Short-range Becke exchange (ITYH)

`becke_sr` uses the attenuation of Iikura, Tsuneda, Yanai and Hirao
(JCP 115, 3540 (2001)). For one spin density $\rho_\sigma$ the Becke energy
density is

$$
e_\mathrm{B} = -\rho_\sigma^{4/3} K(x),\qquad
K(x) = C_x + \frac{\beta x^2}{1 + 6\beta x\,\sinh^{-1}x},\qquad
x = \frac{|\nabla\rho_\sigma|}{\rho_\sigma^{4/3}},
$$

with $\beta=0.0042$. ITYH treats $e_\mathrm{B}$ as the exchange of a
uniform gas with the local Fermi wavevector

$$
k_\sigma = \sqrt{\frac{9\pi}{2K(x)}}\;\rho_\sigma^{1/3},
$$

which is the true Fermi wavevector when $K = C_x$ (LDA). The erfc
attenuation of a uniform gas multiplies its exchange by $F(a)$ with
$a = \omega/(2k_\sigma)$:

$$
F(a) = 1 - \frac{8}{3}a\left[\sqrt{\pi}\,\mathrm{erf}\!\left(\frac{1}{2a}\right)
+ (2a - 4a^3)\,e^{-1/(4a^2)} - 3a + 4a^3\right],
$$

$$
\frac{dF}{da} = -\frac{8}{3}\left[B(a) + a\left(12a^2\left(1-e^{-1/(4a^2)}\right) - 3\right)\right],
$$

where $B(a)$ is the bracket in $F$. The energy density is
$e = e_\mathrm{B}\,F(a)$.

For large $a$ (low density, or large $\omega$) the bracket is a difference
of nearly equal large numbers, so for $a>2$ $F$ is summed from its series
in $u = 1/a^2$ instead:

$$
F = \frac{u}{36} - \frac{u^2}{960} + \frac{u^3}{26880} - \frac{u^4}{829440}
+ \frac{u^5}{28385280} - \frac{u^6}{1073479680}
+ \frac{u^7}{44590694400} - \frac{u^8}{2021444812800} .
$$

At $a=2$ the eight-term series is exact to a relative $6\times10^{-17}$,
while the closed form evaluated in double precision has already lost
about $3\times10^{-12}$ to cancellation.

The KS potential needs $\partial e/\partial\rho_\sigma$ and
$\partial e/\partial|\nabla\rho_\sigma|$. Both are analytic: $a$
depends on $\rho_\sigma$ directly (through $\rho_\sigma^{1/3}$) and through
$K(x)$:

$$
\frac{\partial a}{\partial\rho_\sigma} = a\left(\frac{K'}{2K}\frac{\partial x}{\partial\rho_\sigma} - \frac{1}{3\rho_\sigma}\right),
\qquad
\frac{\partial a}{\partial|\nabla\rho_\sigma|} = a\,\frac{K'}{2K}\,\frac{1}{\rho_\sigma^{4/3}} .
$$

### VWN correlation, RPA form

QTP17 uses the VWN correlation fitted to the RPA, not the VWN5 fit that
ACES calls `vwn`. The new `vwn_rpa` uses the VWN form [4.4] with the RPA
parameters (energies in hartree):

| | $A$ | $x_0$ | $b$ | $c$ |
|---|---|---|---|---|
| paramagnetic | 0.0310907 | −0.409286 | 13.0720 | 42.7198 |
| ferromagnetic | 0.01554535 | −0.743294 | 20.1231 | 101.578 |

The spin interpolation is the plain one,

$$
E_c(r_s,\zeta) = e_P + f(\zeta)\,[e_F - e_P],\qquad
f(\zeta) = \frac{(1+\zeta)^{4/3} + (1-\zeta)^{4/3} - 2}{2^{4/3}-2},
$$

without the spin-stiffness term of VWN5. This is the functional that libxc
calls `LDA_C_VWN_RPA`.

## The hybrids

All named hybrids are defined in a single table, `aces2/intgrt/hybparm.F`.
For a range-separated hybrid the table gives $\alpha$, $\beta$, $\omega$
and the correlation, and `hybparm` expands them into the four exchange
terms above.

| Name | $\alpha$ | $\beta$ | $\omega$ (bohr⁻¹) | Exchange | Correlation | Reference |
|---|---|---|---|---|---|---|
| `b3lyp` | – | – | – | 0.20 hf + 0.08 lda + 0.72 becke | 0.19 vwn + 0.81 lyp | Becke III LYP |
| `cam-b3lyp` | 0.19 | 0.46 | 0.33 | 0.19 hf + 0.46 hf_lr + 0.35 becke + 0.46 becke_sr | 0.19 vwn + 0.81 lyp | Yanai, Tew, Handy, CPL 393, 51 (2004) |
| `cam-qtp-00` | 0.54 | 0.37 | 0.29 | 0.54 hf + 0.37 hf_lr + 0.09 becke + 0.37 becke_sr | 0.20 vwn + 0.80 lyp | Verma, Bartlett, JCP 140, 18A534 (2014) |
| `cam-qtp-01` | 0.23 | 0.77 | 0.31 | 0.23 hf + 0.77 hf_lr + 0.77 becke_sr | 0.20 vwn + 0.80 lyp | Jin, Bartlett, JCP 145, 034107 (2016) |
| `cam-qtp-02` | 0.28 | 0.72 | 0.335 | 0.28 hf + 0.72 hf_lr + 0.72 becke_sr | lyp | Haiduke, Bartlett, JCP 148, 184106 (2018) |
| `lc-qtp` | 0 | 1 | 0.475 | hf_lr + becke_sr | lyp | Haiduke, Bartlett, JCP 148, 184106 (2018) |
| `qtp17` | – | – | – | 0.62 hf + 0.38 lda | 0.20 vwn_rpa + 0.80 lyp | Jin, Bartlett, JCP 149, 064111 (2018) |

`vwn` is VWN5. The choice of VWN in each QTP functional (VWN5 for
CAM-QTP-00/01, VWN_RPA for QTP17) follows the definitions in libxc. It was
found by fitting the libxc functionals, at spin-polarized density points,
to linear combinations of their components; the fits reproduce libxc to
about $10^{-16}$. If the original papers or the original ACES II
implementation used another VWN, change the one line in `hybparm.F` and
regenerate the test references.

To add a hybrid, add its index to `dftfunc.com` (raise `fun_num_hyb`) and a
branch to `hybparm.F`. `initintgrt`, the energy printout and VMOL all read
the table, so nothing else needs to change.

## Program flow

```text
xjoda
xvmol      TWOEL         -> IIII, IJIJ, IIJJ, IJKL     (1/r12)
           RSOMEGA       omega from *INTGRT
           RSFILES + TWOEL (RSOMEG=omega)
                         -> IIIILR, IJIJLR, IIJJLR, IJKLLR   (erf(omega r12)/r12)
 or, with INTEGRALS=ERD:
xerdints   ERD_SO2E      -> IIII, IIJJ, IJIJ, IJKL     (1/r12)
           RSOMEGA       omega from *INTGRT
           ERD_SO2E (/ERDRSO/ RSOMEG=omega)
                         -> IIIILR, IIJJLR, IJIJLR, IJKLLR   (erf(omega r12)/r12)
xvscf      each KS iteration:
             MKRHFF/MKUHFF   J + alpha K                (IIII, ...)
               LRXFOCK       + beta K_lr                (IIIILR, IJIJLR)
             INTEGXC -> NUMINTEFF -> INTEG_EFF1
                             V_xc of becke, becke_sr, lyp, vwn, ... on the grid
xintgrt    ANLYTC -> EXCHANGE  E_x(hf)
                  -> LREXCH    E_x(hf_lr)               (IIIILR, IJIJLR)
           NUMINTINT -> INTEGINT
                         E_x(becke), E_x(becke_sr), E_c(...) on the grid
           E = E_Hartree + sum_i cx(i) E_x(i) + sum_j cc(j) E_c(j)
```

### VMOL: the long-range integrals

`vmol/molecu.F` runs `TWOEL` as before, with `RSOMEG=0`. It then calls
`RSOMEGA` to find $\omega$. If $\omega\neq 0$ it prints

```text
  @MOLECU-I, Long-range integrals with omega =  0.330000 bohr**-1.
```

and does the following:

1. `RSFILES` closes the integral files that `READIN` opened (units 2, 81,
   82 and 83) and opens `IIIILR`, `IJIJLR`, `IIJJLR` and `IJKLLR` on the
   same units. As in `READIN`, the number of files depends on `NSYMOP`:
   one file without symmetry, three with one generator, four otherwise.
2. `RSOMEG` (common `/RSINT/`) is set to $\omega$, and `TWOEL` runs a
   second time.
3. In `conlor.F` and `conlox.F`, right after `CALL CMPR1(TMUL,E1)` (where
   `E1` holds $\rho$ and `ROOT` the prefactor of each quartet),
   `RSATTN(E1,ROOT,NCQRS,RSOMEG)` applies the two replacements above.
4. `RSOMEG` is reset to 0.

The long-range files have the same format as the ordinary ones (a
`TWOELSUP` label followed by records of packed indices and values), so the
existing readers work on them unchanged.

`intgrt/rsomega.F` (in the `intgrt` library, so that VMOL and `xerdints`
share it) reads `*INTGRT` with `nl_init('INTGRT',err,.false.)` and
does not call `nl_term`, because `nl_term` would reject every `*INTGRT`
keyword that the integral program does not read. The default $\omega$ is taken from the
first entry of `func`, else of `kspot`, if it names a hybrid in
`hybparm`; `rs_omega` overrides it. A negative $\omega$ is an error.
`initintgrt` applies the same rule, so the integral program, VSCF and
INTGRT always agree on $\omega$.

The link line for `xvmol` (`aces2/Makefiles/GNUmakefile`) now includes
`intgrt`, because `RSOMEGA` calls `HYBPARM`.

### ERD: the long-range integrals with `INTEGRALS=ERD`

With `INTEGRALS=ERD` the driver runs `xerdints` instead of `xvmol` (and
skips `xvmol2ja`). KS calculations already worked on this path; B3LYP
gives the VMOL energy to $10^{-12}$. For range-separated hybrids:

1. `erdints/erd_driver.F` writes the ordinary integral files as before
   (`Erd_so2e`). It then calls `RSOMEGA`. If $\omega\neq0$ it prints

   ```text
    @Erd_driver: long-range integrals  0.330000 bohr**-1
   ```

   and opens `IIIILR`, `IIJJLR`, `IJIJLR` and `IJKLLR` on the same units,
   for the same symmetry cases as the ordinary files. It writes the
   `********` and `TWOELSUP` labels, sets `RSOMEG` in the ERD common block
   `/ERDRSO/` (`erd/erd__rsomega.inc`) to $\omega$, runs `Erd_so2e` again
   and resets `RSOMEG` to 0. The SO transformation, the symmetry
   selection, the thresholds and the file layout are the same as for
   $1/r_{12}$, so `LRXFOCK` and `LREXCH` read these files unchanged.
2. In the general Rys path, `erd__e0f0_pcgto_block.f` calls
   `ERD__RS_ATTENUATE` twice when `RSOMEG≠0`:
   - after it has set `TVAL`, `PQPINV` and `SCALEPQ` for every exponent
     quadruple, to replace $T$ by $fT$ and the prefactor by
     $\sqrt f$ times it;
   - after `ERD__RYS_ROOTS_WEIGHTS`, to multiply every root (ERD stores
     $t^2$) by $f$.

   $\rho$ is recomputed there as $p\,q\,/(p+q)$ from `P`, `Q` and
   `PQPINV`.
3. ERD has dedicated closed-form routines for quartets of s and p shells
   (`erd__ssss_pcgto_block.f` ... `erd__pppp_pcgto_block.f`, called from
   `ERD__1111_CSGTO`), with $1/r_{12}$ built into their formulas. Only
   (ss|ss) is attenuated in place: in `erd__ssss_pcgto_block.f`
   $F_0(T)$ becomes $\sqrt f\,F_0(fT)$. When `RSOMEG≠0`,
   `erd__gener_eri_batch.f` and `erd__memory_eri_batch.f` send every other
   s/p quartet to the general path (`ERD__CSGTO`), which handles them
   correctly; (ss|ss) cannot go there, because the general path assumes
   that the bra or the ket has angular momentum.

`RSOMEG` is zero (the default value of a common block) in every other
program that links ERD (`xerdgrad`, `libint`, `a3ecp`, `dkh`), so they are
unchanged. The primitive screening of ERD bounds the $1/r_{12}$ integrals,
which are larger than the erf ones, so it stays safe.

`xerdints` now links `intgrt` and `sbtools` for `RSOMEGA`.

#### Could the existing ERD have done this?

Not as it was, and with a small change yes.

- As it was: `aces2/erd` and `aces2/erdints` compute $1/r_{12}$ only. The
  ACES III `erd_F12` routines (F12, K12, Yukawa operators) have no erf
  operator and are not used by ACES II. OED computes one-electron
  integrals and is not involved.
- With a change: Rys quadrature makes the erf attenuation simple (three
  scalings, above), about as much code as `RSATTN` in VMOL. The work is
  the attenuation, the s/p routing, and a second pass in `xerdints`.

| | VMOL | ERD |
|---|---|---|
| Change to the integral code | `RSATTN`, two calls | `ERD__RS_ATTENUATE`, (ss\|ss), routing |
| When it is used | always (default) | `INTEGRALS=ERD` |
| Long-range derivative integrals (gradients) | none in VDINT; new work | the same three scalings in `erd__derv_pcgto_block.f` would give them |
| KS gradient program | `xvdint` | `xerdgrad` is HF-only; KS (grid) terms would have to be added |

For energies the two routes are equivalent. For future analytic
gradients of range-separated hybrids, ERD is the better starting point,
because its derivative integrals use the same Rys quadrature; the
`becke_sr` gradient terms and KS support in `xerdgrad` would still be
needed.

### VSCF: the Kohn–Sham Fock matrix

`mkrhff.F` and `mkuhff.F` already build $J + c_\mathrm{hf}K$, with
$c_\mathrm{hf}$ in `coef_pot_nonlocal`. For range-separated hybrids,
`initintgrt` also sets `coef_pot_lrx` $=\beta$ (`coef_pot_exch(fun_exch_hflr)`),
and the Fock builders call

```fortran
      CALL LRXFOCK(F,D,BUF,IBUF,NTOTAL,NBAST,IMAP,ILNBUF,LUINT,
     &             coef_pot_lrx*FOURTHM)
```

(for UHF, once per spin with `coef_pot_lrx*HALFM`). The factors follow the
conventions of the exchange terms in the same routines. `libr/lrxfock.F`
reads `IIIILR`, and `IJIJLR` when `NIRREP>1`, and adds only the exchange
contractions:

$$
F_{ik} \mathrel{+}= f\,D_{jl}(ij|kl)_\mathrm{lr},\quad
F_{jl} \mathrel{+}= f\,D_{ik}(ij|kl)_\mathrm{lr},\quad
F_{il} \mathrel{+}= f\,D_{jk}(ij|kl)_\mathrm{lr},\quad
F_{jk} \mathrel{+}= f\,D_{il}(ij|kl)_\mathrm{lr} .
$$

The diagonal factors of ½ (for $i=j$, $k=l$, $ij=kl$) are the ones used for
`IIII`. `IIJJLR` and `IJKLLR` have no exchange contribution within one
symmetry-blocked density and are not read.

The long-range term is added only in KS iterations: not on the last
iteration (`scfkslastiter`) and not with `kshf`, just like the scaling by
`coef_pot_nonlocal`. Both builders stop with an error if long-range
exchange is requested together with SEWARD integrals (`petite_list`) or
direct SCF (`dirscf`), which do not produce the `*LR` files (VMOL and
`xerdints` do).

The exchange–correlation potential on the grid is assembled in
`integ_eff1.F`. A hybrid other than B3LYP sets `int_kspot=fun_special`, so
the potential is summed term by term from `coef_pot_exch` and
`coef_pot_corr`; the new terms `pot_exch_becke_sr` and `pot_corr_vwnrpa`
are called when their coefficients are nonzero. `hf_lr` has no grid
potential: its slot is zeroed.

### INTGRT: the final energy

`anlytc.F` sets `ehflx=0`. For the SCF density (`idns=0`) and
$\omega\neq0$ it calls `LREXCH`, which reads `IIIILR` and, with symmetry,
`IJIJLR`, and contracts them exactly as `exchange.F` does for the full
exact exchange. It returns $E_x^\mathrm{lr}$ in `ehflx`, which is stored in
`/exch/` (`energ.com`).

`numintint.F` and `numinteff.F` loop over the grid exchange functionals up
to `fun_num_xgrid` (lda, becke, pbe, pw91 and becke_sr; `hf_lr` and `hf`
are not grid functionals), then set

```fortran
      tot_exch(fun_exch_hf)=ehfx
      tot_exch(fun_exch_hflr)=ehflx
```

The hybrid energies are now a generic loop over `hybparm`:

$$
E_\mathrm{hyb} = \sum_i c_{x,i}\,E_{x,i} + \sum_j c_{c,j}\,E_{c,j} .
$$

The long-range integrals exist only for the current `rs_omega`, so a
range-separated hybrid whose $\omega$ differs from `int_rs_omega` cannot be
evaluated. It is skipped, and it is not printed in the table of hybrid
energies.

## Code changes

### New routines

| File | Purpose |
|---|---|
| `vmol/rsattn.F` | Attenuates $\rho$ and the prefactor for a batch of quartets |
| `vmol/rsfiles.F` | Switches units 2/81/82/83 to the `*LR` files |
| `intgrt/rsomega.F` | Determines $\omega$ from `*INTGRT` (VMOL and `xerdints`) |
| `erd/erd__rs_attenuate.F` | Attenuates $T$, the prefactor and the Rys roots |
| `erd/erd__rsomega.inc` | ERD common block `/ERDRSO/` with $\omega$ |
| `libr/lrxfock.F` | Long-range exchange in the Fock matrix |
| `intgrt/hybparm.F` | The table of named hybrids |
| `intgrt/lrexch.F` | Long-range exact-exchange energy |
| `intgrt/ityh_attn.F` | $F(a)$ and $dF/da$ |
| `intgrt/becke_sr_spin.F` | Short-range Becke energy density of one spin and its derivatives |
| `intgrt/func_exch_becke_sr.F`, `pot_exch_becke_sr.F` | Energy and potential of `becke_sr` (α and β spins; RHF copies α) |
| `intgrt/func_corr_vwnrpa.F`, `pot_corr_vwnrpa.F` | Energy and potential of `vwn_rpa` |

### Changed routines and include files

| File | Change |
|---|---|
| `vmol/conlor.F`, `vmol/conlox.F` | `COMMON /RSINT/ RSOMEG`; call `RSATTN` when `RSOMEG≠0` |
| `vmol/molecu.F` | Second `TWOEL` pass for the long-range files |
| `libr/mkrhff.F`, `libr/mkuhff.F` | Call `LRXFOCK`; error with SEWARD or direct SCF |
| `intgrt/include/dftfunc.com` | New indices (see below) |
| `intgrt/include/intgrtflags.com`, `include/intgrtflags.com` | `/intgrtflagsd/` gains `coef_pot_lrx` and `int_rs_omega` |
| `intgrt/include/energ.com` | `ehflx` at the end of `/exch/` |
| `intgrt/initintgrt.F` | Names of the new terms; hybrids read from `hybparm`; `rs_omega`; checks |
| `intgrt/integ_eff1.F`, `integint.F`, `integ.F` | Potential and energy of `becke_sr` and `vwn_rpa`; zero all slots |
| `intgrt/numinteff.F`, `numintint.F` | `ehflx` in `tot_exch`; generic hybrid loop |
| `intgrt/anlytc.F` | Call `LREXCH` |
| `vksdint/integAG.F`, `vksdint/numintAG.F` | Zero all slots; error for range-separated gradients |
| `erd/erd__e0f0_pcgto_block.f` | Call `ERD__RS_ATTENUATE` before and after the roots |
| `erd/erd__ssss_pcgto_block.f` | (ss\|ss): $\sqrt f\,F_0(fT)$ |
| `erd/erd__gener_eri_batch.f`, `erd__memory_eri_batch.f` | With $\omega\neq0$, s/p quartets other than (ss\|ss) go to `ERD__CSGTO` |
| `erdints/erd_driver.F` | Second `Erd_so2e` pass for the long-range files |
| `Makefiles/GNUmakefile` | `xvmol` links `intgrt`; `xerdints` links `intgrt` and `sbtools` |

The indices in `dftfunc.com` are now:

| Kind | Indices |
|---|---|
| exchange (`fun_num_exch=7`) | 1 lda, 2 becke, 3 pbe, 4 pw91, 5 becke_sr, 6 hf_lr, 7 hf |
| grid exchange (`fun_num_xgrid=5`) | 1–5 |
| correlation (`fun_num_corr=7`) | the existing six, 7 vwn_rpa |
| hybrids (`fun_num_hyb=7`) | 1 b3lyp, 2 cam-b3lyp, 3 cam-qtp-00, 4 cam-qtp-01, 5 cam-qtp-02, 6 lc-qtp, 7 qtp17 |

`hf` stays the last exchange index and the grid functionals come first, so
loops "over the grid functionals" run to `fun_num_xgrid`.

The build does not track include dependencies. After changing
`dftfunc.com`, `intgrtflags.com` or `energ.com`, run `ppclean` in every
library that includes them (`intgrt`, `vmol`, `oeptools`, `intpack`,
`vksdint`, `vscf`), rebuild, and relink the programs.

## Input

In `*INTGRT`:

| Keyword | Meaning |
|---|---|
| `func=name` | Functional(s) for the final energy |
| `kspot=name` | Functional for the KS potential |
| `rs_omega=x` | $\omega$ in bohr⁻¹; default: the $\omega$ of the hybrid named in `func`, else in `kspot`; otherwise 0 |

A named hybrid is normally given in both, for example

```text
*INTGRT
 kspot=cam-qtp-01
 func=cam-qtp-01
```

Other range-separated functionals are written term by term with
`hf_lr` and `becke_sr` and an explicit `rs_omega`. Test `lcblyp_h2o_rks`
builds LC-BLYP ($\omega=0.33$) this way.

Errors stop the run in these cases:

| Message | Cause |
|---|---|
| `rs_omega must not be negative` | $\omega<0$ (INTGRT, VMOL, `xerdints`) |
| `becke_sr and hf_lr need rs_omega` | `hf_lr` or `becke_sr` with $\omega=0$ |
| `func and kspot name hybrids with different rs_omega` | e.g. `func=cam-b3lyp`, `kspot=lc-qtp` |
| `range-separated hybrids need the long-range integral files of VMOL or xerdints` | SEWARD or direct SCF |
| `no analytic gradients for range-separated hybrids` | gradients or geometry optimization (`numinteff`, `numintAG`) |

## Limits

- **Energies only.** VDINT has no derivative integrals of
  $\mathrm{erf}(\omega r_{12})/r_{12}$, and `becke_sr` has no gradient
  terms, so gradients stop with an error.
- **VMOL or ERD integrals only.** SEWARD and direct SCF do not write the
  `*LR` files.
- **One $\omega$ per run.** The hybrid table printed by INTGRT omits
  range-separated hybrids with another $\omega$.
- **Symmetry is supported.** The tests run in C2v, and an H₂O run without
  symmetry gives the same energy. Point groups larger than C2v have not
  been tested separately.

## Verification

- **Long-range integrals.** An HF calculation with every integral file
  replaced by its long-range version matches a PySCF HF calculation with
  `omega=0.33` exactly, both without symmetry and for C2v/cc-pVTZ
  (−113.3001788926).
- **becke_sr derivatives.** `becke_sr_spin` was checked against finite
  differences, and against a 40-digit derivative for the
  $|\nabla\rho|$ term, where finite differences are too noisy.
- **ITYH series.** At $a=2$ the eight-term series matches a 40-digit
  evaluation of the closed form of $F(a)$ to a relative $6\times10^{-17}$.
- **VWN_RPA.** It agrees with libxc `LDA_C_VWN_RPA` to $10^{-15}$.
- **Total energies.** cc-pVDZ, the test-case grids, compared with PySCF;
  differences in hartree, consistent with the different grids:

| Functional | H₂O RKS | CH₂ triplet UKS |
|---|---|---|
| CAM-B3LYP | 1.5e-7 | 2.3e-7 |
| LC-BLYP | 6.8e-8 | – |
| CAM-QTP-00 | 2.7e-8 | 1.8e-7 |
| CAM-QTP-01 | 3.0e-8 | 2.2e-7 |
| CAM-QTP-02 | 6.6e-8 | 2.1e-7 |
| LC-QTP | 9.8e-8 | 1.5e-7 |
| QTP17 | 9.0e-8 | 1.9e-7 |

- **ERD against VMOL.** The same inputs with `INTEGRALS=ERD` give the
  VMOL energies to $5\times10^{-11}$ hartree or better, so the PySCF
  comparison above carries over:

| Functional | VMOL | ERD − VMOL |
|---|---|---|
| CAM-B3LYP, H₂O RKS | −76.391854413937 | 3.7e-11 |
| CAM-B3LYP, CH₂ UKS | −39.121407481873 | −2.7e-11 |
| LC-BLYP, H₂O RKS | −76.265013867952 | −4.8e-11 |
| CAM-QTP-00, H₂O RKS | −76.396362527926 | 2.6e-11 |
| CAM-QTP-01, H₂O RKS | −76.352853274522 | 2.2e-11 |
| CAM-QTP-02, H₂O RKS | −76.292658514275 | −2.5e-11 |
| LC-QTP, H₂O RKS | −76.266733480239 | 3.9e-11 |
| QTP17, H₂O RKS | −76.151936309928 | 2.8e-11 |
| QTP17, CH₂ UKS | −38.968733793801 | 1.0e-12 |
| CAM-QTP-00, H₂O C1, cc-pVTZ | −76.4312085945 | < 1e-10 |

  The cc-pVTZ run (f functions, no symmetry) checks the general Rys path
  at higher angular momentum and the single-file (C1) layout.
- **Regression tests.** The new `aces2/tests/fast/cases` entries are
  `b3lyp_h2o_rks`, `b3lyp_ch2_uks`, `camb3lyp_h2o_rks`, `camb3lyp_ch2_uks`,
  `lcblyp_h2o_rks`, `camqtp00_h2o_rks`, `camqtp01_h2o_rks`,
  `camqtp02_h2o_rks`, `lcqtp_h2o_rks`, `qtp17_h2o_rks`, `qtp17_ch2_uks`,
  and, with `INTEGRALS=ERD`, `camb3lyp_h2o_rks_erd` and
  `camb3lyp_ch2_uks_erd`. The aces2 suite passes 171/171, pyaces 60/60,
  and py-api 8/8.
