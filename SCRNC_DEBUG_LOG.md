# scrnc / GW manuscript debug log

Plain-text resume document, kept in the repo so it survives terminal
freezes / HPG outages and gets backed up on every `git push`. Read the
**RESUME POINT** section first on any restart.

---

## RESUME POINT (as of 2026-08-14, latest)

**Immediate next step**: dig into `aces2/scrnc/scrnc_jj_rcc_lrspn.F`'s
`Iside=2` (L-side) branch specifically. This is now the conclusively
isolated location of the JJ-matrix asymmetry bug (see "Where the bug is
isolated to" below). No decision needed to start this — just go read/trace
`Scrnc_form_mubar_s(Iside=2)` and the `Hbar^(-1)^T` transpose logic in that
routine.

**Working tree state**: uncommitted changes exist in (`git status` in this
repo confirms these, nothing lost):
- `aces2/ringcc/rcc_driver.F` (+ preprocessed `.f`)
- `aces2/scrnc/scrnc_driver.F` (+ `.f`)
- `aces2/scrnc/scrnc_form_rcc_hbar.F` (+ `.f`)
- `aces2/scrnc/scrnc_jj_cc_lrspn.F` (+ `.f`)
- `aces2/scrnc/scrnc_jj_rcc_lrspn.F` (+ `.f`)
- `aces2/scrnc/scrnc_lrspn_rcc_main.F` (+ `.f`)

Also 7 local commits ahead of `origin/master`, not yet pushed.

**Standing rule, repeated across this whole saga — ALWAYS do this before
testing any change**: `make` alone does NOT update `aces2/bin/x*`. You must
run `make install` in the changed module's directory AND in any module that
depends on it (e.g. ringcl depends on ringcc's library). Verify with
`md5sum aces2/lib/libX.a aces2/<module>/libX.a` if in doubt. Also: the
`aces2/bin/x*` executables in this tree are `ifx`-built and need
`module load intel/2025.1.0` (NOT the default `intel/2020.0.166`) loaded in
the SAME shell command as the run, or you get
`libimf.so: cannot open shared object file`.

---

## The actual goal

Ajith's own 1994-era idea (from NMR spin-spin coupling response theory,
generalized from one-particle/spin to two-particle/Coulomb quantities) for a
coupled-cluster-level GW / screened-Coulomb method: **W = v + v·χ·v**,
computed via CC linear-response machinery, to produce "a new set of
correlated two-electron integrals" (the user's own framing of the goal).
Module: `aces2/scrnc/` (executable `xscrnc`). Not yet a manuscript — physics
identified, code being debugged first.

**Framing correction the user gave (important, don't drift from this):**
- The **real production target is CCSD-level GW** (general coupled-cluster
  screened Coulomb, via the `xvcc`+`xlambda` path). This is **already
  confirmed exactly symmetric** — see "CCSD control case" below.
- The **ring-CC/DRCCD path inside scrnc is NOT a separate production
  method** — it's the RPA *limit* of the same machinery, kept specifically
  because RPA-GW is externally published, so it's the one place scrnc's
  output can be checked against a known answer.
- scrnc gets its RPA by routing through **ringcc's own** Hbar/T2 machinery
  — NOT by an externally re-derived textbook TDHF A/B matrix. Building an
  A/B matrix from scratch (tried once) is the WRONG approach — it bypasses
  the real architecture and gave a misleading "45% different" result.

Literature anchor: Bintrim & Berkelbach, JCTC,
https://pubs.acs.org/doi/full/10.1021/acs.jctc.6c00378 — establishes GW ≡
EOM-drCCD, i.e. scrnc's W=v+v·χ·v-via-drCCD-response is a real published
identity, not something being invented from scratch.

---

## Where the bug is isolated to (the key finding)

**CCSD (non-ring) path — control case, exactly symmetric:**
H2/STO-3G/UHF/CCSD/DENSITY=RESPONSE
(`/blue/ufhpc/perera/scrnc_test/h2_scrnc_ccsd_uhf_sto3g`, pipeline
xjoda...xvcc,xlambda,xscrnc). Needed `CC_EXTRAP=OFF` + `sym=off` to route
around H2/STO-3G's 1-amplitude DIIS degeneracy (unrelated to scrnc, happens
inside xvcc/xlambda). Final 4x4 JJ matrix is **exactly symmetric to 11
digits**. This path reuses the already-solved R-side amplitude and applies a
bra-side transform for the L-side — no second matrix build/inversion
anywhere.

**DRCCD/ring-CC path — the buggy one:**
H2/STO-3G/UHF/DRCCD/DENSITY=RESPONSE/sym=off
(`/blue/ufhpc/perera/scrnc_test/h2_scrnc_drccd_uhf_sto3g`). The final JJ
matrix is asymmetric with a very clean signature: the diagonal diff
(JJ(1,1) vs JJ(4,4)) and the off-diagonal diff (JJ(1,4) vs JJ(4,1)) are
**exactly equal** — looks like a single additive offset tied to atom-2's
AO-pair index, not generic numerical noise.

Two bugs were found and fixed along the way (both real, both kept,
neither fully explains the asymmetry alone):
1. **Hbar double-build bug (FIXED, uncommitted)**: `scrnc_form_rcc_hbar.F`
   destructively overwrites the same list it reads raw integrals from
   (`PUTALL(...,LISTDUMP)` where `LISTDUMP=LISTAA`). Since R-side
   (`scrnc_lrspn_rcc_main.F`) and L-side (`scrnc_jj_rcc_lrspn.F`) each
   independently call `Scrnc_rcc_hbar_inv`, the second call read back the
   first call's *finished* Hbar instead of pristine raw integrals, double
   counting one full orbital-energy gap. Fixed by refactoring so Hbar⁻¹ is
   built ONCE in `scrnc_driver.F` and shared via new `Hbarinv`/`Iloc`
   arguments into both routines. Confirmed correct via instrumentation
   (build now fires once, not twice; diagonal comes out right).
   **Fixing this alone made the asymmetry WORSE, not better** — proof a
   second, independent bug was partially canceling it out.
2. **Decisive elimination test**: routed T2 through the *independent*
   native RPA module (`aces2/rpa`, executable `xrpa`) instead of ringcc's
   own iterative DRCCD solve — full new pipeline
   `xjoda,xvmol,xvmol2ja,xvprops,xvscf,xvtran,xintprc,xrpa,xringcc,xringcl,
   xscrnc` (note: `xrpa` must be inserted manually between `xintprc` and
   `xringcc`; ACES2's automatic pipeline builder doesn't do this even
   though the machinery to consume it exists). Needed two more real fixes
   in `ringcc/rcc_driver.F`'s `RPA_T2VECS` bypass branch (both uncommitted,
   see git diff):
   - lists 61/62/63 (where ringcl reads its starting T2) were left zeroed
     by `rcc_setlst.F` because the bypass never wrote into them — fixed by
     copying 44→61, 45→62, 46→63 with correct per-list `IDSYMSZ` sizing.
   - `drcc_tener.F` needs raw `<AB||IJ>` integrals on lists 1/2, normally
     built by the (skipped) `RCC_DRIVE_D2T2` — fixed by inlining just the
     list-1/2 reorder snippet from `rcc_make_t2.F` (NOT the whole routine,
     which would clobber the real T2 with a throwaway MBPT2 guess).
   - Test case: `/blue/ufhpc/perera/scrnc_test/h2_scrnc_rpa_t2_uhf_sto3g/`.
   - **Result: asymmetry persisted, same signature, even with T2 sourced
     entirely from the independent RPA module.** This conclusively rules
     out ringcc/ringcl's own iterative solve as the cause.

**Conclusion**: the bug is isolated to **scrnc's own L-side JJ-matrix
construction**, `aces2/scrnc/scrnc_jj_rcc_lrspn.F`'s `Iside=2` branch
specifically. Its architecture has no counterpart in the CCSD path: it does
its own fresh `Scrnc_form_mubar_s(Iside=2)` computation (not reusing a
previously solved amplitude) and then separately applies the `Hbar^(-1)^T`
transpose (`Xgemm("T","N",...)`) to list `Mu_s=482`. This architectural
mismatch is the prime suspect. **Not yet traced line-by-line** — the actual
next step.

Also found and confirmed genuinely dead/orphaned (not related, don't touch
unless deliberately picking this up separately):
`ringcl/rcl_check_rcchbar.F` + `ringcl/rcl_drvrpa.F` — a fully-written but
never-called Hbar-diagonalization self-check (`#ifdef _DEBUG_LVL0`), zero
callers anywhere. Could be wired into `rcl_driver.F` cheaply later if
useful as an independent cross-check, but not required for the current
line of investigation.

The `xrpa`-side `-D_RCC` compile flag (whether `rpa/` builds in
ringcc-compatible or dringcc-compatible mode) was flagged as a possible
mystery earlier but is **confirmed intentional by the user** — a deliberate
compile-time switch, not a bug or accidental entanglement. Don't re-flag
this.

---

## Practical test-running notes

- `scrnc_set_vars.F` hard-requires `REF=UHF` (RHF explicitly rejected) and
  `DENSITY=RESPONSE` — both load-bearing ZMAT keywords for any scrnc test.
- `CC_EXTRAP=OFF` (real keyword, confirmed via `include/flags.h`'s
  `h_IOPPAR_cc_extrapol` alias table) routes around H2/STO-3G's inherent
  1-amplitude DIIS/MINV degeneracy in xvcc/xlambda/xringcc — unrelated to
  scrnc itself, just a minimal-system numerics issue.
- Test artifacts kept on disk (not in git) at:
  - `/blue/ufhpc/perera/scrnc_test/h2_scrnc_ccsd_uhf_sto3g/` (symmetric
    control case)
  - `/blue/ufhpc/perera/scrnc_test/h2_scrnc_drccd_uhf_sto3g/` (asymmetric
    reproduction case)
  - `/blue/ufhpc/perera/scrnc_test/h2_scrnc_rpa_t2_uhf_sto3g/` (RPA-sourced
    T2, decisive elimination test)
  - `/blue/ufhpc/perera/scrnc_test/h2_rpa_verify_uhf_sto3g/` (built for an
    earlier planned Hbar-eigenvalue cross-check, not yet run through)
  - 11 further scrnc_test dirs from earlier sessions, see memory for detail
    if needed.

---

## Separately open, NOT part of this thread (don't conflate)

- `ringcc/ringcl` bug-hunt saga (List-56, NSZIE typo, dodiis0 ifx-O2 bug,
  RHF=UHF DRCCD self-consistency) — **fully CLOSED and committed**
  (`fc32434a`, `f9ce968b`, `6361f7f6`). Not related to the current
  scrnc-specific bug hunt above.
- JCP 2020 EOM-RCCD paper recovery — closed, committed `660f0e73`, separate
  ~0.3 eV Fock-dressing discrepancy noted but not investigated.
- Native RPA module (`aces2/rpa`) test suite — separate task, H2O/Ne crash
  bugs deferred there, unrelated to scrnc.

---

*This file is meant to be committed/pushed as a durable backup whenever
updated — ask before committing if unsure, but the file itself can always
be freely edited/appended.*
