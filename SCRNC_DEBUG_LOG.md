# scrnc / GW manuscript debug log

Plain-text resume document, kept in the repo so it survives terminal
freezes / HPG outages and gets backed up on every `git push`. Read the
**RESUME POINT** section first on any restart.

---

## SESSION 2026-08-16 (evening, continued) -- IRREPZR gate bug found+fixed (RESOLVED for its own scope); SEPARATE, bigger, non-symmetry-specific bug found (OPEN)

**Fix found and verified**: `scrnc_form_mubar_d_uhf.F`, the
`IRREPX.EQ.1 .AND. ISIDE.EQ.2 .AND. IEOMPROP.EQ.0` correction inside
the `DO 300 IRREPZR=1,NIRREP` loop (adds `HALF*MBAR00*T1` to the
doubles-like Z-block being built for `IRREPZR`) was missing the
`IRREPZR.EQ.IRREPX` gate that the structurally identical correction
just above it (lines 398-405) already has. Confirmed via matched
R-side/L-side checksums first (`Scrnc_form_mubar_00`'s
`Mexpect`/`Mbar00` and `Scrnc_form_mubar_d_uhf`'s own Iside=1 checksum
both agree exactly between the R-side driver's own call and the
L-side's internal redundant Iside=1 recompute, for every one of the 49
perturbations -- ruled out the "computed twice" double-build pattern
that caused earlier bugs in this project) before reading code further
and finding this specific missing gate by direct comparison against
the adjacent, correctly-gated block. Fixed by adding
`IRREPZR.EQ.IRREPX .AND.` to the condition.

**Regression-tested and confirmed NOT the cause of a separate finding
below**: git-stashed the fix (reverting to committed HEAD) and reran
the sym=off H2O CCSD control case fresh -- the asymmetry found below
was IDENTICAL with and without this fix, proving the two issues are
independent. Restored the fix afterward (`git stash pop`).

**Committed**: `scrnc_form_mubar_d_uhf.F`.

---

## SESSION 2026-08-16 (evening) -- SEPARATE bug: H2O/STO-3G/CCSD asymmetric even with sym=off, NOT a convergence artifact (OPEN, bigger than the symmetry bug above)

While regression-testing the IRREPZR fix above, parsed the
`h2o_scrnc_ccsd_uhf_sto3g` control case's own "Response-only JJ matrix"
(the SAME system the 2026-08-15 evening session's RESUME POINT claimed
was "exactly symmetric, max diff = 0.0") with the same rigorous
programmatic parser used for the sym=on case. **It is not symmetric**:
600 of 1176 nonzero off-diagonal pairs differ, up to 5.17e-4
(JJ(9,41)=0.36434921829 vs JJ(41,9)=0.36486641546, both independently
verified by hand against the raw printed matrix, not a parser
artifact). The worst-offending pairs disproportionately involve AO-pair
index 9 (=natural AO pair (2,2), O's 2s orbital self-pair) and 41
(=natural AO pair (6,6), one H's 1s self-pair) -- both DIAGONAL
self-pairs. This is almost certainly why every prior "exactly
symmetric" verification in this whole project's history never caught
it: EVERY SINGLE ONE used H2/STO-3G (1 occupied orbital per spin --
the OO block is a trivial 1x1 scalar, so an OO-block off-diagonal bug
is structurally invisible there). H2O has 5 occupied orbitals per spin
-- a genuinely non-trivial 5x5 OO block -- and is apparently the first
time this project has checked JJ symmetry for a system where that
block isn't degenerate.

**Ruled out, decisively**: numerical convergence noise. Reran with
`cc_conv=12,scf_conv=12` (vs the control's `cc_conv=9,scf_conv=10`) --
the asymmetry is bit-identical to the last printed digit
(5.172e-04 both times). A real convergence-noise effect would shrink
when tightening convergence by 3+ orders of magnitude; this didn't
move at all, meaning it's a deterministic, reproducible bug, not
numerical residue.

**Ruled out**: the same R-side/L-side "computed twice" double-build
check used to isolate the IRREPZR bug above (`Scrnc_form_mubar_00`
Mexpect/Mbar00, `Scrnc_form_mubar_d_uhf` Iside=1 checksums) -- both
agree exactly between R-side and L-side's internal recompute for this
sym=off H2O case too. Whatever is wrong is NOT in the "same quantity
computed twice inconsistently" pattern that explains most of this
project's prior bugs -- it's something that's WRONG THE SAME WAY both
times it's computed (invisible to an R==L cross-check), which then
shows up downstream when L-side and R-side get combined asymmetrically
in the final JJ Ddot.

**Live lead, not yet confirmed**: `scrnc_form_mubar_s.F`'s own header
comment on the Iside=2 OO-block construction says "L(ji) which is what
actually being evaluated is stored in L(ij) ... Moreover, L(ij)
[MIJBAR, Iside=2] has no use in the rest of the program." This appears
to be STALE/WRONG: `scrnc_form_mubar_d_uhf.F`'s own DO-300 loop DOES
consume `MIJBAR` for ISIDE=2 (`CALL XGEMM('N','T',...,MIJBAR(IOFFM),
NCOL,...)`, transposed). If MIJBAR(Iside=2) is stored pre-transposed
(per the comment, holding "L(ji)" labeled as "L(ij)") AND consumed with
an explicit 'T' flag downstream, that could be a correct
double-negative (deliberately stored transposed, deliberately
un-transposed on read) OR a genuine transpose-convention bug -- this
needs the SAME careful adjoint-style re-derivation used to find the
`Scrnc_bcktran_left.F` bugs earlier today, not a guess. Whether this is
really the OO-block issue, or something else entirely (VV-block MABBAR
has an analogous "no use" comment worth checking too, and the doubles
G-intermediates from lists 191/192 haven't been examined at all yet)
is still open.

**Why this matters more than the symmetry bug above**: this affects the
CCSD path's correctness for ANY system with a non-trivial OO or VV
block, independent of point-group symmetry -- i.e. it's not really a
"symmetry" bug at all, despite being found while chasing one. Turning
on Abelian symmetry (this session's actual goal) can't be fully
validated as "done" until this is fixed too, since even the sym=off
H2O baseline is now known-broken.

**Not yet done**: isolate which specific term (OO, VV, or the doubles
G-intermediate contributions) is responsible, via the same
backward-trace/instrumented-checksum methodology as always. Once fixed,
recheck the sym=on C2v case (its remaining asymmetry, if any, after
combining with the already-fixed IRREPZR bug) to see whether it's then
clean, or whether a third, genuinely symmetry-specific issue remains.

---

## SESSION 2026-08-16 (evening) -- real Abelian symmetry turned on for the first time ever, NEW bug found (OPEN)

Ajith's direction: "let's get abelian sym as rest of ACES II" -- every
single scrnc test case in this entire project (all sessions above) has
used `sym=off`. Built the first-ever real-symmetry scrnc test:
H2O/STO-3G/UHF/CCSD/DENSITY=RESPONSE, same ZMAT as the trusted
`h2o_scrnc_ccsd_uhf_sto3g` control case but with `sym=off` simply
removed (`/blue/ufhpc/perera/scrnc_test/h2o_scrnc_ccsd_uhf_sto3g_symon`).
ACES2 auto-detects C2v correctly (confirmed in xjoda output: "The
computational point group is C2v", `Iaopop(1..4)=(4,1,2,0)`,
`Irpdpdao(1..4)=(21,8,16,4)`, summing to Nao=7/Nao_pairs=49 as expected).
Pipeline runs to completion, no crash.

**Real bug found, immediately, as expected**: parsed the pre-`Symmet2`
"Response-only JJ matrix" block with a small Python script (`file`
reports this pipeline.log as "data" not "ASCII text" -- some stray
non-UTF8 byte trips plain `grep`'s binary-file heuristic silently;
`grep -a` or Python's `encoding='latin-1'` are needed, NOT a race
condition -- this explains several earlier "grep returned nothing" false
alarms this whole session, worth remembering). Found **203 of 210
nonzero off-diagonal pairs asymmetric**, magnitudes ~1e-4 to ~4e-4 (not
noise) -- e.g. JJ(6,21)=-0.08588314 vs JJ(21,6)=-0.08631766. Every
single asymmetric pair has BOTH indices in AO-pair range 1-21 (the
totally-symmetric A1 irrep block, per `Irpdpdao(1)=21`) -- this is not a
cross-irrep leak, it's an internal ordering/indexing inconsistency
*within* the symmetric block, which is exactly why it was invisible on
every `sym=off` (Nirrep=1) test run so far: there was only ever one such
block, trivially self-consistent with itself.

**One candidate ruled out, decisively, empirically** (not by
inspection): isolated the BARE (uncorrelated) AO Coulomb integral part
by subtracting the "Response-only JJ matrix" from the "CI-like JJ
matrix" (=bare+response, both pre-`Symmet2`) -- the bare part alone is
symmetric to ~1e-11 (parsing-precision noise from 11-digit text
formatting, not a real discrepancy), while 100% of the 4.3e-4-scale
asymmetry lives in the response part alone (identical whether measured
in the response-only or CI-like block). This rules out
`scrnc_prep_2ints.F` / `Scrnc_reord_cols.F` / `Scrnc_sympack_aoints.F`
(the AO-pair irrep-grouping/symmetry-packing machinery that builds
List_xxxx AND the `Irrep_ao_pairs` array both `scrnc_lrspn_cc_main.F`
and `scrnc_jj_cc_lrspn.F` rely on for their own `Npert_count`/`Irrepr`
bookkeeping) -- that machinery is confirmed CORRECT even under real
symmetry, not the bug source, despite an initial hand-derivation
(below, kept for the record since the reasoning pattern may still be
useful) that looked very plausible before this empirical check.

**Investigated and ruled out on paper** (`Scrnc_reord_cols.F`,
2026-08-16): worth recording briefly since the reasoning was concrete,
even though the empirical bare-integral check above shows it isn't the
active bug. `Scrnc_reord_cols` reorders G2's COLUMNS into irrep-grouped
order via a permutation array (`Iordr_ao`) built by scanning raw
(k,l) AO-pair positions and collecting matches per irrep -- then
SEPARATELY bubble-sorts the `Irrep_ao_pairs` VALUE array in place (the
data-side column swap matching that same bubble sort is guarded by
`#ifdef _NOSKIP`, apparently never defined -- dead code, but harmless
since the array values still get the correct final sorted values either
way). `Scrnc_sympack_aoints.F` then independently re-derives the ROW-side
packing per column using a completely different traversal (Nu/nu-outer,
Mu/mu-inner, via `Iaopop`/`IoffAO`) -- a real, structurally distinct
construction from the column-side `Iordr_ao` scan (k-outer, l-inner).
Flagged as a plausible row/column within-block-ordering mismatch --
turned out NOT to be the active bug (bare integrals are provably
correct), so this specific theory is dead, but the general lesson
("verify empirically before fixing a plausible-looking theory") is worth
keeping.

**Next, not yet done**: the bug is isolated to the CCSD response
construction itself (`scrnc_lrspn_cc_main.F`'s R-side amplitude solve
via `Lineqy`/`Plineqy`, and/or `scrnc_jj_cc_lrspn.F`'s L-side Mubar
construction/bra-transform) -- the same general territory as every
previous bug in this whole project, but this is a genuinely new
instance, first exposed by symmetry. `Parteom`/`Nodavid`'s effect on
whether the Davidson-iterative vs partitioned solve path is taken for
this ZMAT (no `PARTEOM` keyword set) hasn't been checked yet -- worth
confirming which path actually ran before assuming either is the
culprit. **On resume: do NOT re-derive a "read the code" plan from
scratch -- follow the same backward-trace-from-JJ-assembly methodology
that cracked every prior bug in this log** (matched checksums at each
stage, L-side vs R-side, per [[feedback_stay_on_agreed_debug_path]]).

Not yet checked: whether the ring-CC/DRCCD path shows the same or a
different signature under symmetry -- worth doing once the CCSD bug is
understood, since it may or may not be the same root cause.

---

## SESSION 2026-08-16 (later) -- Scrnc_bcktran_left MO->AO transform kernel fixed+verified

Ajith: "I intended to do a correct bcktrn... let's do it correctly." Scope
note: this session deliberately did NOT re-decide what scrnc_bcktran's
real caller/consumer will be (the AO-vs-MO self-energy direction question
from earlier this session is still open) -- it treats "finish the
transform kernel itself, correctly, verified" as a self-contained,
valuable unit of work regardless of how it eventually gets called.

**Reference established first**: read aces2/bcktrn/ (the general ACES II
MO->AO two-particle-density back-transform used for CCSD/MBPT gradients,
J.F. Stanton & J. Gauss 1991) as the trusted precedent. Confirmed via
sytraa.f ("Q(XX) = C(X,A) Z(A,I) C^T(X,I)") and a `grep -l AOOVRLAP` over
the whole module (zero hits) that the correct, established convention for
this class of transform is Gamma_AO = C * Gamma_MO * C^T with NO AO
overlap-matrix (S) correction -- this is NOT the algebraic inverse of the
AO2MO integral transform used elsewhere in scrnc (scrnc_i1tran.F's
ICASE=3, "M_MO = C^T M_AO C", which WOULD need an S-metric sandwich to
invert exactly in a non-orthogonal AO basis), it's the separate, standard
convention for re-expressing a density/response-like tensor in a
different one-particle basis. Settled a real worry (a naive AO-integral
round-trip test would NOT have been the right validation method, since
forward and backward here are genuinely different operations, not
inverses of each other).

**Two real bugs found in scrnc_bcktran_left.F**, by re-deriving it as the
explicit adjoint of scrnc_i1tran.F's ICASE=3 block (the trusted,
already-working local AO2MO reference, same AOSYM/SYM/SYMPOP offset
conventions):
1. The per-irrep offset into the input MO-pair array incremented by
   `Nsum` (=Pop(Irrep_j,Ispin) alone) instead of the full block size
   `Numi*Numj` -- would misread every irrep block after the first once
   more than one irrep contributes to a pair symmetry.
2. The column offsets used to slice the eigenvector matrix ("Ioff_j"/
   "Ioff_i") were raw MO population COUNTS (Pop(Irrep,Ispin)) instead of
   actual MO column OFFSETS (Ioffo(Irrep,Ispin) -- confirmed via
   scrnc_getaoinf.F's own construction, "ioffo(1,1)=1;
   ioffo(irrep+1,1)=ioffo(irrep,1)+pop(irrep,1)", a genuine 1-based
   cumulative offset, and scrnc_i1tran.F's own direct usage
   "EVEC1(IOFFL1,IOFFR1)" with no extra +1).
Both bugs are silently invisible on any single-irrep (sym=off) system --
Nirrep=1 makes every irrep-block offset trivially 0 regardless -- which
is why they'd escaped notice: every scrnc test case run so far this whole
project has used sym=off. They only bite once real point-group symmetry
is turned on.

**Also found and fixed**: the routine was hardcoded alpha-spin-only
(SCFEVECA, Pop(*,1) unconditionally) with no way to get beta -- wrong for
the UHF-only systems scrnc requires (Scrnc_set_vars.F rejects RHF
outright). Added explicit Ispin/Iuhf arguments, threaded through to the
SCFEVECA/SCFEVECB record selection (same LABELX(ISPIN) pattern
scrnc_i1tran.F already uses) and the Pop(Irrep,Ispin)/Ioffo(Irrep,Ispin)
lookups.

**Driver-level bugs found in scrnc_bcktran.F** (the (IJ,KL)->(XX,XX)
caller that drives Scrnc_bcktran_left twice per irrep, "transform left
pair, transpose, transform right pair" -- same pattern
aces2/bcktrn/bxaaaa2.f uses):
1. `I040 = Iend` used Iend before it was ever assigned (undefined-value
   read) -- that whole I040 buffer was never actually referenced again
   in the routine, so removed outright rather than patched forward.
2. `JJ_respn_ijkl(1,Icol)` always read starting at GLOBAL row 1, local
   column Icol=1..Nij_cols, for EVERY Irrep_r -- correct only for the
   very first irrep block; any irrep past the first silently re-read the
   first irrep's rows/columns instead of its own. Fixed with explicit
   cumulative per-irrep offsets (Nij_cumoff(Irrep), built once before the
   main loop, indexed directly by irrep number so it stays correct
   regardless of loop order or Irrepx). NOTE: this assumes
   JJ_respn_ijkl's GLOBAL occ-pair index is laid out as consecutive
   per-irrep blocks in irrep order -- scrnc's own invented convention (no
   other scrnc routine stores a dense all-irrep occ-pair matrix to check
   this against), flagged as a documented assumption, not an established
   fact, until a real caller exists to confirm it.
Also added Ispin as an explicit driver-level argument (using
Irpdpd(Irrep,20+Ispin) instead of hardcoded column 21 -- confirmed
correct convention via scrnc_xxxx_pqxx.F's own identical
"NUMIJ=IRPDPD(IRREPIJ,20+ISPIN)" for IFLAG=3).

**Verification**: added a permanent, cheap, self-contained brute-force
cross-check directly inside Scrnc_bcktran_left.F (guarded by
_DEBUG_LVL0, kept as a standing diagnostic, not exploratory cruft) --
independently recomputes each per-block result via direct nested-do-loop
evaluation of Result(X1,X2)=Sum_I Sum_J C(X1,I)*Z(I,J)*C(X2,J), sharing
nothing with the Xgemm code path except the raw offsets, and reports
max|brute-xgemm|. Also added a small, explicitly-temporary standalone
harness (new file Scrnc_bcktran_test.F, called once from Scrnc_driver.F
right after Scrnc_getaoinf, guarded by _DEBUG_LVL0) that feeds a
synthetic deterministic occ-pair matrix through Scrnc_bcktran_left for
both spins, since nothing wires real data into it yet -- header comment
marks it "safe to delete once a real caller exists."

Ran on H2O/STO-3G/UHF/sym=off (fresh full pipeline, xjoda through
xscrnc) -- a genuinely non-trivial case (25x25 occ-occ pair block per
spin, 49 AO pairs, unlike H2's degenerate 1x1 case):
```
Ispin=1 (alpha): max|brute-xgemm| = 0.166533D-15
Ispin=2 (beta):  max|brute-xgemm| = 0.333067D-15
```
Both at machine precision -- confirms the Xgemm call convention
(argument order, leading dimensions, transposes) and the UHF spin
threading are correct. **Still NOT independently verified**: the
multi-irrep Ioffo/block-size fixes (Nirrep=1 for every scrnc test case
run so far, including this one, can't exercise them) and the
driver-level scrnc_bcktran.F itself (only the Scrnc_bcktran_left kernel
was exercised by this test; the (IJ,KL)->(XX,XX) driver's global-offset
bookkeeping and its assumed Npairs ordering convention remain unverified
end-to-end).

**Committed**: scrnc_bcktran.F, scrnc_bcktran_left.F, new
scrnc_bcktran_test.F, scrnc_driver.F (temporary harness call), this log.

**Not yet done**: this is still just the transform kernel -- the earlier
open question (does the self-energy pipeline need G0 backtransformed
MO->AO, or W_c transformed the other way AO->MO, or something else
entirely -- see the "bcktran target" discussion above) is UNRESOLVED and
should be settled before building the real caller. Once decided, the
driver's Npairs-ordering assumption should get a real end-to-end check
(not just the kernel's synthetic-data test), and a genuine multi-irrep
(symmetry-on) test case would close the remaining verification gap on
the Ioffo/block-size fixes.

---

## SESSION 2026-08-16 -- imaginary-frequency axis, both paths (WORKING, committed)

Picked up per Ajith's direction ("do 1" = build the GW self-energy,
"use your own judgement"). Both 2026-08-15 dynamic-W(omega) sessions
independently flagged the same open item: "imaginary-frequency axis not
yet tested." Investigated first rather than diving straight into
self-energy code, since any spectral/contour-based self-energy
construction needs W(i*omega) as its actual numerical input (standard
numerical-GW contour-deformation technique) -- so this was a genuine
prerequisite, not a detour.

**Real bug found on the CCSD path (not just "untested," actually
broken)**: `vcceh/secprp.F` -- the existing, working, presumably
decades-validated dynamic-polarizability code this whole scrnc dynamic-W
effort mirrors -- has its own explicit convention at line 124:
`IF(.NOT. REALFREQ) FREQ = FREQ*FREQ`, i.e. the imaginary-axis magnitude
must be SQUARED by the caller before being passed into `Lineqy`/
`Plineqy`. Confirmed why: `vcceh/lineqy.F` (lines 143,242,277) divides by
`(denom**2+FREQ)` on its `.NOT.REALFREQ` branches -- it expects `FREQ` to
already be `omega**2`, not raw `omega`. `scrnc_lrspn_cc_main.F` (added
2026-08-15) passes `Freq`/`Realfreq` straight through to `Lineqy`/
`Plineqy` unchanged but never did this squaring -- so any imaginary-axis
scrnc run was silently dividing by `(denom**2+omega)` instead of the
correct `(denom**2+omega**2)`. Confirmed numerically: baseline (buggy)
H2/STO-3G/CCSD run at `omega=0.3i` gave response(1,1)=0.04522885127;
after the fix (below), the same input gives 0.05449824469 -- a real,
non-trivial difference, not a rounding-level discrepancy.

**Fix (CCSD path, `scrnc_lrspn_cc_main.F`)**: added
`If (.NOT. Realfreq) Freq = Freq * Freq`, once, right after the
`Itop`/`Fact` determination and before the `Do Icount` loop -- same
placement `secprp.F` uses relative to its own `ITOP`/`FACT` logic.
Confirmed no other consumer in scrnc needs the un-squared magnitude
(`Freq` is used only for the `Lineqy`/`Plineqy` calls and the real-axis
`Freq = -Freq` negation; `Xfreq(:)` array, unaffected, is unused apart
from being populated at read time).

**Ring-CC path had no imaginary-frequency handling AT ALL** (not a typo
this time -- a missing branch): `scrnc_form_rcc_hbar_inv.F` unconditionally
did `Hbar_diag += Freq` on any nonzero `Freq`, regardless of `Realfreq`,
then a single `MINV`. That is the correct object for the REAL-axis
single-pass cases (static `Freq=0`, or one leg of the already-verified
+/-Freq two-pass sum) but is not the right matrix for the imaginary axis,
because ring-CC builds Hbar^(-1) as a genuine full (non-diagonal) matrix
inverse rather than `Lineqy`'s per-element denominator shortcut.

**Derivation (linear algebra identity, not an approximation)**: the
already-verified real-axis two-pass sum is
`(Hbar-Freq)^(-1) + (Hbar+Freq)^(-1) = 2*Hbar*(Hbar**2-Freq**2)^(-1)`.
Substituting `Freq -> i*sqrt(Freq)` (purely imaginary) and summing the
+i*omega/-i*omega poles the same way gives the real-valued
`2*Hbar*(Hbar**2+Freq)^(-1)` (`Freq` already squared by the caller, same
convention as the CCSD fix). This collapses to plain `Hbar^(-1)` at
`Freq=0`, exactly matching the existing static case
(`Hbar*(Hbar**2)^(-1) = Hbar^(-1)`). As with the real-axis branch, the
existing `Fact=2` in `Scrnc_jj_rcc_assemble.F` supplies the overall
factor of 2 -- only `Hbar*(Hbar**2+Freq)^(-1)` itself needed building
here.

**Fix (ring-CC path)**: `Scrnc_driver.F` squares `Freq` once (same
placement/convention as the CCSD fix) before the dynamic Icount loop.
`scrnc_form_rcc_hbar_inv.F` now branches three ways: `REALFREQ` nonzero
(unchanged, already-verified `Hbar+-Freq` single-shift+MINV), imaginary
nonzero (new: build `Hbar**2` via `Xgemm`, add `Freq` to its diagonal,
`MINV` that, then `Xgemm` the original (untouched) `Hbar` against the
inverse to get `Hbar*(Hbar**2+Freq)^(-1)`, `Dcopy` into the output slot),
and `Freq=0` (unchanged plain `Hbar^(-1)`). New scratch (`Hbar**2`
buffer, product buffer, MINV pivot array) allocated from whatever `IEND`
already is at that point in the routine -- correct under both `_DEBUG_LVL0`
on and off, since those conditionally-compiled diagnostic blocks already
advance `IEND` before this code runs.

**Verified** (H2/STO-3G/UHF, fresh full pipelines from ZMAT+GENBAS every
time, `/blue/ufhpc/perera/scrnc_test/h2_scrnc_{ccsd,drccd}_*imagfreq*`):
- CCSD path (response(1,1), before bare-v add): static
  0.05974824293; omega=0.001i: 0.05974817899 (~static, continuous);
  omega=0.3i: 0.05449824469; omega=1.0i: 0.02886296549.
- Ring-CC/DRCCD path: static 0.03008849944 (fresh regression rerun,
  bit-identical to the already-committed baseline -- confirms the
  three-way branch split didn't disturb the real/static path); omega=0.001i:
  0.03008847497; omega=0.3i: 0.02800114343; omega=1.0i: 0.01520543542.
- Both paths: smooth, monotonically DECREASING response magnitude as
  |omega| grows on the imaginary axis (opposite trend from the real
  axis's increasing-toward-resonance behavior, exactly as expected --
  the imaginary axis never crosses a pole) -- continuous through
  omega->0, no crashes, no NaNs, full W(omega) matrix ("LR JJ matrix" /
  "RPA JJ matrix") exactly symmetric at every frequency tested on both
  paths.
- This is a strong but NOT independently-cross-validated check (no
  external reference imaginary-axis polarizability was compared against
  -- unlike the real axis, which had the ~0.55 a.u. RPA excitation energy
  as an external landmark to approach). The qualitative signature
  (smoothness, monotonic decrease, correct static limit, symmetry) is
  exactly what a correct implementation must produce and what several
  classes of bugs (sign errors, missing terms, wrong shift target) would
  very likely have broken, but it is not a substitute for an independent
  numerical cross-check if one becomes convenient later (e.g. comparing
  against `vcceh/secprp.F`'s own imaginary-axis dynamic polarizability
  for some ordinary property, on a system where that path is already
  trusted).

**Committed**: `scrnc_driver.F`, `scrnc_form_rcc_hbar_inv.F`,
`scrnc_lrspn_cc_main.F`, plus this log.

**Not yet done / actual next step**: the imaginary-frequency machinery
validated here is the necessary numerical input for a contour-deformation
GW self-energy (`Sigma_c(p,q;i*omega)` built from `W_c(p,q;i*omega')`
integrated along the imaginary axis then analytically continued to the
real axis, e.g. via Pade) -- this session did NOT build any of that:
no AO-to-MO transform of the correlation part of `W(omega)` into the
particle-hole/quasiparticle index pattern a self-energy contraction
needs, no frequency quadrature, no analytic continuation, no quasiparticle
equation solve. `W(omega)` itself is still purely in the AO-pair basis
(`JJ_lrespn`, `Nao_pair x Nao_pair`) on both paths. Two genuinely
different architectures were considered for the self-energy step itself
and NOT decided between yet:
1. Numerical contour-deformation using `W(i*omega')` sampled along the
   imaginary axis (what this session's fix directly enables) -- avoids
   ever diagonalizing Hbar, stays closest to "W computed via CC linear-
   response machinery" as stated in the goal, but needs a frequency
   quadrature grid, analytic continuation (Pade or similar) from
   imaginary to real axis, and the AO-to-MO transform of `W_c` either
   way.
2. Sum-over-poles/spectral representation using ring-CC's own Hbar
   eigenvalues+eigenvectors as the poles of `W_c` (the standard
   molecular-GW residue formula, e.g. Bruneval/van Setten style,
   `Sigma_c^pq(omega) = sum_s sum_i (pi|s)(qi|s)/(omega-e_i+Omega_s-i*eta)
   + ...`) -- `ringcl/rcl_drvrpa.F` already computes Hbar eigenvalues AND
   left/right eigenvectors via `DGEEV` (confirmed 2026-08-16 while
   investigating this session's bug) but currently discards the
   eigenvectors immediately after computing them (only eigenvalues are
   ever printed) -- would need those captured/stored first.
Route 1 is closer to the architecture already built and to the stated
goal's own framing ("computed via CC linear-response machinery"); route 2
reuses the already-external-validated (EOM(Sf)-DRCCD vs ringcc's own
eigenvalues) Hbar diagonalization but requires wiring up eigenvector
storage that doesn't exist yet and moves further from a "pure CC
response" construction toward a conventional RPA-diagonalization one.
**Ask Ajith which he prefers before starting the self-energy contraction
itself** -- this is a real fork in the physics/architecture, not an
implementation detail, and picking wrong means redoing real work.

---

## SESSION 2026-08-15 -- backward logic trace on H2O/STO-3G (CLOSED -- see the two "dynamic W(omega)" entries below for what followed)

Per Ajith's direction: use H2O/STO-3G directly (already confirmed
asymmetric, no need to escalate basis), and do a full backward LOGIC
review (not just instrumentation) from JJ assembly back through the call
chain. Traced and read in full: `scrnc_jj_rcc_lrspn.F` (JJ assembly) ->
`scrnc_lrspn_rcc_main.F` (R-side, list 373) -> `scrnc_form_mubar_s.F` +
`scrnc_form_mubar_00.F` -> `scrnc_driver.F` (top-level, Hbar-once
build/Iloc) -> `scrnc_form_rcc_hbar.F` + `scrnc_form_rcc_hbar_inv.F` ->
`scrnc_prep_culomb_ints.F` -> `scrnc_form_lrspn_rcc_lists.F` ->
`scrnc_built_respn_den.F`. Not yet reached: `scrnc_dens_oo/vv/vo.F`,
`scrnc_xxxx_pqxx.F`, `scrnc_ao2mo.F`, `scrnc_getaoinf.F`,
`scrnc_make_lists.F`, `scrnc_prep_2ints.F`, `scrnc_form_rccd_g2.F`,
`scrnc_i1tran.F`, `scrnc_bcktran*.F` -- keep working backward from here on
resume, in that rough order (density-build routines next, then the AO2MO
transform chain).

**Two real (but probably non-causal) bugs found and NOT yet fixed** (do
not fix without deciding whether to bundle with the eventual root-cause
fix, or fix now separately -- ask Ajith which he prefers when resuming):
1. `scrnc_driver.F:144`, inside the `Ndim` accumulation loop (used only
   to size the Hbar^(-1) scratch-padding budget, `Ineed=Ndim*Ndim+Ndim`,
   `Ipad=10000+10*Ineed`): `Ndim = Ndim + Irpdpd(Irrepr-1,10)` should read
   `Irpdpd(Irrepr,10)` -- off-by-one, reads `Irpdpd(0,10)` (out of bounds)
   on the Irrepr=1 iteration. Protected by the 10x Ipad margin, so almost
   certainly harmless in practice, but worth fixing for correctness/
   robustness regardless. Does NOT affect `Iloc(:)`, which is built
   correctly in the loop just above it (verified: `Iloc(Irrepr) =
   Iloc(Irrepr-1) + Length` with `Length` from `Irpdpd(Irrepr-1,9)[+10]`
   squared -- that recurrence is right, each step adds the PREVIOUS
   irrep's own block size, which is the correct prefix-sum pattern).
2. `scrnc_form_lrspn_rcc_lists.F:61`, sizing lists 373/374/375/376:
   `Length(Irrep) = Length(Irrep) + Irpdpd(Irrep,9)` for the UHF case
   should add column `10` (beta PH-pair count), not `9` (alpha) a second
   time -- same copy-paste-off-by-one-column shape as the already-fixed
   `scrnc_form_mubar_s.F` I020 typo. Numerically inert for any genuinely
   closed-shell UHF reference (H2 and H2O both), since alpha/beta
   occ+virt counts match per irrep there, so `Irpdpd(Irrep,9) ==
   Irpdpd(Irrep,10)` always holds and the list ends up the right size
   anyway. Real risk only for open-shell UHF references, out of scope
   for this thread.

**Key structural finding, narrows the search a lot**: for RCCD/DRCCD
(ring-CC, no singles), `scrnc_rccd_dummy_inits.F` hard-zeros lists 90
(T1), 190 (L1), 191 (G_im), 192 (G_ea) every run. `scrnc_form_mubar_s.F`
reads exactly these lists to drive its OO-block, VV-block, and the
G-intermediate part of its VO-block (Iside=2 lines ~232-256) -- all zero,
so those terms contribute nothing for ring-CC regardless of Iside. This
means, for ring-CC specifically:
- `MIJBAR`/`MABBAR` are identically zero on BOTH sides (no singles at
  all) -- not a source of L/R asymmetry.
- The only nonzero content of `Scrnc_form_mubar_s`'s output, either side,
  is the T2-contraction term at lines ~273-318: `MAIBAR(ai) +=
  T2(ai,em)*MAI(em)` [+AA spin term], using list `38-ISPIN`/`33+ISPIN`
  (regular T2, right-hand amplitude, Iside=1) vs list
  `138-ISPIN`/`133+ISPIN` (offset +100, "L2"/Lambda-T2, left-hand
  amplitude, Iside=2).
- Also confirmed the `Mbar_00-M_expect` SAXPY term (line 321-324,
  `IF(ISIDE.EQ.2 .AND. IRREPX.EQ.1 .AND. IEOMPROP.EQ.0)`) is a genuine
  no-op for ring-CC too -- it scales list 190 (L1), which is zero. So a
  bug in `Scrnc_built_respn_den.F`'s density (which only feeds
  `M_expect`) CANNOT be the cause here; ruled out without needing to read
  `scrnc_dens_oo/vv/vo.F` in detail for this purpose (may still be worth
  a pass for its own sake, but not on this critical path).

**Where this points next**: the L/R asymmetry, if it's a code bug and not
a deeper physics/algebra gap, most likely lives in whether list
`133+ISPIN`/`138-ISPIN` ("L2", populated by `ringcl`'s own Lambda solve,
NOT by scrnc) uses the exact same DPD storage/ordering convention as list
`33+ISPIN`/`38-ISPIN` ("T2", populated by `ringcc`). This is a
cross-module convention-matching question (ringcc vs ringcl), same shape
as the original list-374 bug but one level removed -- NOT yet
investigated. Concretely: find where `ringcc` writes lists 33+ISPIN/
38-ISPIN (likely `ringcc/rcc_rng.F` or `ringcc/rcc_rng_4rhf_rccd.F`,
grep-confirmed both reference `LISTT=33+ISPIN`/`38-ISPIN`) and where
`ringcl` writes 133+ISPIN/138-ISPIN (`ringcl/rcl_l2rng.F`,
grep-confirmed `LISTT=133+ISPIN`/`138-ISPIN`), and diff their storage
conventions (row/col ordering, symmetry-packing, sign conventions)
line-by-line. If they match, the next fallback is the raw-value-dump plan
from yesterday (H2O AO-pairs 9/20/38/41/49) but now narrowed specifically
to the T2AB/L2AB values feeding the line-296 XGEMM in
`scrnc_form_mubar_s.F`, since everything else in that routine is now
confirmed zero for ring-CC.

**Raw-value-dump result (2026-08-15, DONE) -- no code bug found, this now
looks like a physics/design question.** Instrumented
`scrnc_form_mubar_s.F` (raw T2AB/T2AA/MAI/MAIBAR dumps, both sides) and
`scrnc_jj_rcc_lrspn.F` (raw dump of `Work(I050)` right after the Iside=1
call returns, before any Iside=2 setup runs). Rebuilt, reran fresh on
`/blue/ufhpc/perera/scrnc_test/h2o_scrnc_drccd_uhf_sto3g_rawtrace2/`.

**First pass at reading the output was WRONG** -- a self-inflicted parsing
bug, not a code bug: the `"L-routine Npert_count=N"` marker prints
*between* AO-pair N's Iside=1 trace and its Iside=2 trace (see code at
`scrnc_jj_rcc_lrspn.F:140-146`), so a naive "split the log by marker
line" parse pairs pair N's Iside=2 block with pair (N+1)'s Iside=1 block
by accident. That false pairing produced an apparent ~5-10x magnitude
mismatch between "Iside=1's MAIBAR output" and "Iside=2's MAI input" that
looked exactly like the original list-374 collision bug's signature --
worth remembering for next time: always sanity-check parsed line ranges
against `grep -n` + manual inspection before trusting a derived
"mismatch," especially when the file has interleaved multi-block prints
like this one.

**Corrected pairing (hand-verified against `xscrnc.out` line numbers,
AO-pair 9, lines 5915-5975)**: Iside=2's "MAI" input (`Work(I050)`,
read via `MAI(IOFFVO2)`) matches Iside=1's own `MAIBAR` output for the
SAME AO-pair EXACTLY, as the code intends (Iside=2 legitimately consumes
Iside=1's dressed Mbar as its own "bare" operator input -- confirmed
intentional, matches the doc-comment's overbarred `M(ei)` notation for
Iside=2 in the routine header). No aliasing/collision bug here.

The only real difference between the two sides is T2 (ringcc, list
`38-ISPIN`/`33+ISPIN`) vs L2 (ringcl, list `138-ISPIN`/`133+ISPIN`)
differing at the 6th-7th significant digit for the same AO pair (e.g.
AO-pair 9: T2AB(1)=-0.0007489825 vs L2AB(1)=-0.0007489964). This is
almost certainly REAL physics (T2 and L2 are genuinely different
amplitudes for a non-symmetric Hbar, not a bug), and it's plausible on
its face that this tiny per-component difference, propagated through the
`MAIBAR` contraction (~4e-5 level difference in the dressed value) and
then the `Hbar^-1`/`Hbar^-T`-multiply + 373-dot-374 sum over
~10-20 components, is fully sufficient to produce the observed ~1e-4
level JJ asymmetry -- i.e. this may not need a "bug" explanation at all.

Also confirmed (same trace, both sides, not new): for several AO pairs
the AA and AB xgemm contributions to `MAIBAR` exactly cancel (T2AA reads
bit-identical to T2AB at the `dis=1` slice), leaving `MAIBAR` equal to
the raw undressed `MAI`. Happens identically on R-side (trusted) too, so
NOT a source of the asymmetry -- likely a real closed-shell-UHF
degeneracy at this particular DPD slice, not a bug. Don't re-chase this.

**Bottom line after this full backward pass**: exhausted static logic
review of every routine in scrnc's own ring-CC call chain (see full list
above) plus a corrected, hand-verified raw-value trace at the one
remaining line with real physics content, and found NO further coding
bug. What's left is a genuine physics/design question, needs Ajith's
judgment, not more code reading:

- Is there supposed to be a specific mathematical relation between T2,
  L2, and Hbar's left/right eigenvectors that forces JJ to be EXACTLY
  symmetric (in which case something enforcing that relation is
  currently missing -- e.g. maybe ringcl's L2 solve needs an extra
  biorthogonality/normalization constraint against ringcc's T2 that isn't
  currently being applied before scrnc consumes both)?
- OR was H2's earlier "exact symmetry to 11 digits" itself a special
  case (H2/STO-3G's minimal, highly-degenerate amplitude structure)
  rather than a general guarantee, and a small residual asymmetry at the
  H2O scale (~1e-4, small relative to the ~0.5-1.3 magnitude of the JJ
  elements themselves) is actually the EXPECTED, correct behavior of
  this method?

Whichever answer, the next concrete step depends on it: if a missing
constraint, need to find where/how ringcl's L2 solve should be tied back
to T2 (biorthogonal Davidson normalization is the likely place to look,
`ringcl/rcl_*.F`, not yet investigated). If expected residual, this
thread is effectively DONE and the number itself (or its convergence
behavior vs cc_conv/basis size) becomes the thing to characterize instead
of chased as a bug.

**PHYSICS CALL, 2026-08-15 -- RESOLVED, no bug: this asymmetry is
expected, not a defect.** Ajith's ruling: like CC (and worse, since
ring-CC/RPA-type EOM has no explicit exp(T) similarity-transform form
constraining it), left and right EOM eigenvectors are genuinely different
even though they share the same excitation energy -- T2 (right/ringcc)
and L2 (left/ringcl) are a biorthogonal pair, NEVER identical in general.
H2's earlier "exact symmetry to 11 digits" was the special case (H2/
STO-3G's minimal, highly-degenerate amplitude structure let T2≈L2 there
by accident), not a general guarantee. The ~1e-4-level H2O JJ asymmetry
found this session is the CORRECT, EXPECTED behavior of the method, not
a code defect -- no further chasing needed on this specific thread.

**SUPERSEDED, 2026-08-15 -- there WAS still a real fix, just not the one
guessed above.** The "expected, no bug" verdict was based on a wrong
premise: I had assumed `ringcl`'s L2 was a state-specific EOM left
eigenvector (biorthogonal to T2, genuinely never equal, no exact-symmetry
guarantee). Ajith corrected this directly: `ringcl` is ring-CC's
GROUND-STATE Lambda solver (same role as `xlambda` for CCSD), not an
excited-state EOM solve -- confirmed by reading `ringcl/rcl_driver.F`
directly: no diagonalization anywhere in the real solve path (`RCL_FORM_
HBAR_AB/JI`, `RCL_DORNG`, `DRCL_DWMBEJ`, DIIS, `RCL_NEWT2`, convergence
test -- a plain iterative Lambda-amplitude solve, `'L'`-tagged). The only
`DGEEV` calls anywhere in `ringcl/` are in `rcl_drvrpa.F`/`rcl_heffxc.F`,
which are the orphaned diagnostic *we* wired in this project, not part of
production ringcl at all.

Given that, the CCSD-vs-ring-CC asymmetry difference needed a real
explanation, not a physics shrug. Found it by direct comparison with
`vcceh/secprp.F` (Stanton 1993 EOM-effective-Hamiltonian property code --
Ajith's own prior work in this exact area, predating/contemporaneous with
Koch & Jorgensen 1990). `secprp.F` builds its polarizability tensor
`XPROP` via essentially the same L/R-response-vector Ddot construction as
scrnc's JJ_lrespn (same `OPRBAR`/372+ISIDE list-pairing architecture) --
and then, at line 484, explicitly calls `Call Symmet2(XPROP,MAXPRT)`
before ever printing or using the result. `Symmet2` (`aces2/libr/
symmet2.f`, J. Gauss 1990) is a plain `A(pq)=0.5*(A(pq)+A(qp))`
force-symmetrization. This is the ACES2-established, decades-old standard
practice for exactly this class of L/R response construction: don't rely
on the raw construction to be exactly symmetric, symmetrize the finished
property matrix as the final step (same conceptual move as the complex-
conjugate symmetrization of a propagator, per Ajith's own framing).
scrnc's JJ routines were clearly modeled on `secprp.F`'s architecture but
never carried over this final step.

**FIX (2026-08-15, VERIFIED)**: added `Call Symmet2(JJ_lrespn,Nao_pair)`
to both `scrnc_jj_rcc_lrspn.F` (right before the final `_DEBUG_LVL0`
print, after the bare-integral add) and `scrnc_jj_cc_lrspn.F` (right
after the optional `Calcquad1` quadratic-contribution call), matching
`secprp.F`'s placement exactly. Rebuilt, reran fresh full pipelines for
all four available control/test cases:
- `h2o_scrnc_drccd_uhf_sto3g`: max asymmetry 0.0 (was 1.104E-04)
- `h2_scrnc_drccd_uhf_sto3g`: max asymmetry 0.0 (still exact, as before)
- `h2o_scrnc_rccd_uhf_sto3g`: max asymmetry 0.0 (was 6.653E-04)
- `h2o_scrnc_ccsd_uhf_sto3g`: max asymmetry 0.0 (still exact, as before)

Symmetrized values are sane (same order of magnitude/structure as the
pre-fix asymmetric pairs, consistent with a genuine average rather than
corruption) -- e.g. H2O DRCCD JJ(1,1)=0.79387968734,
JJ(1,3)=JJ(3,1)=0.30939310680.

**This closes out the entire scrnc JJ-asymmetry investigation**, this
time for real: both the H2 list-374 collision bug from 2026-08-14 (real,
fixed, committed) and the missing-Symmet2-call bug found and fixed today.

**Cleanup items also CLOSED (2026-08-15, committed `77aa9fc4`)**: the two
minor indexing bugs found today (`scrnc_driver.F:144` Ndim off-by-one,
`scrnc_form_lrspn_rcc_lists.F:61` Length beta-column typo -- both
described in detail above) are now fixed and committed. Verified by a
fresh full-pipeline rerun of `h2o_scrnc_drccd_uhf_sto3g`
(`/blue/ufhpc/perera/scrnc_test/h2o_scrnc_drccd_uhf_sto3g_idxfix/`): JJ
matrix bit-identical to the pre-fix baseline
(`JJ(1,1)=0.79387968734`, `JJ(1,3)=0.30939310680`), confirming both fixes
are numerically inert for closed-shell UHF as expected -- real effect
only for open-shell UHF references, not currently tested by any scrnc
case.

**Working tree status**: clean on the scrnc thread -- both `Symmet2`
calls (`35cdb853`) and both indexing-typo fixes (`77aa9fc4`) are
committed. The raw-trace instrumentation added earlier today in
`scrnc_form_mubar_s.F`/`scrnc_jj_rcc_lrspn.F` is guarded by
`_DEBUG_LVL0` (already the project's standing debug flag, see
`GNUmakefile.src`) and was left in as permanent diagnostics, also
committed as part of `0cddc3eb`/`35cdb853` -- not exploratory cruft to
clean up.

---

## SESSION 2026-08-15 evening -- dynamic (frequency-dependent) W(omega), CCSD path (Phase 1, WORKING, committed 83c8980d)

Goal shift: from static W=v+v*chi*v toward an actual GW-like quasiparticle
pipeline, per Ajith's direction ("mimic dynamic pol[arizability]").
`vcceh` already has a working dynamic-polarizability solver
(`secprp.F`, driven by a `frequency` file, FREQ/REALFREQ threaded into
`Lineqy`/`Plineqy`) that was never imported into scrnc -- this session
imported it into scrnc's CCSD path (`Rccd`/`Drccd` ring-CC path NOT yet
touched, still static-only).

**Key finding, secprp vs secprp2**: `secprp2.F` (the NMR spin-spin
coupling-tensor code, `JJSPNCUP` call -- direct ancestor of scrnc's own
"JJ_lrespn" naming) is static-only, `FREQ=0.0d0`/`FACT=TWO` hardcoded,
no omega loop at all. `secprp.F` (general dynamic polarizability) has
the real ITOP/FACT/+-omega machinery. scrnc's response-list
architecture (lists 373/374/375/376, common `/Freqinfo/`) was already
modeled on this exact secprp scheme -- `Freq`/`Realfreq` were already
being threaded into `Lineqy`/`Plineqy` calls in
`scrnc_lrspn_cc_main.F`, and `scrnc_form_lrspn_cc_lists.F` already
allocated all 4 lists -- but nobody ever varied `Freq` away from 0 or
did the +-omega double-solve-and-sum. That's what this session added.

**Implementation** (mirrors `secprp.F` almost line-for-line):
- `scrnc_lrspn_cc_main.F`: added `Itop = 1 if Freq=0-or-imaginary else
  2` and wrapped the existing per-perturbation solve loop in
  `Do Icount = 1, Itop`, negating `Freq` between passes. On the first
  pass of a 2-pass (+omega) run, the converged R-side amplitude is also
  saved to list 375 (secprp's own convention) before the second
  (-omega) pass overwrites list 373 in place.
- `scrnc_jj_cc_lrspn.F`: generalized the final assembly Ddot loop to
  the same `Do Icount = 1, Itop` structure, reading list 375 (+omega)
  on the first pass of a 2-pass run and list 373 (-omega, or the single
  static solve) otherwise, with `X = Fact*Ddot(...)` (`Fact=2` for
  `Itop=1`, `Fact=1` for `Itop=2`) instead of the old plain
  unweighted `Ddot`. The L-side Mbar/list-374 build (lines ~141-230,
  unaffected by frequency -- property operators like the perturbation
  Mbar don't depend on omega, only the solved amplitude does) is
  untouched.
- Two incidental bugs fixed while in this code: (1)
  `scrnc_form_lrspn_cc_lists.F:208`, the same alpha/beta
  copy-paste-off-by-one-column typo as the two fixed earlier today,
  except this one sizes the actual 373-376 response-amplitude storage
  (not just scratch padding) -- for open-shell UHF this could have
  under-allocated and corrupted memory, not just been numerically inert.
  (2) `scrnc_set_vars.F`: `Inquire(File="frequency",...)` vs
  `Open(...,File="Frequency",...)` -- case mismatch on a
  case-sensitive filesystem, meant ANY `frequency`-file run would crash
  immediately (`forrtl: file not found, unit 10`). Both fixed, both
  numerically inert for every test run so far (closed-shell UHF /
  never-before-exercised frequency path), same caveat as before.

**MAJOR FINDING, changes default static-path numbers**: rebuilding
after these changes (note: `make` alone did NOT relink `aces2/bin/xscrnc`
this time either, needed `make install` -- same trap as always, caught
via stale timestamp before wasting a test cycle) and rerunning the
static control case (`h2_scrnc_ccsd_uhf_sto3g`, no `frequency` file,
fresh full pipeline) gives `JJ(1,1)=JJ(4,4)=0.83435406798`,
`JJ(1,4)=JJ(4,1)=0.50949591159` -- DIFFERENT from the previously
committed/verified `0.80447994652`/`0.53936933731`. Isolated the cause:
the new `Fact=2.0` (introduced by mirroring secprp's `ITOP=1,FACT=TWO`
static branch) exactly doubles the *response-only* part of JJ (bare-v
part untouched) -- confirmed by hand: old response(1,1) =
0.80447994652 - (bare v, backed out) = 0.02987412147; new response(1,1)
= 0.05974824293 = exactly 2x. Still exactly symmetric either way (H2's
known special case). This means the OLD "verified" static CCSD-path
JJ was always self-consistent (symmetric, matched scrnc's own
independent Hbar-eigenvalue diagnostic) but a factor of 2 off from
vcceh's own established response-normalization convention in its
response term -- invisible to every prior check because those checks
only tested internal consistency (symmetry, no crashes), never the
absolute normalization against the vcceh convention. Matches exactly
what Ajith described: the EOM-effective-Hamiltonian route (secprp2/
CI-like, what scrnc's JJ was modeled on) is known to differ from
genuine linear response by factors like this; the quadratic/
Perera-Nooijen-extended correction (`IEOMPROP=2`, already wired via
`Calcquad1` in `scrnc_jj_cc_lrspn.F`) is the documented fix for the
gap, "small but can be significant."

**RESOLVED, 2026-08-15**: Ajith confirmed `Fact=2` on the static branch
is correct -- committed as `83c8980d`, superseding the old static
numbers everywhere.

**Validated, dynamic mechanism itself works correctly**: real-frequency
test cases (`h2_scrnc_ccsd_uhf_sto3g_freqtest`/`_freqtest2`, H2/STO-3G/
CCSD, `frequency` file format `No_freqs` / real-or-imaginary flag (0/1)
/ freq value per line, one value per line after the first two)
confirm:
- omega=0.001 a.u.: response(1,1)=0.05974830686 vs omega=0's
  0.05974824293 -- tiny, smooth, continuous shift as expected.
- omega=0.3 a.u.: response(1,1)=0.06611559494 -- clear, sane,
  non-divergent dispersion (~11% growth), still exactly symmetric
  (`Symmet2` a no-op both times, LR JJ = CI-like JJ).
- No crashes, no memory errors, in either case.
- One test-harness lesson (not a code bug): copying a completed test
  directory wholesale and rerunning only `xscrnc` on top of its
  leftover list files (instead of a fresh full pipeline from ZMAT+
  GENBAS) throws `@ACES_LIST_TOUCH: list already exists and is of
  different dimensions` -- confirms the existing "not idempotent,
  always rerun fresh" rule applies to xscrnc's own list files too, not
  just JOBARC-level quantities.

**Not yet done**:
- Ring-CC/DRCCD/RCCD path has NO frequency-shift capability yet
  (`scrnc_form_rcc_hbar_inv.F`'s Hbar^-1 build via plain `MINV` has no
  omega argument at all) -- CCSD path only so far, which is fine since
  that's "the real production target" per the earlier goal-framing
  entry, but ring-CC's RPA-vs-published-GW validation role means it'll
  need this eventually too.
- Imaginary-frequency axis not yet exercised (only tested REALFREQ=
  .TRUE. so far) -- the file-format/REALFREQ machinery claims to
  support it (mirrors vcceh's own IMAG_FREQ handling) but untested.
- No actual self-energy/quasiparticle-equation construction yet --
  this session only gets W(omega) itself working dynamically, which is
  the necessary input, not the full GW pipeline.
- Committed as `83c8980d` (four source files + this log).

---

## SESSION 2026-08-15 evening -- dynamic W(omega), ring-CC/RPA path (Phase 2, WORKING, committed)

Extended the same dynamic-response machinery to the ring-CC (RCCD/
DRCCD) path, per Ajith's direction ("do ring-CC first, so we have
something to verify our machinery" -- ring-CC is the one path with a
known-good external cross-check, both against ringcc/ringcl's own
independent Hbar eigenvalues and against the published GW=EOM-drCCD
identity).

**Architecture**: ring-CC builds Hbar explicitly as a matrix (unlike
CCSD's implicit Lineqy/Plineqy iterative solve) and inverts it once via
MINV, sharing the result (`Hbarinv`/`Iloc`) between the R-side
(`Scrnc_lrspn_rcc_main.F`) and L-side (`Scrnc_jj_rcc_lrspn.F`) -- and,
per the code's own comment, ring-CC's L-side genuinely needs Hbar^(-1)
too (unlike CCSD's L-side, which is a fixed Freq-independent operator).
So going dynamic here means building `(Hbar+Freq)^(-1)` instead of
`Hbar^(-1)`, once per +/-Freq pass, and sharing it with BOTH sides for
that pass -- both sides end up Freq-dependent, unlike CCSD.

**A real, non-trivial bug found and fixed (not a guess-and-check
sign issue)**: initial implementation (wrap the Hbarinv-build+R-call+
L-call in a Do Icount=1,Itop loop, negate Freq between passes, same
shape as the CCSD path) gave a response that *decreased* with
increasing Freq -- backwards from the CCSD path's correct increasing-
toward-resonance trend. Flipping the shift sign (Hbar-Freq vs
Hbar+Freq) changed the number but made it WORSE in the wrong direction,
which is a logical contradiction (summing +Freq and -Freq shifts should
be invariant to which pass is labeled which, by commutativity of
addition) -- proof this wasn't a sign-convention question at all but a
real bug.

**Isolated via direct diagnostic printing** (Hbar(1,1) before/after
shift, both passes): calling `Scrnc_rcc_hbar_inv` (the top-level ONE
that `Scrnc_driver.F` actually calls) twice in a row, even with
Freq=0.0 throughout (fully isolating away from anything Freq-related),
gave DIFFERENT raw Hbar(1,1) values each time
(1.40743653968732 -> 2.65424468244473 -> 3.90105282520215 on a third
call), with the EXACT SAME increment (1.24680814275741) added every
time -- a constant, deterministic, additive corruption, not noise.

**Root cause**: `Scrnc_rcc_hbar_inv.F` (the actual routine
`Scrnc_driver.F` calls -- easy to miss, named almost identically to
`Scrnc_form_rcc_hbar_inv.F`, which is the one this whole thread had
been reading/editing/instrumenting) is a thin WRAPPER:
```
Call Scrnc_form_rcc_hbar(Work,Memleft,Iuhf,Imult)
Call Scrnc_form_rcc_hbar_inv(Work,Memleft,Iuhf,Iloc)
```
`Scrnc_form_rcc_hbar` -- previously believed to be dead/orphaned code
(a grep for its call sites turned up nothing, see the 2026-08-15
"actual goal" entry above) -- IS in fact called, every time, via this
wrapper. The grep that missed it filtered matches with
`| command grep -v "_inv\|_INV"` to exclude irrelevant "_inv" hits, but
that filter matches on the whole line INCLUDING the filename, and
`scrnc_rcc_hbar_inv.F`'s own filename contains "_inv" -- silently
excluding the one file that actually mattered. Lesson: a second grep
pass "to clean up noise" can eat the real match if it filters on
whole-line text rather than being scoped to what's actually noise.

`Scrnc_form_rcc_hbar` destructively reads the raw AA/BB Hbar blocks
(lists 54/55) via `GETALL`, adds the one-particle Fock diagonal
(`F(ab)d(ij)-F(ij)d(ab)`, via `Scrnc_addoneh`), and writes the result
back via `PUTALL` to the SAME list it read from -- exactly the
"destructive read-modify-write" bug class already fixed once before
(2026-08-14, the original Hbar-double-build bug) in a *different*
routine. Calling the wrapper twice (needed for the dynamic Icount loop)
added the Fock diagonal a second time on the second call, corrupting
the "raw" Hbar that gets shifted+inverted -- which is exactly why
flipping the shift sign changed the number in the wrong way: both
sign choices were being applied to an already-corrupted second-pass
Hbar, so neither could give the right answer.

**Fix** (clean architectural split, not a patch): new
`Scrnc_rcc_hbar_read.F` -- pure GETLST-based read of the raw Hbar
blocks (RHF or UHF AA/BB/AABB/BBAA), no dressing, no PUTALL, safe to
call as often as needed since it never writes. `Scrnc_form_rcc_hbar_inv.F`
reworked to take the pristine Hbar as an INPUT array (`WHBAR`) plus a
separate OUTPUT array (`W`) -- copies `WHBAR`->`W` via `DCOPY`, applies
the Freq shift and MINV to `W` only, leaving `WHBAR` untouched.
`Scrnc_driver.F` now calls `Scrnc_form_rcc_hbar` (the dressing/PUTALL
step) exactly ONCE, then `Scrnc_rcc_hbar_read` once into a dedicated
pristine buffer (`Ihbar_bgn`, same `Ineed+Ipad` sizing as `Ihbarinv_bgn`
since `Scrnc_form_rcc_hbar` needs it as scratch too) -- both before the
dynamic Icount loop. Inside the loop, only
`Scrnc_form_rcc_hbar_inv(Work(Ihbar_bgn),Work(Ihbarinv_bgn),...)` gets
called per pass, which only copies+shifts+inverts, never re-touches the
underlying lists.

Also moved `Scrnc_form_lrspn_rcc_lists` (list creation) out of
`Scrnc_lrspn_rcc_main.F` and into `Scrnc_driver.F`, called once before
the Icount loop instead of once per pass -- same "don't repeat
list-touching operations across a dynamic loop" lesson, cheap insurance
even though it turned out not to be the actual bug this time.

List 376 (already allocated by `Scrnc_form_lrspn_rcc_lists.F`, unused
until now) is used as the L-side's +Freq companion to list 374, mirroring
375's role as the R-side's +Freq companion to 373. New
`Scrnc_jj_rcc_assemble.F` (split out of the old `Scrnc_jj_rcc_lrspn.F`,
which now only builds the per-pass L-side Mbar/374(+376)) runs once
after the Icount loop, combining {375,376} (+Freq pass) and {373,374}
(-Freq pass, or the single static solve) with the same Fact=2(static)/
1(dynamic) convention as the CCSD path's assembly.

**Verified** (H2/STO-3G/UHF/DRCCD, fresh full pipelines each time,
`/blue/ufhpc/perera/scrnc_test/h2_scrnc_drccd_finaltest_*`):
- Static (no `frequency` file): response(1,1)=0.03008849944 (before
  bare-v add) -- exactly reproduces the earlier post-Fact=2-fix
  baseline, confirming the read/dress/invert split changed nothing
  about the static answer.
- omega=0.1 a.u.: response(1,1)=0.03045810992 -- correct direction
  (increasing), small.
- omega=0.4 a.u.: response(1,1)=0.03666255361 -- clearly increasing,
  accelerating as omega approaches the known ~0.55 a.u. (15.09 eV)
  RPA excitation energy from below, exactly the expected dynamic-
  response shape. Exactly symmetric (`Symmet2` a no-op) at every
  frequency tested.
- CCSD path regression: h2_scrnc_ccsd_uhf_sto3g static case still
  gives the exact committed baseline (0.83435406798/0.50949591159) --
  confirms the shared `Scrnc_driver.F` changes didn't disturb the CCSD
  branch.

**Committed**: `scrnc_driver.F`, `scrnc_form_rcc_hbar_inv.F`,
`scrnc_lrspn_rcc_main.F`, `scrnc_jj_rcc_lrspn.F`, plus new files
`scrnc_rcc_hbar_read.F` and `scrnc_jj_rcc_assemble.F`.

**Not yet done**: imaginary-frequency axis still untested for either
path (file-format support exists, mechanism untested). No self-energy/
QP-equation construction yet on either path -- this and the CCSD-path
session both only get W(omega) itself working, which is the necessary
input, not a full GW pipeline. `scrnc_form_rcc_hbar.F`'s companion copy
in `aces2/rpa/scrnc_form_rcc_hbar.F` (found during the dead-code
search) not investigated -- may be a vestigial duplicate, may be
genuinely used by the native RPA module, not part of this thread.

---

## SESSION END, 2026-08-14 evening -- RESUME TOMORROW HERE

Stopping for the day. State to pick back up with, in one paragraph: H2's
list-374 collision bug (both L-side and CCSD R-side instances) is fixed
and committed, confirmed still correct. H2O's ring-CC path (DRCCD/RCCD)
is a SEPARATE, still-open bug -- JJ is genuinely asymmetric
(1.104E-04 DRCCD / 6.653E-04 RCCD max diff, both deterministic, confirmed
via independent reruns), while H2O's CCSD control path stays exactly
symmetric. Two candidate causes have been tested and RULED OUT today (see
below): the `scrnc_form_mubar_s.F` I020 sizing typo (A/B tested, zero
effect on this system) and the Hbar^-1 transpose identity (verified clean
to 1E-15). The next concrete step, not yet started: raw-value-dump trace
of `scrnc_form_mubar_s.F`'s `Mubar_s(Iside=2)` XGEMM/`MAIBAR` formation
(~lines 287-320) for the H2O worst-offending AO-pair indices
(9/20/38/41/49), the same method (not just checksums) that caught the
original H2 bug. Start there tomorrow -- don't re-derive from scratch,
and don't re-litigate the two ruled-out candidates above.

---

## RESUME POINT (as of 2026-08-14, latest -- CORRECTED, NOT RESOLVED)

**STATUS CORRECTION, 2026-08-14 evening.** The "RESOLVED" claim directly
below (from earlier the same day) is WRONG for H2O -- it was only ever
verified against H2/STO-3G's tiny 4x4 JJ matrix. When actually tested on
H2O/STO-3G (49x49 JJ matrix), **the ring-CC path (DRCCD/RCCD) is
genuinely asymmetric**, with the two H2-derived fixes (`b2443051`,
`82f3f2fe`) both already in place. The CCSD (non-ring) control path is
still exactly symmetric for H2O too, so this is confined to the same
ring-CC L-side architecture as before, but it is evidently a DIFFERENT or
ADDITIONAL bug from the list-374 collision -- that one is confirmed fixed
(H2 case still exact), this is new.

**Hard numbers (2026-08-14 evening), from `/blue/ufhpc/perera/scrnc_test/`:**
- `h2o_scrnc_ccsd_uhf_sto3g` (production CCSD-GW path): exactly symmetric,
  max|JJ(r,c)-JJ(c,r)| = 0.0. Still clean.
- `h2o_scrnc_drccd_uhf_sto3g` (ring-CC, orig run) and
  `h2o_scrnc_drccd_uhf_sto3g_rerun` (completely fresh pipeline rerun,
  same ZMAT/GENBAS, run independently end-to-end): **max asymmetry =
  1.103811E-04 at JJ(9,20)=-0.00496532118 vs JJ(20,9)=-0.00507570228**,
  bit-identical between the two independent runs -- deterministic bug,
  not a race/numerical-noise artifact. Next-worst: JJ(9,41) diff
  1.027E-04, JJ(9,49) diff 6.32E-05, JJ(33,38) diff 5.98E-05.
- `h2o_scrnc_rccd_uhf_sto3g` (ring-CC, regular RPA variant): **max
  asymmetry = 6.653099E-04 at JJ(9,41)=1.28137548208 vs
  JJ(41,9)=1.28204079195**. Same general pattern (row/col 9 recurring on
  the worse-offending pairs), larger magnitude than DRCCD.
- No obvious simple signature yet (the old H2 bug's "diagonal-diff exactly
  equals off-diagonal-diff" signature has NOT been checked for H2O yet --
  worth checking before assuming it's unrelated to the fixed collision
  bug's mechanism).

**Process failure worth noting for next time**: my own first attempt to
check H2O symmetry (same session, human collaborator flagged it was
wrong) used a regex parser that assumed one matrix column per
`COLUMN` header block. `scrnc`'s real output (verified 2026-08-14 evening
via `sed`/raw inspection) uses the standard ACES2 `output` subroutine
format -- **4 columns per header line**
(`COLUMN 1  COLUMN 2  COLUMN 3  COLUMN 4`), with up to 4 values per `ROW`
line. The broken parser silently matched almost nothing (11 of 49 column
blocks) and reported a false "exactly symmetric, diff=0.0" for ALL THREE
H2O tests including the CCSD control -- which is why the CCSD 0.0 result
above should be trusted (it independently reproduced after the parser fix)
but the original DRCCD/RCCD "symmetric" claims earlier this session were
never real. Also note: the file additionally contains a SEPARATE, smaller
"Response-only JJ matrix (before bare v added)" debug block (from the new
uncommitted debug print in `scrnc_jj_rcc_lrspn.F`/`scrnc_jj_cc_lrspn.F`,
see below) printed just before the real final "RPA JJ matrix" block --
grab the LAST matrix in the file (search for `RPA JJ matrix`), not the
first one that parses.

**Uncommitted work found sitting in the tree this session** (from an
earlier, un-logged session that hit a terminal freeze before writing
back here -- NOT written by Ajith by hand): `git diff` on
`aces2/scrnc/{scrnc_form_mubar_s,scrnc_form_rcc_hbar_inv,scrnc_jj_cc_lrspn,
scrnc_jj_rcc_lrspn}.F` shows:
1. A real fix in `scrnc_form_mubar_s.F`: `I020=I010+MAX(NT2DSZAA*NT2DISAA,
   NT2DSZAB*NT2DSZAB)` corrected to `...,NT2DSZAB*NT2DISAB)` -- a
   copy-paste typo (squared instead of size*dis, inconsistent with the
   correct AA term right above it). Numerically inert when
   `NT2DSZAB==NT2DISAB` (true for H2/STO-3G's small UHF case), so may
   matter for H2O where AB block dimensions can differ -- **not yet
   confirmed whether this is related to the new H2O asymmetry, worth
   checking directly** (e.g. temporarily reverting just this one line and
   rerunning to see if the H2O asymmetry changes).
2. New `_DEBUG_LVL0` checksums around the `MAIBAR` XGEMM calls in the same
   file.
3. In `scrnc_form_rcc_hbar_inv.F`: an independent `(Hbar^-1)^T ==
   (Hbar^T)^-1` identity cross-check (builds Hbar^T, inverts it
   separately, diffs against the transpose of the production Hbar^-1,
   prints `Xpmax`). Guarded by `_DEBUG_LVL0`, same guard as the
   pre-existing DGEEV diagnostic; did NOT print `Xpmax` in the H2O DRCCD
   run's output even though `_DEBUG_LVL0` is clearly on (the sibling
   "drCCD Hbar inverse matrix" print under the SAME `#ifdef` block did
   fire) -- **unexplained, worth checking** (possible early RETURN/branch
   skips past it, or a build staleness at the time that run happened).
4. New "Response-only JJ matrix (before bare v added)" debug print in both
   `scrnc_jj_cc_lrspn.F` and `scrnc_jj_rcc_lrspn.F` (this is the extra
   matrix block noted above that a naive "first matrix in the file" parse
   will grab by mistake).

All of this is STILL uncommitted (`git status` in
`~/Develop/ACES_PROGRAM` shows these 4 `.F`/`.f` pairs modified, not
staged). Do not commit yet -- Ajith has not confirmed the mubar_s fix or
the new diagnostics are correct/wanted, and the H2O asymmetry means the
debugging isn't finished regardless.

**Both "transformation" suspects ruled out (2026-08-14 evening, DONE):**
1. **I020 A/B test**: reverted just the `scrnc_form_mubar_s.F` I020 fix
   back to the buggy `NT2DSZAB*NT2DSZAB`, rebuilt (confirmed via md5sum,
   different binary), reran `h2o_scrnc_drccd_uhf_sto3g_i020ab` fresh from
   ZMAT+GENBAS. Result: **bit-identical** to the fixed build -- same
   1.103811E-04 max asymmetry, same worst pairs, same values to the last
   digit. `NT2DSZAB==NT2DISAB` for this H2O/STO-3G/UHF system too (same
   as H2), so the AB-block-squareness this fix addresses never actually
   triggers here. Fix is still correct/worth keeping for other systems,
   but it is NOT the cause of the H2O asymmetry. Reverted the test back
   to the fixed source immediately after (md5sum-verified restored,
   `f69468720b0b61fe32af2208b8a6674f`) -- tree is back in the
   I020-fixed state.
2. **`(Hbar^-1)^T == (Hbar^T)^-1` identity check**: DOES fire (my first
   grep for it used the wrong search string, "Xpmax" -- that's the
   Fortran variable name, not the printed text; the real printed line is
   `max|(Hbarinv)^T - (Hbar^T)^-1| = `). Found at
   `h2o_scrnc_drccd_uhf_sto3g_rerun/xscrnc.out:3418`:
   **`0.111022D-15`** -- clean to machine epsilon. The Hbar^-1 transpose
   shortcut is NOT the bug.

**Immediate next step**: with both current "transformation" suspects
cleared, narrow into the actual `Mubar_s(Iside=2)` formation math itself
in `scrnc_form_mubar_s.F` (the XGEMM calls building `MAIBAR` from T2 and
`MAI`, lines ~287-320) -- trace with RAW VALUE dumps for a couple of
AO-pair indices spanning the worst-offending H2O pairs (9/20/38/41/49),
not just checksums (checksums alone missed the original H2 bug's
sign/collision signature; this is the same lesson, method already
validated). Compare Iside=1 vs Iside=2 intermediate values for the same
underlying quantity, the same cross-checkpoint approach that worked
before.

---

## OLDER, INCORRECT "RESOLVED" entry (kept for its real content on the
## H2 bug trace, which IS still believed correct -- just don't trust the
## top-line "done" framing below)

**RESOLVED, 2026-08-14. Root cause found and fixed -- TWO independent
instances of the same collision bug, both fixed. JJ_lrespn is now exactly
symmetric for DRCCD, RCCD, AND the CCSD control case.** See "Step 4
result / ROOT CAUSE" and "Second instance / CCSD path" below for the full
trace. Steps 1-3 (Hbar proven correct, twice over) turned out to be red
herrings for the actual bug -- it was never in the CC/Hbar physics at
all, it was an ACES2 list-storage bookkeeping collision, independently
present in both the ring-CC L-side (`scrnc_jj_rcc_lrspn.F`) and the CCSD
R-side (`scrnc_lrspn_cc_main.F`). The CCSD control case's original
"symmetric to 11 digits" claim was accidental -- it had BOTH bugs
(L-side in `scrnc_jj_cc_lrspn.F`, R-side in `scrnc_lrspn_cc_main.F`)
simultaneously, and they happened to cancel for this specific H2 system.
Fixing only the L-side first (prompted by Ajith noting the CCSD path
"has not been tested thoroughly") broke that accidental cancellation and
exposed the R-side bug too -- both are now fixed and verified together.
**This H2-only claim is confirmed still true (re-checked with the corrected
parser); it's the generalization to H2O that fails, see correction above.**

**Immediate next step — DO NOT diverge from this into generic code-reading.**
The bug is isolated to `aces2/scrnc/scrnc_jj_rcc_lrspn.F`'s `Iside=2`
(L-side) branch (see "Where the bug is isolated to" below), but the METHOD
for chasing it must stay consistent with the debugging path already built
and validated this project (the "3-route" comparison: native RPA / ringcc
/ CCSD) — not ad-hoc source reading in isolation. Concretely, in order:

1. Instrument the L-side (`Iside=2`) with checksums, but checkpoint them
   **against the R-side's own already-correct intermediates** for the same
   raw integrals (R-side is proven correct via the route-1 RPA-sourced-T2
   test) — not in isolation.
2. Wire in the existing orphaned dead-code cross-check
   (`ringcl/rcl_check_rcchbar.F` + `rcl_drvrpa.F`, `#ifdef _DEBUG_LVL0`,
   zero current callers) to get ringcc's own independent Hbar
   diagonalization/eigenvalues for this system. This is step 1 of the
   original agreed two-step plan ("verify scrnc's Hbar is genuinely
   ringcc's Hbar") that was never executed — cheap, already written, just
   needs a call site added in `rcl_driver.F`.
3. Run the already-built, never-yet-executed test case
   `h2_rpa_verify_uhf_sto3g` (`excite=EOM(Sf)-DRCCD`) to get ringcc's own
   independently-reported EOM(Sf)-DRCCD excitation energy, and compare
   against scrnc's Hbar eigenvalues from step 2 — confirms whether scrnc's
   shared Hbar itself is trustworthy before blaming the L-side's
   Mubar/transpose logic specifically.
4. Only if Hbar checks out clean, narrow into `Mubar(Iside=2)` +
   `Hbar^(-1)^T` transpose as the remaining suspect.

**On every resume, re-read this numbered list before proposing a plan.**
Do not re-derive a "go read the suspect file" plan from scratch — it skips
the route-comparison methodology and has already caused one correction
this project (see feedback memory `feedback_stay_on_agreed_debug_path`).

**Working tree state (corrected 2026-08-14, post-freeze resume)**: the
Hbar single-build refactor + RPA_T2VECS bypass fixes listed above as
"uncommitted" were actually committed just before the terminal froze —
`2e8d79f0` ("scrnc/ringcc: Hbar single-build refactor + RPA_T2VECS bypass
fixes for GW debug"). `git status` now shows those files clean. Working
tree is 8 local commits ahead of `origin/master`, still not pushed. Only
stray leftovers from the freeze: a `.scrnc_driver.F.swp` (untracked) and a
deleted `.asv_ctl.c.swp` — both harmless vim swap noise, not touched.
`aces3/config.log`/`config.status` also show as modified (unrelated build
byproducts, not part of this thread).

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

**Step 1 result (2026-08-14, post-freeze resume, DONE)**: added matched
`@CHECKSUM` instrumentation (label + `Npert_count`) to both routines'
`Iside=1`-stage `Mbar_00`/`Mubar_s` intermediates — `scrnc_lrspn_rcc_main.F`
(R-side, proven correct) and `scrnc_jj_rcc_lrspn.F`'s own internal
`Iside=1` stage (computed before its `Iside=2` stage and the
`Hbar^(-1)^T` transpose). Rebuilt (`make install` in `aces2/scrnc/` with
`intel/2025.1.0` loaded — note: `which xscrnc` resolves to a STALE separate
binary at `~/Develop/bin/xscrnc`, md5-different, not a symlink; had to
invoke `~/Develop/ACES_PROGRAM/aces2/bin/xscrnc` by full path to test the
real rebuild) and reran `h2_scrnc_drccd_uhf_sto3g`. Confirmed the
documented asymmetry signature is unchanged (JJ(1,1) vs JJ(4,4) diff =
JJ(1,4) vs JJ(4,1) diff = 0.00448453013, both routines' checksums for all
4 `Npert_count` values bit-identical between R-side and L-routine.

**Conclusion**: the shared `Iside=1`-stage machinery (raw integral prep,
`Scrnc_form_mubar_00`, `Scrnc_form_mubar_s` for `Iside=1`) is NOT the bug
source — ruled out with a positive, decisive checksum match. This further
narrows suspicion onto exactly what step 4 already named: the L-routine's
own `Iside=2`-stage `Scrnc_form_mubar_s` call (using `Mbar_00-M_expect`
and the `Iside=1`-derived buffers as input) and/or the `Hbar^(-1)^T`
transpose block. Per the agreed plan, steps 2-3 (independent ringcc Hbar
eigenvalue cross-check) should still run before diving into that code,
since they test whether scrnc's *shared* Hbar itself is trustworthy
(orthogonal to this L/R intermediate check).

Instrumentation is uncommitted in `aces2/scrnc/scrnc_lrspn_rcc_main.F` and
`aces2/scrnc/scrnc_jj_rcc_lrspn.F` (checksums only, guarded by
`_DEBUG_LVL0`, no logic changes) — kept in the tree, not reverted, since
it's reusable for future checkpoints on this bug.

**Step 2/3 result (2026-08-14, DONE)**: `RCL_DRIVER.F`'s own call site
(where the plan originally assumed the check belonged) turned out to
never fire for this test — `RCL_DRIVER` is the lambda-equation solver,
gated on `Solve_4lambda=.True.`, which an EOM-only run like this never
sets. The pathway actually exercised by `excite=EOM(Sf)-DRCCD` is
`rcl_hbar.F`, which builds Hbar lists 54-59 and returns with nothing ever
diagonalizing them — confirming the log's "never executed" description
of this self-check was accurate for the reason originally guessed
(no caller) but the *fix* location was wrong; reverted the `RCL_DRIVER.F`
call site (added defensively, never verified, would have been an
unverified corruption risk if that path is ever exercised) and wired the
call into `rcl_hbar.F` instead, right after `rcl_set_hbar2zero`.

Caught a real repeat of the exact same double-build bug class before it
shipped: `Rcl_check_rcchbar` → `Rcl_heffxc` destructively overwrites the
same same-spin list(s) it reads from (54/55 UHF, 56 RHF) by adding the
Fock diagonal on top — needed by `Rcl_drvrpa`'s own diagonalization, but
would silently corrupt those lists for whatever downstream production EOM
step consumes `rcl_hbar`'s raw output next, same failure mode as the
scrnc Hbar bug. Added explicit save/restore of those lists around the
diagnostic call (`GETALL` before, `PUTALL` after) in `rcl_hbar.F`.
Verified side-effect-free: reran the full fresh pipeline
(`xjoda,xvmol,xvmol2ja,xvprops,xvscf,xvtran,xintprc,xringcc,xringcl`) on
`h2_rpa_verify_uhf_sto3g` and confirmed the raw/final `Wprime`/list-54-57
checksums are bit-identical to the uninstrumented baseline
(`0.170975410565E+00` raw, `0.160628396930E+00` final on all of
54/55/56/57) — the check's own list-54/55 mutation is fully undone.

**Result: ringcc's own independent EOM(Sf)-DRCCD Hbar diagonalization
gives real eigenvalues 42.66923610 eV and 33.92739023 eV** (imaginary
parts both exactly zero, as expected for a real, non-defective Hbar in
this 2-state system). This is ringcc's own number, computed via a
completely separate code path (`Rcl_heffxc`+`Rcl_drvrpa`, full `DGEEV`
non-symmetric diagonalization) from anything scrnc touches.

**Not yet done, needs a physics call before proceeding**: comparing this
against scrnc's own Hbar. scrnc never reports raw eigenvalues — only the
JJ=W=v+v·χ·v matrix — so "confirms whether scrnc's shared Hbar is
trustworthy" (per step 3's original framing) needs Ajith to say what
scrnc-side number should match 42.669/33.927 eV, or whether these should
instead be cross-checked against the native RPA module (`xrpa`)'s own
independent output for the same system as a third, orthogonal source.
Test dir: `/blue/ufhpc/perera/scrnc_test/h2_rpa_verify_uhf_sto3g/`
(rebuilt from scratch, all 9 pipeline steps' `.out` files present).

**2026-08-14 follow-up, `_RCC` flag clarified**: Ajith asked to add
Hbar-diagonalization code to scrnc itself (see new TODO below), and
separately to switch to "Sf-RCCD, not DRCC" via recompiling `rpa` with
its CPP directive off. Tried `_RCC` off first — but that selects the
`#else` branch in `do_uhf_rpa.F`/`mkb_uhf.F` (prints "dr-CCD Correlation
energy", list-4+ISPIN numbering), i.e. DRCCD, the opposite of what was
asked. Ajith clarified: **"regular RPA" = ring-CCD (RCCD), keeps BOTH
Coulomb (J) AND exchange (K) integrals** — that's the `#ifdef _RCC`
branch ("r-CCD" label, list-18+ISPIN numbering), not `#else`. **Reverted
`rpa` back to `_RCC` ON** (`make ppclean && make install` in
`aces2/rpa/`, no override) — matches the tree's actual global default in
`Makefiles/GNUmakefile.src` (never touched), confirmed via the same
label/list-number check. Net effect: no change from before this
excursion, `xrpa` is back to its normal build. Direct-vs-exchange
terminology note for future reference: DRCCD/dRPA = direct-only
(Coulomb/J), RCCD/regular-RPA = full (Coulomb+exchange, J+K) — `_RCC`
selects the latter when defined.

**Step 2/3/scrnc-Hbar-diag result (2026-08-14, DONE)**: added independent
Hbar diagonalization directly to scrnc, in
`scrnc_form_rcc_hbar_inv.F`, right before the existing `CALL MINV(...)`
(which destroys its input in-place) -- snapshots Hbar into scratch
(`DCOPY`) first so the diagnostic (`DGEEV`, same non-symmetric
diagonalization ringcl's own check uses) can never perturb the
production `Hbar^(-1)` build MINV computes right after. Guarded by
`_DEBUG_LVL0` like everything else here. Rebuilt cleanly (`make install`
in `aces2/scrnc/`).

Built a new RCCD (regular RPA, both Coulomb+exchange kept) scrnc test,
`/blue/ufhpc/perera/scrnc_test/h2_scrnc_rccd_uhf_sto3g/` (H2/STO-3G/UHF,
`CALC=RCCD`, same 10-step pipeline including `xscrnc`) -- **scrnc's own
new diagnostic gives real eigenvalues 25.26716566 eV and 15.09016923 eV**
(imaginary parts zero). Also built the RCCD counterpart of the ringcc
verify test, `h2_rpa_verify_rccd_uhf_sto3g/` (`excite=EOM(Sf)-RCCD`,
`eomref=RCCD` -- the "switch to Sf-RCCD, not DRCCD" Ajith asked for) --
**ringcc's independent diagonalization gives the exact same
25.26716566 eV and 15.09016923 eV.** Bit-identical between the two
completely separate code paths (scrnc's own `Scrnc_form_rcc_hbar`+`DGEEV`
vs ringcl's `Rcl_heffxc`+`Rcl_drvrpa`+`DGEEV`).

Also confirmed the JJ-matrix asymmetry reproduces in RCCD too, same
exactly-equal-offset signature as DRCCD (diagonal diff = off-diagonal
diff = `0.06020595989` here, vs `0.00448453013` for DRCCD) -- physics-
method-independent, consistent with an architectural bug in scrnc's own
L-side code rather than anything CC-method-specific.

**`aces2/rpa`'s `_RCC` flag**: tried turning it off per Ajith's initial
phrasing, got DRCCD instead of the intended RCCD (see follow-up above) --
reverted to `_RCC` ON (`make ppclean && make install`, no override),
matching the tree's normal default. No net change to `rpa/` from before
this thread. `_RCC` is scoped to `rpa/` only (grep-confirmed: no other
module's `.F` files reference it), so the scoped `DEFINES_EXTRA=`
command-line override used during the experiment never touched the
shared `Makefiles/GNUmakefile.src` or any other module's build.

**Step 4 result / ROOT CAUSE (2026-08-14, RESOLVED)**: Ajith confirmed
Hbar being exactly symmetric for H2/CCSD-derived Hbar is a known special
case unique to H2 -- not itself suspicious -- which redirected the
search away from the `Hbar^(-1)^T` transpose (mathematically a no-op on
a symmetric matrix, so it could not have been the bug regardless).

Traced the actual mechanism by dumping raw vector values (not just
checksums) at each stage of `scrnc_jj_rcc_lrspn.F`'s `Iside=2` handling,
for the 4-AO-pair H2/STO-3G system:
1. Raw `Mu(482)` (pre-`Hbar`-multiply): correctly antisymmetric between
   AO-pair 1 and 4 (`+0.113.../-0.113...`).
2. `Work(Iwrk)` immediately post-`Xgemm`: still correctly antisymmetric
   (`+0.1218.../-0.1218...`) -- confirms the `Hbar^(-1)^T` multiply
   itself was never the problem.
3. List 374's record 1, read back immediately after being written: still
   correct at that instant.
4. The FINAL value of list 374 record 1, read later in the JJ-assembly
   `Ddot` loop: **wrong** -- came back as AO-pair 4's value, not AO-pair
   1's own.

Root cause: the code wrote the `Hbar^(-1)^T`-multiplied result to a
"scratch" slot -- **record 1 of list 374** -- via
`Putlst(Work(Iwrk),1,1,0,Irrepx,374)`, then immediately relocated it to
its real per-perturbation slot via
`Getlst(...,Ioffset=1,...)`/`Putlst(...,Npert_count,...)`. But record 1
of list 374 is ALSO the permanent final-storage slot for `Npert_count=1`
itself (since perturbations are numbered 1..Nao_pair and stored at
matching record numbers). So every iteration after the first
(`Npert_count=2,3,4,...`) reused that same record 1 as scratch, and each
one silently overwrote whatever `Npert_count=1`'s finished value had
been left there. By the time the JJ-assembly loop read record 1 back, it
held the LAST perturbation's value instead of the first's -- explaining
both the "additive offset tied to atom-2's AO-pair index" signature
(record `Nao_pair`'s value leaking into record 1's slot) and why fixing
the earlier, unrelated Hbar double-build bug made the asymmetry *worse*
(that fix changed what ended up in the last-written record, changing the
leaked value, not removing the leak).

**Fix** (`aces2/scrnc/scrnc_jj_rcc_lrspn.F`): write the `Xgemm` result
directly to its real record --
`Call Putlst(Work(Iwrk),Npert_count,1,1,Irrepx,374)` -- removing the
scratch `Putlst(...,1,...)` and the subsequent
`Getlst(record 1)`/`Putlst(record Npert_count)` relocation pair entirely
(`Npert_count` is already known at that point in the loop, so the
round-trip was never necessary). Also removed the now-unused `Ioffset`
variable.

**Verified**: rebuilt (`make install` in `aces2/scrnc/`), reran full
fresh pipelines for both control cases from clean directories:
- `h2_scrnc_rccd_uhf_sto3g`: `JJ(1,1)=JJ(4,4)=0.80470880500`,
  `JJ(1,4)=JJ(4,1)=0.53913978308` -- exactly symmetric.
- `h2_scrnc_drccd_uhf_sto3g`: `JJ(1,1)=JJ(4,4)=0.78965007478`,
  `JJ(1,4)=JJ(4,1)=0.55419851331` -- exactly symmetric.

scrnc's own Hbar-eigenvalue diagnostic (added earlier this session,
`scrnc_form_rcc_hbar_inv.F`) still reproduces ringcc's independent
25.26716566/15.09016923 eV exactly -- confirming this fix didn't touch
(and didn't need to touch) anything Hbar-related, consistent with the
root cause being pure list-bookkeeping, not physics.

Left in `scrnc_jj_rcc_lrspn.F` as ongoing diagnostics (kept, not
exploratory cruft): the "R-side/L-routine Npert_count" cross-checkpoint
checksums from step 1, and new "373 raw, Iao_pair=N" / "374 raw,
Jao_pair=N" raw-vector prints (these are what caught the sign/collision
bug -- checksums alone, being sum/sumsq, don't preserve enough
sign/structure information to have caught this one directly).

**Second instance / CCSD path (2026-08-14, RESOLVED)**: Ajith flagged
that the CCSD path "has not been tested thoroughly," so rather than trust
the proven-symmetric-for-H2 claim, applied the same defensive fix to
`scrnc_jj_cc_lrspn.F`'s own scratch-record-1-of-374 pattern (identical to
the ring-CC bug -- same `Putlst(Work(I000),1,1,0,Irrepx,374)` then
`Getlst(record 1)`/`Putlst(record Npert_count)` relocation). That alone
**broke** the CCSD control case's symmetry
(`JJ(1,1)=0.74473239934` vs `JJ(4,4)=0.80447994652`, asymmetric) --
meaning the "proven symmetric to 11 digits" result had been masking a
SECOND, independent collision bug on the R-side that was accidentally
canceling the L-side one for this specific H2 geometry (same "one bug
partially cancels another" pattern already seen once for the earlier
Hbar double-build fix).

Traced the R-side bug to `scrnc_lrspn_cc_main.F` (CCSD's R-side
companion to `scrnc_jj_cc_lrspn.F`, analogous to
`scrnc_lrspn_rcc_main.F` on the ring-CC side). Its `Iside=1` branch does
the same `Putlst(Work(I000),1,1,0,Irrepx,373)` scratch write, but the
"solve" in between isn't a simple matrix multiply like ring-CC's Hbar^-1
-- it's a genuine iterative Davidson-type linear-equation solve via
`Lineqy` (`vcceh/lineqy.F`) or `Plineqy` (`vcceh/plineqy.f`, in-core
variant), general CCSD response-property infrastructure shared well
beyond scrnc. Read both: `Lineqy`/`Plineqy` read their RHS from, and
write their converged solution back to, list `372+Iside` record
`Ioffset` -- but `Ioffset` is a fully generic pass-through parameter to
both (confirmed via `vcceh/lineqy.F:107,398` and
`vcceh/plineqy.f:89,328`) -- nothing requires it to be a fixed scratch
record. `scrnc_lrspn_cc_main.F` just happened to always pass
`Ioffset=1`, recreating the exact same record-1-collides-with-
Npert_count=1 problem, one level removed through the solver.

**Fix** (`aces2/scrnc/scrnc_lrspn_cc_main.F`, scrnc-local only --
`vcceh/lineqy.F`/`plineqy.f` untouched, so no blast radius beyond scrnc):
write the initial RHS directly to record `Npert_count` instead of record
1 (both the `Iside=1`→373 and dead `Iside=2`→374 branches, for
consistency), and set `Ioffset=Npert_count` (both the `Incore`/`Plineqy`
and non-`Incore`/`Lineqy` branches) instead of hardcoding `1`. Removed
the now-bug-causing relocate `Putlst`, but kept a plain `Getlst` after
the solve returns -- `Lineqy`/`Plineqy` write the converged solution into
their OWN scratch array (a different memory region than this routine's
`Work(I000)`) and to disk at record `Ioffset`, so a re-fetch is still
needed to get the result back into `Work(I000)` for the existing debug
checksum print.

**Verified**: rebuilt, reran the CCSD control case
(`h2_scrnc_ccsd_uhf_sto3g`, fresh full pipeline
`xjoda,xvmol,xvmol2ja,xvprops,xvscf,xvtran,xintprc,xvcc,xlambda,xscrnc`)
from a clean directory: `JJ(1,1)=JJ(4,4)=0.80447994652`,
`JJ(1,4)=JJ(4,1)=0.53936933731` -- exactly symmetric again, with BOTH
fixes (L-side + R-side) in place together.

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
  - `h2_rpa_verify_rccd_uhf_sto3g/` (RCCD counterpart of the above,
    ringcc's own Hbar eigenvalues for regular RPA)
  - `h2_scrnc_rccd_uhf_sto3g/` (RCCD counterpart of the drccd control
    case, scrnc's own Hbar eigenvalues + JJ matrix for regular RPA)
  - 11 further scrnc_test dirs from earlier sessions, see memory for detail
    if needed.

**PATH trap, bigger than previously documented (2026-08-14)**: `which
xjoda`/`which xscrnc`/`which xaces2` etc. all resolve to a STALE, separate
`~/Develop/bin/` tree by default (not a symlink, genuinely different
binaries, confirmed via `md5sum`) -- `aces2/bin` is NOT on `$PATH` by
default in a fresh shell. Always
`export PATH=/home/perera/Develop/ACES_PROGRAM/aces2/bin:$PATH` (and
`module load intel/2025.1.0`, same shell command, not piped -- piping
`module load | tail` runs it in a subshell and silently loses the PATH
change) before testing any rebuilt module, or you're silently running old
code. `xaces2`'s automatic pipeline builder also doesn't know to route
`EOM(Sf)-DRCCD`/`EOM(Sf)-RCCD` through the manual ring-CC pipeline -- it
sends these through the standard `xvcc`/EOM-CC machinery instead, which
can crash on H2/STO-3G edge cases unrelated to anything in this thread
(hit an integer-divide-by-zero in `xvcc`'s `e4s`/`t2t1aa1` this way).
Always build these test pipelines manually:
`xjoda,xvmol,xvmol2ja,xvprops,xvscf,xvtran,xintprc,xringcc,xringcl[,xscrnc]`.
New/empty test dirs also need `GENBAS` copied in by hand (3.5MB reference
basis file, not tracked in git) -- copy from any sibling `scrnc_test/*`
dir. Rerunning a single pipeline step (e.g. just `xringcl`) on top of an
already-completed run's leftover `JOBARC` is NOT idempotent for Hbar-type
quantities -- always clear the directory back to `ZMAT`+`GENBAS` and
rerun the full pipeline fresh rather than resuming partway.

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
