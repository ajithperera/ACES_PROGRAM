# scrnc / GW manuscript debug log

Plain-text resume document, kept in the repo so it survives terminal
freezes / HPG outages and gets backed up on every `git push`. Read the
**RESUME POINT** section first on any restart.

## SESSION 2026-08-23 (later) -- Ajith's own scope call: RCC/DRCC RHF/UHF
response-gap thread SET ASIDE, project focus moves to CCSD. RESUME HERE
if scrnc/RCC/DRCC work resumes; this is now the authoritative top entry.

Ajith reviewed the Nmult=1 re-close (entry immediately below) and made a
final scope decision rather than continue chasing the residual:
**RCC/DRCC's role in this project is comparison-with-literature and
debugging, not the main line of work -- CCSD is.** UHF RCC/DRCC is
already trusted (multiple independent consistency checks in this log
confirm it). The RHF/UHF response-property identity gap (1.88x H2/
STO-3G, up to ~12.5% H2O/STO-3G) is left OPEN AND UNPURSUED, not because
it's understood to be correct physics, but because closing it properly
is now believed to need a genuinely DIFFERENT kind of code for the RHF
identity (not further patching of the current dressing/Hbar/Nmult
pipeline, all of which are now confirmed correct as-is via 8ec738fa/
eea48bcc). **No further autonomous work on this thread unless/until
Ajith explicitly reopens it or commissions the separate RHF-identity
code.** Nothing to fix, nothing pending -- working tree for scrnc/
aces2/scrnc left exactly as eea48bcc committed it.

---

## SESSION 2026-08-23 -- Nmult=1 (singlet-only) decision REVISITED per
Ajith's explicit authorization ("you can revisit... need to wrap this
and move forward"), RE-CONFIRMED as correct on the NEW (post-L-side-fix)
baseline via a fresh empirical test, not just re-argued theoretically:
no single Imult=1/Imult=2 recombination weight reconciles RHF with UHF
-- the weight that would be needed varies from 0.32 to 0.70 across just
three tested diagonal elements (2 systems). **FORMALLY CLOSED as
understood, not-a-bug, not further pursued -- the 1.88x (H2/STO-3G) and
up-to-12.5% (H2O/STO-3G) residuals remain OPEN but are now definitively
NOT a singlet/triplet-decoupling/Nmult problem, on this or any prior
baseline. RESUME HERE if picking the residual back up -- next lead is
elsewhere, see bottom of this entry.**

**Re-derived the Wigner-Eckart argument from first principles, checking
specifically whether it depends on anything the L-side mubar fix
(`67a4ed6b`/`8ec738fa`) touched.** The argument: the perturbation this
whole pipeline builds the response to is a one-body, spin-summed,
density-type operator `V = Sum_pq v_pq (a+_pa a_qa + a+_pb a_qb)` -- a
spin-scalar (S=0) under spin rotations. For a closed-shell (S=0)
reference, Wigner-Eckart forbids ANY matrix element of a spin-scalar
operator between the S=0 ground state and an S=1 (triplet-coupled)
excited configuration, full stop -- this is an exact statement about
the OPERATOR's symmetry, not about how the T-amplitudes that go into
building Hbar or dressing the property vector happen to be computed.
Checked the actual dependency chain to make sure the fix couldn't have
changed this premise anyway: `Scrnc_form_mubar_s.F`'s RHF-only dressing
branch (the thing the L-side fix touched) is gated behind
`IF(IUHF.EQ.0)` at line 298, and the RCCD-specific `(list37-list34)/2`
fix sits INSIDE that branch (line 341) -- entirely separate machinery
from `Scrnc_driver.F`'s `Nmult`/`Imult` channel-selection loop and
`Scrnc_jj_rcc_assemble.F`'s `Wt`-weighted recombination, which is what
actually decides whether the triplet channel contributes at all. The
mubar fix changes the MAGNITUDE of what the (already triplet-excluded)
singlet-only channel computes; it has no mechanism to change whether a
second channel's contribution is physically zero. **Conclusion: the
L-side fix could not have invalidated the Wigner-Eckart premise, by
construction -- but per Ajith's mandate this needed empirical
reconfirmation, not just re-reading the same argument, so did the test
anyway.**

**Empirical test, per the project's "test, don't theorize" doctrine.**
Temporarily edited `Scrnc_driver.F` (`Nmult=2` for `Iuhf.EQ.0 .AND.
Rccd`, clearly marked "TEMP DIAGNOSTIC 2026-08-23, NOT FOR SHIPPING")
and `Scrnc_jj_rcc_assemble.F` (added a `_DEBUG_LVL0`-gated print of the
RAW, unweighted per-channel (`Imult=1`/singlet, `Imult=2`/triplet)
contribution to every diagonal JJ element, before the existing
0.25/0.75 `Wt` combination is applied) -- the existing two-channel
machinery from the 2026-08-19 era is still fully intact in the source
(dead code only because `Nmult` is hardcoded to 1), so no new physics
code was written, just instrumentation. Rebuilt `xscrnc`
(`ccfd4721547d193a14f0aeb57aa50575` after the second, all-diagonal-
elements revision of the print), md5-confirmed the binary actually
changed and contains the new debug strings. Ran four fresh manual
pipelines (`xjoda,xvmol,xvmol2ja,xvscf,xvtran,xintprc,xringcc,xringcl,
xscrnc`, no `xvprops`, per this file's 2026-08-20 correction) via
sbatch: `/blue/ufhpc/perera/scrnc_test/h2_scrnc_rccd_{rhf,uhf}_sto3g_
nmultdiag_0823` and `h2o_scrnc_rccd_{rhf,uhf}_sto3g_nmultdiag_0823`
(plus a `_nmultdiag2_0823` pair for the broadened per-element print),
all exit 0.

**Result: for every element tested, solved for the weight `w` that
WOULD make `Xraw(singlet) + w*Xraw(triplet) = UHF's true response-only
value` (computed from UHF's long-established, fix-independent full-JJ
reference values minus the bare-v piece, which is confirmed
bit-identical between RHF and UHF to 11 digits, both in this run and
every prior session's):**
```
                          Xraw(Imult=1)   Xraw(Imult=2)   UHF target   fitted w
H2/STO-3G   JJ(1,1)=(4,4) 0.032051820282  0.086977276510  0.06020595989   0.324
H2O/STO-3G  JJ(1,1)       0.049479347675  0.101942615519  0.09049168785   0.402
H2O/STO-3G  JJ(4,4)       0.003496540181  0.004931653338  0.00694352740   0.699
```
(UHF targets: H2's `0.06020595989` is this thread's long-standing
established response-only value; H2O's come from `UHF full JJ(1,1)=
0.86509751290`/`JJ(4,4)=0.03016663846` -- both established independently
of any RHF-side fix in this thread, since UHF's own construction never
enters `Scrnc_form_mubar_s.F`'s `IUHF.EQ.0` branch or `Scrnc_driver.F`'s
`Nmult` logic at all -- minus this run's own freshly-measured bare-v
diagonal, back-solved from `RHF_full - RHF_response(Nmult=2,weighted)`
and cross-checked against the raw AO integral print for element (1,1)
to 9 digits.)

**This is decisive, and it's the same verdict as every prior attempt in
this thread to fit a weight, just with new numbers because the baseline
moved: no single weight works.** `0.324`, `0.402`, and `0.699` span more
than a factor of 2 across just 3 diagonal elements from 2 systems -- if
a real triplet-channel admixture were the missing physics, the SAME
weight (a property of the perturbation operator's symmetry, not of
which AO pair you're looking at) would have to work everywhere. It
doesn't, on this baseline any more than it did on the old one (which
gave equally scattered, equally non-physical fits: 0.406, then
`Wt_S~2.2,Wt_T~-1.2`, see the 2026-08-19 entries below). **Both the
theoretical re-derivation and the fresh empirical retest agree: the
Nmult=1/singlet-only architecture is still correct, is NOT what the
L-side fix put in question, and this specific hypothesis is closed for
good -- not to be re-litigated again without a genuinely new reason to
doubt the operator symmetry argument itself (as opposed to "the residual
didn't go away").**

**Verdict for Ajith's "wrap this and move forward" mandate: (B).** The
Nmult=1 architecture is confirmed correct, not changed. The 1.88x
(H2/STO-3G) and up-to-12.5% (H2O/STO-3G) response-ratio residuals are
NOT fixed by this round and are not expected to be fixable by any
Nmult/weight change -- formally closing out the "is it a singlet/triplet
decoupling problem" question as understood, not-a-bug, exhausted from
both the theory side and the empirical side. The residual gap ITSELF
remains open as a separate problem for whoever picks it up next -- per
fork #5's own list, the concrete remaining candidates are (a) the still-
uncharacterized R-side (Iside=1) list37/34 vs list61/63 discrepancy
(traced but ruled a non-driver by fork #5's own substitution test), or
(b) something entirely un-isolated yet, given H2O's JJ(4,4) sector (the
"exactly symmetry-decoupled from components 1-3" one, previously the
23.7x-off pre-e48b5cc0 outlier) needing the largest, most different
fitted weight (0.699) of the three tested -- consistent with the
already-flagged, not-yet-root-caused observation that whatever is wrong
concentrates disproportionately in that one symmetry sector.

**Housekeeping**: both temporary edits (`Scrnc_driver.F`'s `Nmult=2`
line, `Scrnc_jj_rcc_assemble.F`'s diagnostic print) reverted via `git
checkout --`, confirmed empty `git diff` on both files. `xscrnc`
rebuilt and md5-reconfirmed to the exact pre-session value
(`0e188fe2576a8a8eea874a35e7078f5a`), matching fork #5's own
end-of-session checksum -- no drift. New test dirs (not committed, data
only), all under `/blue/ufhpc/perera/scrnc_test/`:
`h2_scrnc_rccd_{rhf,uhf}_sto3g_nmultdiag_0823`, `h2o_scrnc_rccd_
{rhf,uhf}_sto3g_nmultdiag_0823`, `h2_scrnc_rccd_rhf_sto3g_nmultdiag2_0823`,
`h2o_scrnc_rccd_rhf_sto3g_nmultdiag2_0823` (the broadened, all-
diagonal-elements rerun). Also reconfirmed, in passing, the pre-existing
"UHF's own log never prints the RPA JJ matrix" quirk (flagged
2026-08-20, still not chased -- irrelevant to today's conclusion since
all UHF numbers used here came from long-established prior values plus
today's directly-measured bare-v, not from parsing a UHF run's own JJ
print).

---

## SESSION 2026-08-22 (autonomous fork #5, R-side/Iside=1 list37-34
identity chase -- the concrete next step flagged by fork #4 below --
TRACED, CROSS-CHECKED, DIAGNOSTIC-TESTED, NOT A BUG: R-side's list37/34
pair is confirmed self-consistent with its own Hbar (no corruption),
confirmed structurally DIFFERENT from list61/63 (T_S/T_T) by a real
0.94%/2.5%, but an empirical substitution test shows that gap barely
moves the response ratio (1.878x->1.883x, wrong direction) -- so the
already-shipped `(T37-T34)/2` formula is left AS-IS, no fix made,
1.88x/12.5% residual UNCHANGED, RESUME HERE) -- picked up fork #4's
own flagged loose end ("was the R-side formula ever independently
derived, or just symmetric-by-construction") and settled it by tracing
+ checksums + a throwaway diagnostic rebuild, per the project's
"test, don't theorize" doctrine.

**Traced Iside=1 vs Iside=2's physical meaning first** (`Scrnc_form_mubar_s.F`
in full, plus its two callers). `Iside=1` is the R(ket)-side, driven
by `Scrnc_lrspn_rcc_main.F` (`Do ISide = 1, 1` -- this file only ever
calls Iside=1); `Iside=2` is the L(bra)-side, driven by the analogous
`Scrnc_jj_rcc_lrspn.F`. Both build the similarity-transformed property
operator `Mbar = exp(-T) M exp(T)` (the routine's own header) that
becomes the RHS ("Y") of the `Hbar*C=Y` linear-response solve --
Iside=2 additionally folds in G-intermediates (lists 191/192, one-particle
density-matrix-type pieces) that Iside=1 does not, consistent with the
bra side needing extra Lambda-CC-style terms the ket side doesn't.
`LISTT1=90+100*(ISIDE-1)` (T1 list 90 for Iside=1, "L1" list 190 for
Iside=2) and `LISTT2AB/LISTT2B` follow the identical `100*(ISIDE-1)+...`
offset pattern, giving lists 37/34 for Iside=1 and 137/134 for Iside=2 --
same offset convention on both sides, not an accident of the fix.

**Confirmed, via a cross-checksum, that R-side's list37/34 are
self-consistent with the R-side's OWN Hbar (no staleness/corruption
bug -- this codebase's most common failure mode in this exact
pipeline).** `rcc_dwmbej.F`'s RHF-only ABAB block (lines 217-289) is
where the response Hbar itself gets built: `IMULT=1` reads `LISTTA=37`
into singlet Hbar (list 56), `IMULT=2` reads `LISTTB=34` into triplet
Hbar (list 54) -- the EXACT SAME list numbers `Scrnc_form_mubar_s.F`
reads for Iside=1's dressing. That routine already carries an existing
`_DEBUG_LVL0` checksum ("T2-in") at exactly this read, so no new
instrumentation was needed -- just a fresh run with the existing prints
on. Ran H2/STO-3G RHF RCCD fresh
(`/blue/ufhpc/perera/scrnc_test/h2_scrnc_rccd_rhf_sto3g_rside_trace_0822`,
current HEAD 67a4ed6b, `xringcc`/`xscrnc` confirmed up-to-date via
`make install`). Final-iteration checksums:
```
RCC_DWMBEJ RHF IMULT=1 (list 37): T2-in  -0.967238661727E-01   0.935550628739E-02
RCC_DWMBEJ RHF IMULT=2 (list 34): T2-in   0.159334637447E+00   0.253875266903E-01
Scrnc_form_mubar_s.F LEADB_RAW37 Iside=1 LISTT2AB=37 val=-0.0967238662
Scrnc_form_mubar_s.F LEADB_RAW34 Iside=1 LISTT2B =34 val= 0.1593346374
```
Exact match to every printed digit. **Settles the "was list37/34
maybe a different intermediate by the time mubar_s reads it" question
definitively: no, it's the identical, fully-converged (12-digit-stable
over ~25 DIIS-type iterations, watched the full convergence trace) T2
amplitude pair that builds the R-side's own singlet/triplet Hbar.**
Internally coherent -- not a plumbing bug.

**But confirmed list37/34 (R-side) are genuinely, structurally NOT
the same numbers as list61/63 (T_S/T_T, the ground-state-energy
formula's own channel decomposition from `RCC_RHF_CMPENEG`,
independently cross-validated against UHF's raw ABAB amplitude to
~4e-10 in fork #3's Finding 1)** -- list37=-0.0967238661727 vs
T_S=-0.0976373126652 (0.94% off); list34=+0.159334637447 vs
T_T=+0.163485120840 (2.5% off). Both pairs are separately converged to
12+ digit stability, so this is a real, reproducible structural fact,
not iteration noise or a stale read -- it means ACES's own RCCD
machinery independently converges TWO not-quite-identical
representations of what should be "the same" physical singlet/triplet
channel decomposition: one inside `rcc_cmpeng.F`'s energy path (lists
61/63), a separate one inside `rcc_dwmbej.F`'s response-Hbar path
(lists 37/34, via `RCC_RNABIJ`'s `SST002/SST003/SSTRNG` resorts of
primary reference lists 46/44 -- already traced in "continued 3"
below). This is very plausibly the SAME phenomenon fork #3's Finding 2
already flagged approximately (a from-scratch textbook Casida-form
Riccati solve reproduced ACES's own T_S/T_T to only ~0.9-2.5%, "not
tight enough," blamed vaguely on "DPD-packed contraction convention")
-- now pinned down to an exact, named pair of ACES's own internal
lists, not merely a textbook-derivation artifact.

**Consequence for the R-side formula**: applying the same
`(T37-T34)/2` identity used for the L-side fix to the R-side's OWN
pair gives `(-0.0967238661727-0.159334637447)/2 = -0.128029251810` --
which itself is ~1.94% off from UHF's independently-converged TRUE
spin-orbital ABAB amplitude (-0.130561216326), unlike the L-side's
essentially-exact (~4e-10) reconstruction. So the formula's underlying
identity is sound, but on the R-side it's being fed an input pair
that's itself ~1-2.5% off from the most-accurate available
representation.

**Empirical diagnostic test (per the project's own doctrine -- test,
don't theorize): does that ~1-2% input imprecision actually matter for
the final response ratio?** Temporarily edited `Scrnc_form_mubar_s.F`
to overwrite list37/34 with list61/63 for `RCCD .AND. ISIDE.EQ.1` only
(clearly marked "TEMP DIAGNOSTIC, NOT FOR SHIPPING" in the source),
rebuilt `xscrnc` (`a35b4ffc4e20791c071d792a8675d109`), reran H2/STO-3G
RHF RCCD fresh
(`/blue/ufhpc/perera/scrnc_test/h2_scrnc_rccd_rhf_sto3g_rsidediag_0822`).
**Result: the ratio got very slightly WORSE, not better** -- response-only
diagonal 0.032051820282 (list37/34, current shipped) ->
0.031970932279 (list61/63 substituted); UHF/RHF ratio 1.8784x ->
1.8831x, a 0.25% change in the RHF value, in the wrong direction. This
is decisive: **list37/34's ~1-2.5% amplitude-level gap from list61/63
is NOT a meaningful driver of the 1.88x residual** -- ruling out the
one concrete, testable hypothesis this session could point to, despite
initially looking like a plausible, comparably-sized culprit.

**Reverted the diagnostic edit cleanly**: `git checkout --
aces2/scrnc/scrnc_form_mubar_s.F`, confirmed empty `git diff`, rebuilt
`xscrnc` and confirmed its md5 (`0e188fe2576a8a8eea874a35e7078f5a`)
exactly matches the pre-diagnostic, HEAD-committed (67a4ed6b) binary.
Reconfirmed with a plain rerun of the same fresh test dir (the
`_rside_trace_0822` run above, done BEFORE the diagnostic edit) that
the committed binary reproduces the exact known baseline
(`0.032051820282`, 1.88x) -- no drift.

**Verdict: no code fix this session, and this is the right outcome,
not a cop-out.** Fork #4's loose end ("was the R-side formula ever
independently derived, or just symmetric-by-construction") is now
RESOLVED, not left open: the R-side pair is confirmed self-consistent
with its own Hbar (no bug), the shared `(T37-T34)/2` formula is
confirmed to be the theoretically-correct spin-adaptation identity for
whichever self-consistent singlet/triplet pair a side actually uses
(not specifically tied to matching list61/63's numeric values), and an
empirical test rules out "R-side amplitude imprecision" as a
meaningful lever on the residual gap. **The 1.88x (H2/STO-3G RCCD) and
up-to-12.5% (H2O/STO-3G RCCD) residuals from 67a4ed6b are UNCHANGED
and remain genuinely open.** The list37/34-vs-list61/63 discrepancy
itself is real and now precisely characterized, but appears to be a
separate, upstream (`rcc_rngdrv.F`/`RCC_RNABIJ`-level) ring-CC
amplitude-solver inconsistency, not a `scrnc`-level bug and not (per
today's test) the source of the response-ratio gap -- flagged for
whoever eventually wants to understand it, but explicitly NOT the next
lead to chase for closing the residual.

**Concrete next step, sharper than before**: mechanism (b) (Lead B,
the whole list37/39/34-substitution line of investigation across
fork #3, #4, "continued 3", and this session) now looks exhausted --
every testable variant of "which list should the RHF dressing formula
combine" has been tried, and none narrows the gap further than the
already-shipped fix. Fork #3's mechanism (a) -- `Scrnc_driver.F:251`'s
`Nmult=1` hardcode, which uses ONLY the singlet-channel Hbar for the
response and structurally discards the triplet-channel information
that both amplitude representations (list61/63 AND list37/34) agree is
non-negligible (~63% of T_T's weight in the ground-state energy
formula) -- was explicitly closed as "not to be re-litigated" earlier
in this thread on Ajith's own physics argument (spin-scalar
perturbation, Sz=0 triplet manifold, zero overlap), but given mechanism
(b) is now exhausted without closing the gap, that argument (not just
the code) may be worth Ajith revisiting directly: does a spin-scalar
perturbation's response TRULY have zero triplet-channel contribution
once you're computing a MATRIX response property (JJ matrix, built
from Hbar eigenstructure) rather than a scalar (ground-state energy),
or could the matrix-vs-scalar distinction from fork #2's own
"Interpretation" section reopen this specific point? Needs Ajith's
judgment call, not another code-side test -- this session found no new
code-level lead to chase.

**Housekeeping**: `git diff`-clean at end of session (only the
pre-existing, unrelated dirty files noted at session start:
`aces2/vcc/.bwcc.F.swp`/`.drrng.f.swp`/`.rngdrv.f.swp` deleted,
`aces3/config.log`/`config.status` modified, plus untracked
`aces2/scrnc/scrnc_sigma_driver` and two `aces3` make logs -- none of
this session's doing, not touched). `aces2/bin/xscrnc` confirmed
md5-matched back to the exact 67a4ed6b committed state
(`0e188fe2576a8a8eea874a35e7078f5a`). New test dirs (not committed,
data only), both under `/blue/ufhpc/perera/scrnc_test/`:
`h2_scrnc_rccd_rhf_sto3g_rside_trace_0822` (checksum cross-check,
clean binary), `h2_scrnc_rccd_rhf_sto3g_rsidediag_0822` (throwaway
list61/63-substitution diagnostic, code since reverted).

---

## SESSION 2026-08-22 (autonomous fork #4, empirical Lead A/Lead B test
per fork #3's own handoff -- REAL BUG FOUND AND FIXED: RHF's
`Scrnc_form_mubar_s.F` mubar T2-dressing formula was reconstructing the
WRONG spin-orbital quantity; fixing it cuts the H2/STO-3G RCCD
response-only RHF/UHF ratio from 5.18x to 1.88x and improves EVERY ONE
of H2O/STO-3G's 49 diagonal JJ elements, no regressions -- committed,
RESUME HERE) -- picked up fork #3's two concrete leads and tested both
empirically (checksums + real pipeline reruns via sbatch), per the
project's own "test, don't theorize" doctrine, rather than reasoning
about the code further.

**Lead A (`scrnc_driver.F:251`'s `Nmult=1` hardcode) -- NOT a fresh
lead, already tested and closed by Ajith himself; did not re-test
blindly.** Reading the full prior history (this log, 2026-08-19
"continued 6" and "continued 7") shows `Nmult=2` (singlet+triplet
Hbar-diagonalization, 0.25/0.75-recombined) was the ORIGINAL
architecture, was tried, and was explicitly reverted after Ajith's
direct call: *"NO, they should not decouple. We do not do that in CCSD
route"* -- backed by a first-principles Wigner-Eckart argument (a
spin-scalar density-type perturbation has zero overlap with the Sz=0
triplet excitation manifold, so the RESPONSE must come from the singlet
channel alone). A LATER session (2026-08-19 "continued 7") independently
re-derived this exact same hypothesis from scratch, wrote a `Wt=1.0/0.0`
test edit, and got a null result -- because `Nmult` was already
hardcoded to 1, making the weighted-recombination branch dead code.
Both routes converge on the same closed verdict. Also checked the
`RCL_HEFFXC.F` `LISTDUMP`-hardcoded-to-56 bug flagged in the task: it
lives in `ringcl` (a SEPARATE diagnostic executable, `Rcl_check_rcchbar`
+ `Rcl_heffxc`), not in `scrnc_driver.F`'s own production `Imult` loop --
scrnc's own equivalent write-back (`scrnc_form_rcc_hbar.F`) was already
fixed to be `Imult`-aware in this same thread's earlier work, so the two
bugs don't interact. **Verdict: Lead A is a settled, closed physics
question (Nmult=1/singlet-only-diagonalization is correct), not
re-litigated -- consistent with "Stay on Agreed Debug Path."** No code
changes, no test run for this lead.

**Lead B (`scrnc_form_mubar_s.F:286-306`'s RHF-only `2*list37-list34`
mubar dressing) -- tested empirically via `_DEBUG_LVL0` checksum print,
CONFIRMED and FIXED.** Added two guarded prints (`LEADB_RAW37`/
`LEADB_RAW34`) right before the existing SSCAL/SAXPY combination, then
rebuilt `xscrnc` (`make install` in `aces2/scrnc/`, confirmed via
`md5sum` the installed binary actually changed and contains the new
debug strings -- `aces2/bin` is NOT auto-updated by a plain `make`, only
`make install`, a new trap worth remembering alongside the existing
`aces2/bin` PATH one) and reran H2/STO-3G RHF RCCD fresh
(`/blue/ufhpc/perera/scrnc_test/h2_scrnc_rccd_rhf_sto3g_leadb_0822`).
Result, decisive:
```
Iside=2 (L-side, lists 137/134): list37=-0.0976373127, list34=+0.1634851201
  -- EXACTLY matches RCC_RHF_CMPENEG's converged T_S=-0.0976373126652 (list 63)
     and T_T=+0.163485120840 (list 61), to ~8 printed digits.
Iside=1 (R-side, lists 37/34):   list37=-0.0967238662, list34=+0.1593346374
  -- a DIFFERENT, related-but-distinct pair (matches the raw numbers already
     printed in the 2026-08-19 "continued 6" entry -- not T_S/T_T themselves).
```
This settles the "list34~T_T" premise flagged UNCONFIRMED in fork #3's
entry: **exactly true for the L-side (Iside=2) amplitudes actually used
in the property construction**, not just approximately. Combined with
Finding 1's already-proven exact identity `T(ABAB) = (T_S - T_T)/2`
(the SAME quantity UHF's own raw, undressed list37 amplitude directly
IS), this makes a precise, derivable prediction: the physically correct
RHF reconstruction of that spin-orbital quantity is **`(list37 -
list34)/2`**, not the CCSD-style `2*list37 - list34` currently shipped
(a formula empirically tuned to improve H2O in `498542f8`, 2026-08-19,
never actually derived -- and already flagged as suspect once list34's
true identity was traced in "continued 3"). Numerically:
`(T_S-T_T)/2 = -0.130561216753`, matching UHF's raw `T_ABAB =
-0.130561216326` to 9 significant figures (vs. the shipped formula's
`2*T_S-T_T = -0.358759746`, which matches nothing UHF actually has).

**Implemented and verified the fix** in `Scrnc_form_mubar_s.F`: guarded
`IF (RCCD)` (mirroring the existing list-selection guard) to scale by
`HALF`/`SAXPY` by `HALFM` instead of `TWO`/`ONEM`; CCSD's own
`2*T37-T39` formula (for its genuinely antisymmetrized amplitudes) is
completely untouched, still using the original `TWO`/`ONEM` path.
Rebuilt, reran three systems fresh via sbatch:
- **H2/STO-3G RHF RCCD** (`h2_scrnc_rccd_rhf_sto3g_leadbfix_0822`):
  response-only diagonal `0.01161070026` (old) -> `0.03205182028`
  (new), vs UHF's `0.06020595989` -- **ratio 5.18x -> 1.88x**. Final
  plusbare diagonal gap -5.82% -> -3.37%.
- **H2O/STO-3G RHF RCCD** (`h2o_scrnc_rccd_rhf_sto3g_leadbfix_0822` vs
  fresh baseline `h2o_scrnc_rccd_rhf_sto3g_leadbbase_0822`, both against
  `h2o_scrnc_rccd_uhf_sto3g_cleanbaseline`): **every single one of the
  49 diagonal JJ elements improved**, none regressed -- e.g. JJ(1,1)
  -6.16%->-4.74%, JJ(9,9) -7.95%->-5.47%, JJ(17,17) -12.79%->-8.78%,
  JJ(41,41) -13.54%->-9.39%; smallest improvements on the
  already-nearly-converged elements (JJ(7,7) -1.25%->-0.96%), largest on
  the worst-offending ones -- a systematic, monotonic improvement
  pattern, not noise.
- **H2/STO-3G RHF DRCCD** (`h2_scrnc_drccd_rhf_sto3g_leadbfix_0822` vs
  `h2_scrnc_drccd_rhf_sto3g_pilot0822`): **bit-identical**, confirming
  the `IF (RCCD)` guard correctly excludes DRCCD (which keeps list 39,
  never reaches the list-34 branch) -- no regression, as expected.

**Verdict: real, substantial, empirically-verified partial fix,
committed. Residual gap (1.88x on H2, up to ~12.5% on H2O) NOT fully
closed -- genuinely open.** The clean match between this session's
from-scratch-derived prediction (1.88x) and fork #3's independent
"crude proxy" estimate (1.84x, using the WRONG `2*T_S-T_T` combination
by coincidence landing close) is a striking cross-validation that this
mechanism is real and correctly identified, not a fluke. What's left
unexplained: (1) the Iside=1 (R-side) amplitude pair (list37/34, NOT
equal to T_S/T_T) was left using the SAME `(T37-T34)/2` formula as
Iside=2 by symmetry of the guard, but its numeric identity was never
independently confirmed to be the analogous "R-side" spin-orbital
reconstruction -- if it's actually a different physical object (e.g. an
unrelaxed intermediate), applying the same halved-combination formula
there may be only incidentally reasonable, not derived; (2) H2's
residual 1.88x and H2O's up-to-12.5% residuals are NOT accounted for by
this mechanism alone -- something else is still contributing,
consistent with fork #3's own point that Lead B was "a real,
contributing effect, not the whole story."

**Housekeeping**: fix is `git diff`-clean (only
`aces2/scrnc/scrnc_form_mubar_s.F`/`.f`, the debug prints kept
`_DEBUG_LVL0`-gated per convention, matching this thread's existing
instrumentation style). `aces2/bin/xscrnc` and `aces2/scrnc/xscrnc`
rebuilt and md5-matched to the new source (final binary
`0e188fe2576a8a8eea874a35e7078f5a`). New test dirs (not committed, data
only), all under `/blue/ufhpc/perera/scrnc_test/`:
`h2_scrnc_rccd_rhf_sto3g_leadb_0822` (checksum-only probe run),
`h2_scrnc_rccd_rhf_sto3g_leadbfix_0822` + `_finalcheck` (fixed, H2),
`h2o_scrnc_rccd_rhf_sto3g_leadbfix_0822` + `_leadbbase_0822` (fixed +
matched-fresh-baseline, H2O), `h2_scrnc_drccd_rhf_sto3g_leadbfix_0822`
(DRCCD regression check). **Next step, concretely**: instrument the
R-side (Iside=1) list37/34 pair the same way L-side was just settled
here -- identify what physical quantity they actually are (not yet
determined) -- before trusting the fix's Iside=1 branch as fully
derived rather than just symmetric-by-construction; then re-attack the
remaining ~1.88x/~12.5% residual with that clarified.

---

## SESSION 2026-08-22 (autonomous fork #3, from-scratch analytic/symbolic
H2/STO-3G RCCD+DRCCD derivation, the flagged next step off fork #2 below
-- gets a genuine THIRD independent number, RESUME HERE) -- Ajith's own
"Pilot-Code Debugging" doctrine applied for real this time: derived the
single-amplitude RCCD/DRCCD Riccati equations and the A+BT Hbar_SS block
by hand/in a standalone script, using REAL ACES MO integrals (not
invented numbers) pulled from the existing H2/STO-3G RHF/UHF RCCD/DRCCD
run logs. **Result: one exact, unambiguous independent confirmation at
the amplitude level; the response-MAGNITUDE ratio (5.18x/1.15x) is NOT
cleanly reproduced from scratch, but this session pins down two
concrete, line-numbered code mechanisms (not a vague "matrix vs scalar"
argument) that plausibly drive it, one of which was already flagged and
half-investigated in this project's own history (2026-08-20 "continued
3" below) and never followed up. Genuinely still open**, per the task's
explicit allowance for an honest "not fully closed" verdict -- but with
a much narrower, more actionable "what's blocking closure" than before.

**Setup / real numbers used** (pulled from `run_scrnc.log` in
`h2_scrnc_rccd_{rhf,uhf}_sto3g_{verify3,reverify}` and
`h2_scrnc_drccd_{rhf,uhf}_sto3g_pilot0822`, H2/STO-3G, R=0.742 Angstrom,
sym=off): RHF SCF gives 1 occ MO (eps_i=-0.5777270364 Ha) + 1 virt MO
(eps_a=+0.6690811064 Ha), and the printed 4-index MO integrals give
(ii|aa)=0.66336391034632 (J), (ia|ia)=0.18132242419990 (K) -- the only
two two-electron integrals that matter for a 1-occ/1-virt system.
**Confirmed UHF's alpha and beta orbitals are BIT-IDENTICAL to RHF's**
(same eigenvalues to 10 digits, `<S**2>=0.0000000` exactly) -- i.e. UHF
does NOT break spin-symmetry for H2/STO-3G at this bond length (well
inside the Coulson-Fischer stability radius). This is an important,
previously-unstated clarification: the RHF/UHF response gap under study
is NOT caused by different mean-field solutions -- it is purely about
how the two pipelines construct/diagonalize the response Hbar from
IDENTICAL orbitals/integrals.

**Finding 1 (exact, clean success) -- the spin-adaptation amplitude
identity holds to ~4e-10**: RHF's RCCD run internally solves for TWO
spin-adapted "channel" amplitudes (singlet, list 63 in
`RCC_RHF_CMPENEG`, converged T_S=-0.0976373126652; triplet, list 61,
converged T_T=+0.163485120840), used together (Wt=0.25/0.75) ONLY for
the ground-state energy. The standard closed-shell CI/CC spin-orbital
reconstruction identity `T(ABAB) = (T_S - T_T)/2` predicts
-0.130561216753. **ACES's own UHF run, converging the spin-orbital ABAB
amplitude completely independently (different code path, different
process), gives -0.130561216326 -- and RHF's own internal
`Rcl_check_uhf` diagnostic reconstruction gives -0.130561216375.** All
three agree to ~4e-10, a clean, non-ACES-derived, textbook-identity
confirmation that RHF's two channel amplitudes and UHF's spin-orbital
amplitude are exactly consistent representations of the same physical
ground state. Nothing is wrong with the underlying amplitude solve,
RHF vs UHF -- confirms (again, now via an exact match rather than a
plausibility argument) that the already-known ground-state identity
(`rcc_cmpeng.F`) is the correct, complete picture at the amplitude
level.

**Finding 2 (partial, ~1-3% not exact) -- textbook closed-shell RPA/TDHF
(Casida-form) A,B built from these same real integrals**: solving the
scalar Riccati equation `B*t^2 + 2*A*t + B = 0` with the standard
closed-shell singlet (A=Delta+2K-J, B=K), triplet (A=Delta-J, B=-K), and
direct-RPA/dRPA (A=Delta+2K, B=K) forms reproduces ACES's own converged
T_S, T_T, T_dRPA and their ringcl-diagonalized Hbar eigenvalues (singlet
25.50578519 eV, triplet 15.48325067 eV, dRPA 43.23234947 eV) to the
right sign, right order of magnitude, and correct qualitative pattern,
but only to ~0.9-2.5% precision -- NOT machine precision. Good enough to
validate the basic single-amplitude Riccati/A+BT STRUCTURE (Eq 26/27 of
the Rishi/Perera/Bartlett paper, already established as correct in this
log), not tight enough to trust a downstream high-precision numeric
prediction of the exact 5.18x/1.15x response ratio from these
coefficients alone. The residual few-percent mismatch most likely
reflects ACES's actual DPD-packed RCCD contraction (`rcc_newt2.F` et
al.) using some normalization/intermediate convention not identical to
the plain textbook formula -- tracing that exactly, line by line, was
out of scope for this session.

**Finding 3 (the actual response-construction mechanism, read directly
from the Fortran, not guessed) -- TWO concrete, quantifiable asymmetries
between how RHF's and UHF's JJ-matrix pipelines are built**:

(a) `scrnc_driver.F:251` hardcodes `Nmult=1` -- RHF's production
JJ-matrix uses ONLY the singlet-channel Hbar/amplitude (T_S), never
combining in the triplet channel (T_T) the way the ground-state energy
formula does (0.25/0.75 weighted). Given Finding 1's now-EXACT identity
`T(ABAB)=(T_S-T_T)/2`, the TRUE spin-orbital ground-state amplitude that
should dress any T2-based property/transition-moment construction is
NOT recoverable from T_S alone -- it structurally requires T_T too.
This upgrades the prior sessions' vague "scalar vs matrix, no a priori
reason the reconciliation should extend" argument (see the 2026-08-22
fork #2 entry below) into a precise statement of WHAT specific
information is being discarded and why that's expected to matter.

(b) `scrnc_form_mubar_s.F` lines 286-306: the transition-moment (mubar)
T2-dressing contraction is genuinely DIFFERENT for RHF vs UHF, not just
differently-normalized. For RHF (`IUHF.EQ.0`), it builds
`SCR = 2*T2AB(list 37) - T2(list 34)` before contracting into MAIBAR(ai)
(RCCD only; a leftover from the already-shipped 2026-08-19/`498542f8`
fix). For UHF, it uses the raw `T2AB` (list 37) alone, no combination,
no factor of 2. **Critically, this project's OWN prior session
(2026-08-20 "continued 3" below) already traced list 34's construction
via `RCC_RNABIJ`/`RCC_SST003` and found it is NOT a permuted/exchange
view of the singlet-channel ABAB tensor family at all -- it is built
from a completely separate reference list (44, `SPCAS2='AAAA'`) and is
literally "ring-CC's TRIPLET channel's own primary amplitude."** That
session found substituting list34 into the Hbar CONSTRUCTION itself
(not mubar) made a different symptom (2 negative eigenvalues) WORSE,
and explicitly flagged that the ALREADY-SHIPPED `Scrnc_form_mubar_s.F`
fix rests on the same now-suspect premise -- a loose end raised then and
never revisited (log lines ~1373-1393 as of this writing). This session
reconnects that dangling thread to the present investigation: the
"2*T37-T34" RHF-only mubar dressing is, in effect, mixing the triplet
channel's own amplitude into the (supposedly singlet-only, Nmult=1)
response construction, via a different, unaudited code path than the
Hbar block that already rejected the analogous substitution.

**Quantitative probe of mechanism (b)** (script `Section 4`, flagged
UNCONFIRMED -- did not independently verify list 34's actual numeric
value equals T_T bit-for-bit in the same run, would need an instrumented
rerun or a direct checksum comparison): taking list34 ~ T_T at face
value, RHF's effective dressing factor is `2*T_S - T_T = -0.358760`
vs UHF's raw `T_ABAB = -0.130561`. Using `[(1+T_ABAB)/(1+2T_S-T_T)]^2`
as a crude proxy for the mubar-squared response-magnitude ratio gives
**1.84x** -- right DIRECTION (UHF bigger, matching the real 5.18x's
sign) and the right order of magnitude, but undershoots substantially
(~35% of the observed ratio). Consistent with mechanism (b) being A
real, contributing effect, not the WHOLE story -- mechanism (a) (the
Nmult=1 restriction itself, which this crude proxy does not separately
model) plausibly contributes independently on top of it.

**Verdict: genuinely still open, narrower than before.** This session
does NOT close the RHF/UHF response-gap question either direction --
did not reproduce 5.18x/1.15x to convincing precision from first
principles, so this is not a clean "confirmed real physics" result; but
it also identifies two SPECIFIC, real, code-located mechanisms (not
speculative) that plausibly drive at least part of the gap, one of which
(list 34's true identity, and the shipped mubar fix's now-questionable
theoretical grounding) is a previously-flagged, never-resolved loose end
in this project's own history. **Concrete next steps, sharper than "do
the derivation" (now done, partially)**: (1) instrument
`Scrnc_form_mubar_s.F` (or add a `_DEBUG_LVL0` checksum) to print list
34's actual numeric value for H2/STO-3G RHF RCCD and directly compare
against T_T=+0.163485120840 -- settles whether the "list34~T_T" premise
used in the quantitative probe above is exactly right or only
qualitatively suggestive; (2) re-derive what the CORRECT closed-shell
RHF mubar dressing formula should be (a real, from-the-paper-or-
literature spin-adapted transition-moment construction, analogous to
Eq 46/49/65 for the energy) rather than reusing/rejecting CCSD's
antisymmetrized-amplitude pattern by trial and error; (3) trace why the
textbook Casida-form Riccati equation (Finding 2) is off by 1-3% against
ACES's own numbers -- either accept this as expected (real convention
difference) or find the source, since it currently blocks a fully
precise from-scratch numeric prediction of the response ratio.

Script: `aces2/scrnc/independent_pilot/h2_sto3g_analytic_rccd.py`
(matches this directory's existing "Pilot-Code Debugging" convention,
e.g. `reciprocity_pilot.py`). No ACES code changes this session --
derivation/analysis only. Doc + script commit.

---

## SESSION 2026-08-22 (autonomous fork #2, minimal-case pilot test per
Ajith's own approved next step -- confirms the effect is real, gives a
concrete physical/mathematical reason ground-state energy matches but
response doesn't, RESUME HERE) -- Ajith, still skeptical that exact
RHF/UHF response identity is achievable, approved a cheap independent
pilot-case test (this project's own "Pilot-Code Debugging" doctrine)
before committing to a full Lagrangian re-derivation: use the smallest
possible tractable system to either (a) confirm the gap is generic
rCCD physics, or (b) discover it vanishes in the minimal case, reopening
a tractable bug hunt.

**Note on process**: the fork initially assigned this task did ZERO
actual work on its first turn (0 tool calls, just echoed a status
line back after apparently getting confused inheriting the
coordinator's own "fork is running" message) -- had to be explicitly
resumed/corrected before real investigation started. Flagging this as
a recurring failure mode for future forked sub-investigations on this
thread.

**Setup**: H2/STO-3G is about as minimal as this pipeline supports --
Nao_pair=4 (2 AOs), 1 occ + 1 virt spatial MO for RHF; for UHF, 1
occ-alpha/1 occ-beta/1 virt-alpha/1 virt-beta, with AA and BB doubles
channels trivially zero (confirmed in the log output: "The AA/BB
contribution to the correlation energy is: 0.0000000" -- Pauli
exclusion, only 1 electron per spin, so ALL of UHF's rCCD correlation
comes from the single ABAB channel). Ran fresh RHF+UHF pairs for BOTH
RCCD (reused existing validated `h2_scrnc_rccd_{rhf,uhf}_sto3g_verify3`/
`_reverify` dirs, current HEAD 42663467) and DRCCD (new dirs built this
session, `h2_scrnc_drccd_{rhf,uhf}_sto3g_pilot0822`, since no post-fix
H2 DRCCD UHF run existed yet -- needed a fresh GENBAS/ZMAT.BAS copy
into a clean dir, first attempt crashed at xjoda with those files
missing, fixed by copying from an existing validated dir).

**Ground-state confirms EXACT RHF/UHF identity, cleanly, in the
minimal case**: E(SCF) bit-identical (-1.116651278229 a.u., both
methods, both RCCD/DRCCD). Converged correlation energy: DRCCD
bit-identical to all 12 printed digits (-1.137618408411 a.u., RHF AND
UHF); RCCD agrees to 10 sig figs (-1.128488116359 vs -1.128488116293,
diff ~6.6e-11, consistent with cc_conv=9 convergence-threshold noise,
not a real discrepancy). This is the cleanest confirmation yet of the
already-known `rcc_cmpeng.F` ground-state identity, now demonstrated
in the smallest possible system.

**Response property (JJ matrix) does NOT match, and the raw
correlation-only piece differs by a LARGE factor, not just a few
percent**: the `jjlrespn_resp_packed.dat` (response-only, before bare
added) diagonal element for RCCD is 0.01161070026110285 (RHF) vs
0.06020595989104796 (UHF) -- **UHF's raw correlation response is
5.18x LARGER than RHF's**, same system, same converged ground-state
amplitudes' underlying physics, same SCF orbitals. For the final
"plusbare" metric this dilutes down to the historically-quoted ~-5.8%
(RHF-vs-UHF on the 0.786/0.835 diagonal element), because the bare
(uncorrelated) piece is much bigger than the response piece and is
identical between methods -- but the underlying CORRELATION
contribution itself is off by a factor of >5, not a small percentage.
DRCCD shows the same direction, much smaller: response-only piece
0.026120509141200 (RHF) vs 0.030088499441843 (UHF), a 1.15x ratio
(-0.49% on the final plusbare diagonal) -- consistent with the
already-established ~5-7x-smaller-than-RCCD pattern from the H2O full-
diagonal check above, confirming that pattern generalizes down to the
minimal 2-orbital case too.

**A genuine independent check DID pass, and it's informative**: the
response-only matrix's (1,1) and (1,4) [the two AO-pair indices
corresponding to the two on-atom, "gerade"/"ungerade"-type density
combinations for homonuclear H2] are numerically exact negatives of
each other in BOTH RHF and UHF, BOTH RCCD and DRCCD (e.g. RHF-RCCD:
+0.01161070026110285 / -0.01161070026110285). This is a real physical
sum rule -- the response matrix's row for a given AO-pair index must
sum to zero under this specific combination (a particle-number-type
conservation constraint on the density-density response), and it is
satisfied exactly, independently, by both RHF's and UHF's converged
solution. This does NOT resolve the RHF/UHF question (it holds
separately and identically for each method, it isn't a cross-method
identity), but it's a real, passed consistency check confirming
neither method's response calculation is producing nonsense -- both
are internally physically consistent, they just disagree with each
other in magnitude.

**Interpretation / likely mechanism (a genuine reason, not just
another empirical pattern-match)**: the ground-state energy is a
SCALAR contraction (T2 dotted into the integrals) that Ajith's own
paper (`rcc_cmpeng.F`, Eq 46/49/65) proved can be exactly reconciled
between the spin-adapted (RHF, singlet/triplet channels) and
spin-orbital (UHF, AAAA/BBBB/ABAB/BABA channels) parameterizations of
rCCD's non-antisymmetric amplitudes. The JJ-matrix / response property
is NOT a scalar -- it's a full matrix built from diagonalizing the
Hbar operator, i.e. it depends on the amplitudes' and Hbar's detailed
MATRIX/operator structure, not just one contracted scalar. There is no
a priori reason a reconciliation that works for one specific scalar
contraction of a non-antisymmetric-amplitude object should also hold
for every matrix element built from that same object's eigenstructure
-- this is a plausible, mechanistically-grounded reason the ground
state can match exactly while the response doesn't, rather than
requiring either "there's still a missing correction term" or "no
identity should be expected anywhere at all."

**Verdict: strengthens Ajith's skepticism with a mechanism, but does
NOT fully close the case.** The gap is real, large in the fundamental
correlation-only quantity (not a small numerical artifact), present
even in the most minimal possible system, and now has a plausible
mathematical reason to exist (scalar-vs-matrix reconciliation) rather
than just an empirical "smooth spread" observation. It does not
prove no correction is possible -- only that if the ground-state
scalar trick doesn't automatically extend to the response matrix,
finding out whether ANY exact spin-adapted response formula exists
requires the derivation itself, not more numerical fishing. **No
coefficient/correction implemented or attempted** (nothing here
suggests a specific term to test, consistent with the project's
standing rule against unvalidated guessing).

**What would move this forward next, if picked up again**: this
minimal H2/1-occ/1-virt case is now FAR more tractable than H2O ever
was for an actual from-scratch analytic or symbolic derivation of the
ring-CCD/RPA response equations (per this project's own "Pilot-Code
Debugging" doctrine) -- with only a single ABAB doubles amplitude in
play for UHF and one spin-adapted singlet amplitude for RHF, it may be
possible to write out the full A/B RPA-type equations by hand or in a
short standalone script and get a THIRD, fully independent number to
compare against both ACES numbers, which would definitively show which
(if either) is "correct" rather than just observing they disagree.
NOT attempted this session (out of scope/time for a single pass) --
flagged as the concrete next step, not started.

Test dirs: `h2_scrnc_rccd_{rhf,uhf}_sto3g_{verify3,reverify}` (existing,
reused), `h2_scrnc_drccd_{rhf,uhf}_sto3g_pilot0822` (new this session).
Committed doc-only.

---

## SESSION 2026-08-22 (autonomous fork #1, chasing the {JJ11,JJ33}-vs-
{JJ22,JJ44} ~2x-ratio lead from the 2026-08-21/5dbee126 entry below,
RESULT: THAT LEAD DOES NOT HOLD UP, superseded by entry above) -- Ajith expressed
renewed skepticism that exact RHF/UHF response identity is achievable
here at all. Went back to the concrete "clean 2x split" clue flagged as
the next step and tested it rigorously rather than assuming it pointed
to a findable missing term.

**First: sym=off, so there is no code-level symmetry-irrep story to
tell.** Both `_wtcheck` test dirs' ZMATs use `sym=off` (C1) --
confirmed via the `_DEBUG_LVL0` "Isympert(1..Nao_pair)" diagnostic
dump: every one of the 49 AO-pair entries is Irrepx=1. `Iao_pair`
in `Scrnc_jj_rcc_assemble.F`/`Scrnc_lrspn_rcc_main.F` is a flat AO-pair
list position (Nao_pair=49=7x7 for H2O/STO-3G, all (mu,nu) AO index
pairs, not a triangular/unique-pair count), not a molecular-symmetry
irrep index. So "JJ11/22/33/44" are just diagonal entries #1-4 of a
49-entry list -- there was never an Isympert-level irrep to bipartition
by in the first place.

**Extracted all 49 diagonal RHF-vs-UHF %-differences (not just the 4
historically tracked ones), from the `jjlrespn_plusbare_packed.dat`
dumps in the `h2o_scrnc_rccd_{rhf,uhf}_sto3g_wtcheck*` dirs.** Result:
the %-differences do NOT split into two discrete populations at all --
they form a smooth, continuous spread from -1.14% to -14.42% across 22
distinct AO-pair-equivalence classes (repeated values across e.g.
indices {2,8,14,44} or {4,22,28,46} simply reflect water's real,
still-numerically-present C2v symmetry equivalences, unexploited by
`sym=off` but not erased from the integrals). The apparent clean ~2x,
two-group split reported previously was an artifact of the historical
convention of only ever quoting the FIRST four diagonal elements
(JJ11-JJ44) -- a 4-point sample that happened to land near a bimodal
pattern by chance, not a real structural signature. Checked the
obvious confound (large relative error only because of a small
denominator on a small-magnitude element) too: correlation between
log|UHF magnitude| and |%diff| across all 49 points is ~0.10 --
essentially uncorrelated, ruling that out as well.

**Cross-checked against DRCCD (`h2o_scrnc_drccd_{rhf,uhf}_sto3g_postfix`,
whose `jjlrespn_plusbare_packed.dat` files exist even though UHF's own
`xscrnc.out`/`run_scrnc.log` prints no final "RPA JJ matrix" table --
same pre-Symmet2 packed dump is written regardless).** Same result:
a smooth continuum, -0.02% to -1.96%, same 49-AO-pair-equivalence-class
structure (identical index groupings to RCCD, confirming the
underlying AO-pair symmetry bookkeeping is shared and correct across
both methods) -- but critically, the *ranking* of which class has the
worst %-error is DIFFERENT between RCCD and DRCCD. The index-4
class ({4,22,28,46}) is one of the worst offenders in RCCD (-13.12%,
2nd-worst class) but sits mid-pack in DRCCD (-0.21%, actually BETTER
than the index-1 class's -0.34%). A single shared missing term tied to
one AO-pair/symmetry sector should produce the same relative ranking in
both methods, just scaled -- it doesn't. This is a second, independent
line of evidence against the "one findable omission in one sector"
reading of the original lead.

**Verdict: (c), leaning toward (b)/Ajith's skepticism, not (a).** The
specific concrete lead from 5dbee126 (clean {1,3}-vs-{2,4}, ~2x ratio)
does not survive being checked against the full 49-point diagonal or
against DRCCD -- it was a small-sample artifact, not a real bimodal
symmetry-sector split. What actually replaces it is a smoothly graded,
pervasive RHF/UHF difference across essentially the WHOLE AO-pair
response space (not concentrated in one or two components), roughly
5-7x larger in absolute % terms for RCCD than DRCCD but with a
similar-shaped (not identically-ranked) spread in both. A smooth,
pervasive spread is much more consistent with a graded consequence of
rCCD/drCCD's proven non-antisymmetric amplitude structure (paper
Appendix VII.A, `rcc_cmpeng.F` precedent) than with a single missing
Hbar term that would be expected to show up as a sharp discontinuity
isolated to a specific sector. This does not *prove* Ajith's null
hypothesis (a graded spread is also what a genuinely-derivable but
not-yet-derived spin-adapted response correction might look like, the
same way the ground-state energy fix was smooth once found) but it
does remove the one concrete piece of evidence that had been cited
FOR "this is still a findable bug" -- that evidence no longer holds up.

**No code changes.** Per this project's own repeated lesson (CORR_2P's
-1.5 guess, FACT2's borrowed DCCSD value, both later shown not to
work), did not attempt to prototype any correction from this data --
a smooth 49-point spread with no clean discontinuity gives no natural
coefficient or term to test. **What would actually move this forward
next, if picked up again**: the still-not-done real spin-adapted
RHF/UHF *response*-property derivation (the response analog of Eq
46/49/65 that fixed the ground-state energy) is the only lever left
that isn't just more empirical pattern-hunting -- everything findable
by inspecting existing numbers has now been checked (Hbar formula
literature-matched and correct per the paper, downstream JJ-assembly
chain audited clean, pipeline/Hbar-read validated at the eigenvalue
level, and now this full-diagonal/DRCCD cross-check). Doc-only commit.

---

## SESSION 2026-08-21 (autonomous fork, null-hypothesis evaluation per
Ajith's explicit redirect, RESUME HERE) -- Ajith pushed back on the
previous entry's framing: rather than assume a corrective spin-adapted
response formula exists and go hunting for it (the energy case as
template), first seriously evaluate the NULL HYPOTHESIS -- is exact
RHF=UHF response identity for rCCD/drCCD *structurally impossible*,
full stop, regardless of spin-adaptation? Verdict: **NOT supported.**
Found concrete, structured evidence the residual is still an
unexplained bug, not a physical floor -- do not close this thread as
"expected/unfixable."

**Evidence against the null hypothesis:**

1. The non-antisymmetric-amplitude issue (paper Appendix VII.1) breaks
   the *naive* CCD-style RHF/UHF argument, but for the closely-related
   ground-state energy, an EXACT spin-adapted fix was found and
   independently verified to 8 decimal places (`rcc_cmpeng.F`'s
   `E_flip`). That's a real, in-repo existence proof that non-
   antisymmetric amplitudes do NOT by themselves prevent exact RHF/UHF
   agreement when the spin-adaptation is done correctly -- it's not a
   free pass to assume the response case is unfixable too.

2. For DRCCD specifically, `scrnc_jj_rcc_assemble.F`'s own header
   comment (2026-08-19, already in the code, apparently not carried
   forward into recent sessions' framing) names a CONCRETE missing
   piece: UHF's direct-RPA Hbar has a second eigenvalue that RHF's
   single response channel never builds ("DRCCD's RHF response is
   missing real physics (an unbuilt second channel)"), currently
   papered over by an admittedly-wrong patch (`Spntsing`). That's an
   engineering gap with a named fix path, not a physical limit.

3. For RCCD, re-tested the obvious hypothesis (singlet/triplet channel
   decoupling, `Wt=0.25/0.75` in `scrnc_jj_rcc_assemble.F`, matching
   `Check_rhf_rcc.F`'s ground-state-energy formula that `rcc_cmpeng.F`
   later proved WRONG for energy) by flipping it to `Wt=0.25/0.25`
   (matching the verified `E_flip` shape), rebuilt, reran H2O/STO-3G
   RHF RCCD fresh. **Bit-identical to baseline** -- turned out this is
   dead code: `scrnc_driver.F` hardcodes `Nmult=1` (singlet-only, since
   `498542f8`, 2026-08-19) on its own solid first-principles grounds
   (Wigner-Eckart: a density-type perturbation has zero matrix element
   between a singlet ground state and a triplet excited state). Edit
   reverted, binary rebuilt back to original md5
   (`8f4b8dcc5ad7dfb88eac350126bf22cf`), no functional change shipped.
   **This exact test (and null result) was already run before**, per
   this log's own "continued 6"-era entry -- confirms nothing new here,
   but importantly that same entry recorded `498542f8`'s commit message
   as explicitly flagging "gap 3.98%->~1%, residual normalization-scale
   mismatch against UHF... not yet resolved" -- i.e. even the FIRST-
   PRINCIPLES-JUSTIFIED fix only partially closed the gap, and what's
   left was already named a *scale* mismatch, not a physical floor.
   That thread was apparently dropped when the investigation later
   pivoted to the Hbar-term hunt (CORR_2P/FACT2) -- picking it back up
   here.

4. **New structural clue**: today's fresh H2O/STO-3G RHF-vs-UHF RCCD
   numbers (current HEAD, e6c42b62, `_wtcheck_baseline`/`_wtcheck`
   test dirs) --
   JJ11: 0.81143686559 vs 0.86408930080 (-6.09%)
   JJ22: 0.00774638191 vs 0.00888750357 (-12.84%)
   JJ33: 0.16664661741 vs 0.17814706433 (-6.46%)
   JJ44: 0.02642713724 vs 0.03041635736 (-13.11%)
   -- split cleanly into TWO groups: {JJ11,JJ33}~6.1-6.5% low,
   {JJ22,JJ44}~12.8-13.1% low, and the second group's gap is close to
   **2x** the first group's in both pairs (12.84/6.09=2.11,
   13.11/6.46=2.03). A genuine physical floor from non-antisymmetric
   amplitudes would not be expected to produce a clean ~2x split
   between two specific components and not others -- this smells like
   a real, findable factor-of-2-type omission specific to whichever
   symmetry/MO sector components 2 and 4 belong to. This directly
   continues (never followed up) the earlier "continued 6" session's
   own flagged next step: identify AO-pair-index 2/4's symmetry/MO
   origin and trace which routine treats it differently RHF vs UHF.
   Not chased further this session (time-boxed single-pass fork task)
   -- concrete next step for whoever picks this up.

**Conclusion**: do NOT record this as a closed/expected-residual
thread. No fix implemented this session (nothing safely derivable in
scope -- per this project's own discipline, not shipping an unverified
coefficient). Recommend next session start from the {1,3}-vs-{2,4}
~2x-ratio pattern above rather than either (a) more Hbar-formula
guessing (already closed, see entries below) or (b) accepting the
non-antisymmetry argument as a blanket explanation (evidence above
argues against it).

No code changes committed -- edit tested and reverted, binary
confirmed restored to original md5. Doc-only commit.

---

## SESSION 2026-08-21 (Ajith, judgment call) -- resolved
the "FLAGGED FOR AJITH" item from the previous entry below. Ajith's
conclusion: **option (b)** is the answer -- rCCD/drCCD RHF/UHF closed-
shell identity is NOT assured in this code, and isn't expected to be.
The 2020 paper's Appendix VII.A already proves rCCD amplitudes are
non-antisymmetric, so the ordinary-CCD guarantee that RHF and UHF
agree for closed shells doesn't transfer to (d)rCCD automatically --
it required a specific, non-obvious spin-adapted combination just for
the ground-state energy (`rcc_cmpeng.F`, already fixed/committed). The
response-property (JJ-matrix) analog of that same subtlety is the most
likely explanation for the remaining 6-13% (RCCD) / 0.3-2.1% (DRCCD)
residual -- i.e. this is plausibly the physically-correct answer, not
a bug to keep chasing toward exact RHF/UHF equality. Option (a) (JJ-
matrix might not be the same object as the paper's omega_k / might need
the R=A+BT alternate route) is NOT ruled out by this and remains a
loose end if the residual ever needs to be explained more precisely,
but is no longer the leading hypothesis.

**Status**: this closes the "keep hunting for a missing Hbar term"
thread (CORR_2P and FACT2 both confirmed non-starters, previous entry).
Not yet done: an actual derivation of the correct spin-adapted RHF/UHF
response-property formula (the response-quantity analog of Eq 46/49/65
that fixed the ground-state energy) has NOT been worked out -- this
would be needed to know exactly what residual size IS expected/correct,
as opposed to just "some nonzero gap is expected." Left as a future
step, not started this session.

---

## SESSION 2026-08-21 (autonomous fork, literature-derivation session,
per Ajith's explicit invitation while busy with nursing school --
MAJOR REFRAME, RESUME HERE) -- went back to the actual source
literature instead of continuing to guess at borrowed DCCSD
coefficients. Result: the whole "naive Hbar is missing a Lagrangian
term, port DCCSD's fix" premise driving the last several sessions was
built on a misidentified analogy. Found the real governing equations,
and they retroactively EXPLAIN (not just re-confirm) why CORR_2P and
FACT2 both failed -- they were never supposed to help this
calculation. No code changes this session; this is a theory-grounding
result, not a fix. Full derivation trail below; TL;DR at the end.

**Citation correction**: `rcc_driver.F`'s own header cites "JCP 139,
104113 (2013)" for rCCD/drCCD. Fetched it (arXiv:1306.6360, Scuseria/
Henderson/Bulik) -- it's actually about the PARTICLE-PARTICLE
(ladder-CCD/pp-RPA) channel, a different connection. The real rCCD/
drCCD/ph-RPA paper it's a follow-up to is Scuseria, Henderson,
Sorensen, JCP 129, 231101 (2008) (arXiv:0810.3433) -- fetched and read
in full. That 2008 paper only covers the GROUND-STATE energy
equivalence (rCCD's T-amplitude Riccati equation <-> RPA's symplectic
eigenvalue problem, T=YX^-1) -- no Lagrangian, no Hbar, no excited-
state/response formalism at all. Not directly useful for the Hbar
question, but establishes the T=YX^-1, R=A+BT machinery used below.

**The actually-relevant paper is Ajith's own** (already used once in
this codebase for the RCCD ground-state energy formula fix,
`rcc_cmpeng.F`): Rishi, Perera, Bartlett, JCP 153, 234101 (2020)
(arXiv:2008.00336). Fetched and read the FULL text (both PDFs
downloaded via WebFetch, converted with `pypdf` since no
`pdftotext`/poppler is installed on this machine -- `python3 -m pip
install --user pypdf` worked fine, no sudo needed). This paper is
*exactly* about deriving the Lagrangian-based Hbar for (d)rCCD's
excited-state/EOM formalism -- precisely the open question.

**What the paper actually says (Section II.B, Eq 14-27, 30-34)**:

1. The paper explicitly confirms the naive `e^-T H e^T`-mirrored form
   for the one-particle Hbar block -- Eq (30)-(31): `H̄i_m = f_im +
   1/2 t^fa_nm<ni||fa>`, `H̄a_e = f_ae - 1/2 t^fa_mn<mn||fe>` -- is
   **"not precisely correct"** for (d)rCCD. This is structurally the
   SAME quadratic Fmi/Fae dressing that FACT2/PDCC_QUAD2/PDCC_QUAD3
   implements. The paper's own numerical Table I shows this
   "Fock-dressing" makes excitation energies WORSE by
   ~0.3-3.9 eV (their `EOM(S)` vs `EOM(Sf)` comparison) -- dressing
   the one-particle block is empirically harmful, not just
   theoretically suspect. **This is an exact, independent match to
   this session's own FACT2=1.0 test result (update 22 above): clean,
   coherent worsening, no exception.** FACT2 should stay OFF -- now
   confirmed correct by both an independent code test AND the source
   literature, not just empirically.

2. The CORRECT one-particle block, per the paper, is the **BARE,
   UNDRESSED Fock matrix**: `H̄a_e = f_ae`, `H̄i_m = f_mi` (Eq 18) --
   exactly what scrnc's Hbar(i,j)/Hbar(a,b) block already is (already
   noted in an earlier session as "PURE BARE-FOCK, zero correlation
   correction"). **scrnc's existing one-particle block was already
   theoretically correct, unmodified, this whole time.**

3. For the two-particle block, the paper's Lambda-equation derivation
   (Eq 17-18) DOES define an Hbar element with an extra term beyond
   the plain ring-CC W: `H̄mb_ej = W^mb_ej + 1/2 t^nj_fb<mn||ef>`.
   This LOOKS like exactly what CORR_2P was trying to port from DCCSD
   (same T^2-quadratic structure). **But this is the wrong block.**
   Eq 17-18 is deriving H̄ AS USED IN THE Λ2 (ground-state Lagrange
   multiplier) EQUATION -- a separate, ground-state-only solve, not
   the excited-state diagonalization block.

4. When the paper writes out the actual excited-state block that gets
   diagonalized to produce excitation energies -- H̄_SS (Eq 34) and
   the explicit matrix element used in `[H̄SS R]` (Eq 26) -- **it has
   NO extra term**: `[H̄SS R]ai,bj = (εe-εm)δaeδim + <ia||bj> +
   t^ae_im<mj||eb>` (Eq 26), and the paper states outright (Eq 27)
   that this equals **exactly `A + BT`**, matching the canonical RPA
   matrix construction (Eq 11) term-for-term, with T being the plain
   ring-CC amplitude and no additional quadratic dressing anywhere.
   **The extra `+1/2 t<mn||ef>` term belongs only to the Λ2/ground-
   state-multiplier solve, never to the H̄_SS block that excitation
   energies (and, presumably, scrnc's JJ-matrix response quantity)
   are actually diagonalized from.**

**Conclusion: CORR_2P was targeting the wrong Hbar block.** Its
correction term structurally resembles a real term from this paper's
own Lagrangian derivation, but that term lives in the ground-state Λ2
equation, not in the H̄_SS block scrnc actually reads/diagonalizes.
Applying it there was never going to help, regardless of coefficient
-- explains why every CORR_2P value tried (-1.5, -10) either did
nothing or made things worse, and retroactively validates leaving it
at 0. **Both of this multi-session thread's "borrowed DCCSD term"
correction attempts (CORR_2P, FACT2) are now understood, not just
observed, to be non-starters for this specific Hbar block** -- the
paper's own math says the correct H̄_SS for (d)rCCD is *exactly* the
bare-Fock + plain-W construction scrnc had *before any of this
thread's correction attempts began*.

**One more striking, unplanned confirmation**: checked
`scrnc_jj_rcc_assemble.F`'s own singlet/triplet weighting (`Wt=0.25D0`
for Imult=1/singlet, `Wt=0.75D0` for Imult=2/triplet, line 81-82) --
this is an EXACT match to the paper's own spin-adapted RHF rCCD
ground-state-energy weighting formula, Eq (65): `E = 1/4(3*(3B*3T) +
1B*1T)`, i.e. 3/4 triplet + 1/4 singlet (the spin-flip-consistent,
degeneracy-weighted form, not the simpler non-degenerate Eq 49 which
splits 1/4:1/4). scrnc's existing weighting is independently confirmed
correct against this paper too, not a suspect.

**What this DOESN'T resolve**: if H̄_SS is provably just A+BT and both
its blocks (bare Fock, plain W) are already exactly what scrnc
computes, then per this paper's own math there should be NO missing
term left to find in the Hbar formula at all -- yet the RHF/UHF
residual (6-13% RCCD / 0.3-2.1% DRCCD) is real and reproducible. Two
honest possibilities, neither chased further this session (needs
either more of scrnc's own architecture understood, or Ajith's
judgment):
(a) scrnc's "JJ-matrix" is not literally the same object as this
    paper's excitation energies ωk (it's a GW/self-energy-pipeline
    response quantity, not benchmarked in this paper at all) -- it
    may need a DIFFERENT construction than H̄_SS-diagonalization
    entirely, e.g. the paper's own "alternate route" (Eq 28-29:
    reconstruct R=A+BT from the converged T amplitudes directly and
    diagonalize R, bypassing an EOM-CC-style Hbar object altogether --
    the paper itself uses this route whenever ring-CCD fails to
    converge, section III/E). Untested whether scrnc's JJ-matrix
    formula is equivalent to this or genuinely different.
(b) rCCD/drCCD's amplitudes are proven NON-antisymmetric (paper's
    Appendix VII.A, Eq 43-49) -- RHF and UHF rCCD ground-state
    energies are NOT automatically identical for closed shells the
    way ordinary CCD's are; matching them at all required a specific,
    non-obvious correct spin-adapted combination (already worked out
    and committed for the ENERGY via `rcc_cmpeng.F`, per the earlier
    "continued 4" session). It is plausible the RESPONSE-property
    (JJ-matrix) analog of this same non-antisymmetry subtlety is what
    the current 6-13%/0.3-2.1% residual actually reflects -- i.e. a
    real, structurally-expected small RHF/UHF difference for this
    specific derived quantity, not a bug to chase to exactly 0. This
    paper doesn't give the response-property formula directly, so
    this can't be confirmed or ruled out from it alone.

**FLAGGED FOR AJITH**: given (a)/(b) above, the next productive step
isn't more Hbar-term guessing (that avenue is now closed on solid
theoretical grounds, not just empirical failure) -- it's pinning down
exactly what scrnc's JJ-matrix quantity is meant to represent
physically and whether an exact RHF/UHF identity is even the right
expectation for it, the same way it wasn't a given for rCCD's
ground-state energy until Eq 46/49/65's specific spin-adapted
combination was worked out. This is a judgment call about what the
GW/Sigma_c pipeline's ground-truth check *should* be measuring, not
a code bug hunt.

No commits with code changes -- this entry (`SCRNC_DEBUG_LOG.md` only)
committed as a documentation-only commit, same pattern as 64f69be6.

---

## SESSION 2026-08-21 (autonomous fork, downstream-of-Hbar audit per
Ajith's explicit direction to stop guessing at more Hbar-formula
correction terms, RESUME HERE) -- audited the whole Hbar-to-JJ-matrix
chain in `aces2/scrnc/` (`scrnc_form_rcc_hbar_inv.F`,
`scrnc_jj_rcc_assemble.F`, `scrnc_form_mubar_s.F`,
`scrnc_form_mubar_00.F`) for plumbing/indexing/formula bugs,
independent of whether the Hbar formula itself is complete. **No new
code bug found** -- but found and resolved a real, previously-
unexplained puzzle: why did the RHF/UHF gap apparently GROW (1.02-1.41%
on 2026-08-19 -> 6-13% now) after MORE fixes landed, not shrink?

**Answer: the two numbers were never measuring the same thing.**
`git log --follow` on `scrnc_form_mubar_s.F` shows its list34/list39
T2-dressing fix (498542f8, 2026-08-19, "gap 3.98%->~1%") was committed
and validated a full day BEFORE the pipeline-ordering Hbar-
contamination bug was even discovered (9843aa74/e48b5cc0/c4169bec, all
2026-08-20). So that "~1%" agreement was necessarily measured against
the OLD, since-proven-corrupted Hbar (contaminated by xringcl's Lambda
solve, per c4169bec's fix) -- not the current, independently-eigenvalue-
validated one. This is the exact "two bugs partially cancel, looks like
a real partial fix" trap this same log already flagged as a lesson
earlier in the 2026-08-19 thread (see the "Lesson for next time" note
in the entry below) -- just not recognized as applying to itself until
now.

**Confirmed the list34/list39 fix is still correctly present and
untouched**: `scrnc_form_mubar_s.F` was never touched by any of the
2026-08-20 pipeline-fix commits (confirmed via `git log --follow`),
and its `LISTT2B = 34` branch for `RCCD` is still active in the current
build -- its underlying physical justification (RCCD's un-
antisymmetrized T2 makes list39's RNABIJ resort degenerate to an exact
copy of list37, silently zeroing CCSD's inherited "2*T37-T39" dressing
formula; list34 is the genuine distinct exchange partner) never
depended on the pipeline bug at all and is independently sound. So this
fix is not suspect -- it's just that its apparent "closes the gap to
~1%" success was, in hindsight, likely inflated by an accidental
cancellation with the (at-the-time-unknown) pipeline bug.

**Rest of the downstream chain reviewed, clean**:
- `scrnc_jj_rcc_assemble.F`'s singlet/triplet `Wt` combination
  (0.25 singlet / 0.75 triplet) matches BOTH the standard closed-shell
  Casida RPA construction (Ajith-confirmed 2026-08-19, "NO, they should
  not decouple") AND `rcc_cmpeng.F`'s independently-validated
  `E_no_flip` energy-formula shape -- no mismatch found between the two
  already-separately-validated 0.25/0.75 usages.
- `scrnc_form_rcc_hbar_inv.F`'s RHF branch correctly switches
  `LISTDIAG` (56 singlet / 54 triplet) with `IMULT`, and already carries
  its own internal `(A^-1)^T = (A^T)^-1` consistency check (pure linear-
  algebra identity, `_DEBUG_LVL0`-gated) -- no asymmetry found.
- `Symmet2` (the missing-final-symmetrization fix from 35cdb853,
  2026-08-15) is confirmed still wired into the RCCD-specific assembly
  path (`scrnc_jj_rcc_assemble.F:229`), not just the older CCSD path
  (`scrnc_jj_cc_lrspn.F:385`).
- `scrnc_form_mubar_00.F` re-confirmed clean per the existing
  2026-08-19 note (T1 read is naturally zero for RCCD's singles-free
  amplitude structure, no bug).

**No code changes made this session** -- per Ajith's direction, did not
retune/re-enable CORR_2P or FACT2, and found no mechanical bug to fix
in this chain. **Net effect: strengthens, doesn't weaken, the case that
the current 6-13% (RCCD) / 0.3-2.1% (DRCCD) residual is a real, honest
number against the now-correctly-validated Hbar** -- not an artifact of
a plumbing bug hiding downstream, and not explained by the apparently-
better-but-baseline-mismatched 2026-08-19 number. **FLAGGED FOR AJITH**:
the remaining lever is a genuine Lagrangian re-derivation for the
missing Hbar terms (2 of 3 DCCSD-precedent terms already tried and
ruled out, the 3rd structurally inapplicable to RCCD's singles-free
amplitudes, see the b27388de/2026-08-21 entry below) -- this is a
physics-derivation task, not something to keep guessing at in code.

---

## SESSION 2026-08-21 (autonomous fork, coordinator's items 2-4 off the
c4169bec resume point, RESUME HERE) -- closed out CORR_2P (already
decided, just needed committing), found and fixed a real dead-wire bug
that made every prior FACT2 test a no-op, ran the REAL FACT2=1.0 test
for the first time (clean negative result, no pole), confirmed UHF is
untouched by all of today's changes, and confirmed the c4169bec
snapshot fix generalizes to DRCCD (both the eigenvalue-level and
JJ-matrix-level checks).

**Item 2a (CORR_2P, already done, just landing it here):** re-scanned
at -1.5 and -10 against the validated c4169bec baseline. JJ11 nearly
insensitive (-6.16%->-6.22% from 0 to -10); JJ22/JJ33/JJ44 all get
WORSE as CORR_2P goes more negative (JJ22 -12.91%->-14.71%, JJ33
-6.49%->-7.21%, JJ44 -13.20%->-13.63%). CORR_2P=0 (no correction) is
still the best value. Left at 0.

**Item 2b (FACT2, the real work this session) -- found a dead wire,
fixed it, then got a clean negative result:** the working tree already
had FACT2 re-enabled (=1.0D0) from a prior (uncommitted) attempt at
re-testing Part 2 against the validated baseline, per the coordinator's
brief. Rebuilt and ran H2O/STO-3G RHF/UHF JJ-matrix with it -- got
numbers **bit-identical to FACT2=0** (JJ11 -6.16%, JJ22 -12.91%, JJ33
-6.49%, JJ44 -13.20%, matching the shipped baseline to 3+ sig figs).
That's not a coincidence-free pass, it's a red flag: traced it and
found `scrnc_form_rcc_hbar.F`'s read of list 91/92 sub-index 9 (Part
2's actual output, written by `PDCC_QUAD2`/`PDCC_QUAD3`) was **commented
out** -- a dead wire. `PDCC_QUAD2`/`3` only ever write list 91/92 (confirmed
by reading the routines directly: `GETLST`/`PUTLST` on 91, nothing
else) and nothing reads it back for RHF, so every prior FACT2 test,
**including the ORIGINAL pre-pipeline-fix blowup scan** (the one that
found the 2 small Hbar eigenvalues crossing zero and JJ(4,4) exploding
1.34->8.86->33.2), may never have actually exercised Part 2's
contribution to Hbar(i,j)/Hbar(a,b) at all. Unresolved what actually
produced that original blowup, if not this -- flagged, not chased
further, see "NOT YET DONE" below.

Re-enabled the read (`scrnc_form_rcc_hbar.F`, guarded `IF (IUHF.EQ.0)`,
mirroring the pre-existing commented-out code exactly) and re-ran
FACT2=1.0 for real. Verified the fix is safe first: with FACT2=0 (no
write), the now-unconditional read still reproduces the exact baseline
bit-for-bit (`UPDMOI` zero-initializes list 91/92 sub-index 9 before
any write, so reading it back unwritten gives clean zeros, not
garbage) -- confirmed empirically, not assumed. Then the real FACT2=1.0
result: **no blowup, no near-zero eigenvalue crossing -- a clean,
coherent WORSENING of all 4 components** (JJ11 -6.16%->-6.62%, JJ22
-12.91%->-13.32%, JJ33 -6.49%->-6.75%, JJ44 -13.20%->-14.24%).
Unambiguous negative result, same shape as CORR_2P's own scan. FACT2
left at 0 -- Part 2, in its current empirically-tuned form, does not
help either component of the Lagrangian-correction thread. Both
`CORR_2P=0` and `FACT2=0` are now the validated, tested-for-real (not
just inert-by-bug) shipped state.

**Item 3 (UHF re-verification, CONFIRMED, not just assumed):** reran
H2O/STO-3G UHF RCCD JJ-matrix with the current binaries (post-c4169bec
AND post today's `scrnc_form_rcc_hbar.F` read-wiring fix, which is
gated `IUHF.EQ.0` and so should be structurally RHF-only). Result: JJ
diagonal **bit-identical** to the pre-existing `h2o_scrnc_rccd_uhf_sto3g_cleanbaseline`
reference (0.86509751290, 0.00876833114, 0.17666365478, 0.03016663846
-- all 4 components to 11 sig figs). UHF is confirmed completely
unaffected by everything done in this thread since c4169bec.

**Item 4 (generalization beyond H2/STO-3G RCCD):**
- **DRCCD, eigenvalue-level:** H2/STO-3G DRCCD, same-run ringcl-vs-scrnc
  cross-check (the same methodology that validated c4169bec for RCCD):
  ringcl's own independent diagonalization = **43.23234947 eV**;
  scrnc's own Hbar diagonalization = **43.23234947 eV** -- bit-identical.
  The snapshot fix generalizes cleanly to DRCCD, not just RCCD.
- **DRCCD, JJ-matrix-level:** H2O/STO-3G DRCCD RHF-vs-UHF (runs already
  present in the tree from before this session, dated within the
  c4169bec-validated window): JJ11 -0.34%, JJ22 -1.93%, JJ33 -0.80%,
  JJ44 -2.14% -- MUCH tighter than RCCD's -6% to -13% band, and still
  coherent (all undershooting, no sign flips, no blowups). Suggests
  DRCCD's RHF/UHF response agreement is intrinsically better-behaved
  than full RCCD's, which is at least physically plausible (dRPA-type
  ring diagrams vs the fuller RCCD Hbar) but not independently confirmed
  this session.
- **DZ (H2/DZ RCCD): attempted, NOT completed.** SCF itself failed
  (`@MINV Singular Matrix` inside `xvscf`, before ever reaching
  ringcc/ringcl/scrnc) -- an SCF/geometry-basis convergence problem
  unrelated to anything under test here (H2 at R=0.742 with DZ may just
  need a different starting geometry/guess). Not chased further this
  session; a real gap in the generalization evidence, not a finding.

**NOT YET DONE, next steps:** (1) what actually produced the ORIGINAL
FACT2 blowup (session "(continued 2)", CORR_2P/Part-1-only + Part-2 at
various strengths) if the read-side dead wire meant Part 2's Hbar(i,j)/
Hbar(a,b) contribution was never live -- either `PDCC_QUAD2`/`3` have an
undocumented side effect beyond list 91/92 (worth a direct read of
those two routines' full bodies, not just their `GETLST`/`PUTLST`
calls), or that particular historical scan was run against a
DIFFERENT, since-superseded code state where the read WAS wired in and
someone commented it out later without re-testing -- git-blame/history
on `scrnc_form_rcc_hbar.F`'s comment-out would settle this quickly; (2)
H2/DZ RCCD generalization still open (SCF convergence issue, not a
scrnc issue); (3) the residual, now cleanly-measured ~6-13% RCCD /
~0.3-2% DRCCD RHF-vs-UHF undershoot itself is still open and needs
either a real Lagrangian re-derivation or an independent (e.g.
finite-field) ground truth, not more scalar tuning -- this is a
genuine physics question for Ajith, not something further code-reading
will resolve; (4) the still-unfixed `RCL_HEFFXC.F` `LISTDUMP`-always-56
bug (finding 2, "SESSION 2026-08-20 (autonomous fork, independent-
reference cross-check..." entry below) remains a latent risk, not
touched this session either.

---

## SESSION 2026-08-20 (autonomous fork, e48b5cc0's fix corrected and
RIGOROUSLY VALIDATED, RESUME HERE) -- per Ajith's own caution ("confirm
at the eigenvalue level before re-tuning"), re-ran the independent-
reference eigenvalue check against e48b5cc0's committed fix and found
it BROKEN for H2 (gave an exact zero Hbar), then found and fixed the
real architectural problem, then CONFIRMED the new fix bit-identical
against ringcl's own independent diagonalization.

**e48b5cc0 was broken**: its fix re-ran `RCC_RNABIJ`+`RCC_GENINT`
fresh from INSIDE `xscrnc`'s own, separate process. Worked (by
coincidence, not by correctness) for H2O/STO-3G but gave scrnc's own
Hbar as an EXACT ZERO MATRIX for H2/STO-3G (`@CHECKSUM: scrnc-raw-
List-56: 0.0 0.0`, confirmed printed "0ZERO MATRIX." directly) even
though list 46 (T2) itself had correct, non-zero content at that point
-- calling ringcc's own W-construction machinery from a DIFFERENT
process than the one that built its list registrations is not reliably
safe (exact mechanism not fully root-caused, dimension/system-size
dependent).

**First replacement attempt also failed, differently**: tried wrapping
`RCL_DRIVER` (inside xringcl, the actual iterative Lambda-equation
solver -- confirmed the real culprit, see previous session entry below)
with a save-before/restore-after pattern, mirroring the existing
(working) save/restore already used around `Rcl_check_rcchbar`'s own
diagnostic in `rcl_hbar.F`. Failed silently: **`RCL_DRIVER` itself
calls `CALL ACES_FIN` + `CALL ACES_EXIT(0)` internally upon completing
the Lambda solve** (confirmed via `grep` on `rcl_driver.F`, same
pattern `RCC_DRIVER.F` uses in xringcc) -- it never returns to its
caller, so "restore after the call" is architecturally impossible.
Confirmed empirically: a checksum placed immediately after `CALL
RCL_DRIVER` never printed in the log, while a checksum placed
immediately before it did.

**Working fix**: new file `aces2/ringcl/rcl_driver_hbarsafe.F`
(`RCL_DRIVER_HBARSAFE`, wired into `rcl_main.F` in place of the direct
`Call Rcl_driver(...)`) snapshots lists 54(-59) to NEW, dedicated list
numbers 154(-159) via `INIPCK`+`GETALL`+`PUTALL` **BEFORE** calling
`RCL_DRIVER` (not after -- sidesteps the no-return problem entirely,
since the snapshot doesn't depend on the call ever returning).
`Scrnc_form_rcc_hbar.F`'s RHF read (`GETALL(...,156)`/`GETALL(...,154)`
instead of the live `56`/`54`) now reads this protected snapshot
instead of whatever `RCL_DRIVER` leaves behind -- only the READ side
changes, the later `LISTDUMP` write-back still targets `56`/`54`
(scrnc's own internal working copy for the rest of its own pipeline,
unaffected by anything outside its own process). `Scrnc_driver.F`'s
broken rebuild-based attempt was fully reverted.

**Rigorous validation (item 1, Ajith's explicit ask -- confirmed at the
eigenvalue level, not just "looks more coherent"):** H2/STO-3G RHF
RCCD, same-run comparison: scrnc's own Hbar eigenvalue = **25.50578519
eV**; ringcl's own independent diagonalization (computed earlier in the
SAME xringcl run) = **25.50578519 eV** -- BIT-IDENTICAL to all 9
printed significant figures. (For reference, the earlier BROKEN
rebuild-based fix gave scrnc 16.60694108 eV for this same system --
clearly wrong, confirming this new snapshot-based fix is a genuinely
different, correct mechanism, not a lucky coincidence.)

**H2O/STO-3G JJ-matrix result, snapshot fix (even better than the
earlier broken rebuild-based attempt's numbers):**

| | baseline (buggy) | broken rebuild-fix (e48b5cc0) | snapshot-fix (this commit) |
|---|---|---|---|
| JJ11 | +9.58% | -8.57% | -6.16% |
| JJ22 | -7.78% | -14.39% | -12.91% |
| JJ33 | +0.46% | -7.64% | -6.49% |
| JJ44 | +2269.19% | -20.71% | **-13.20%** |

All 4 components now sit in an even tighter, more coherent -6% to -13%
undershoot band (vs the broken fix's -8% to -21%) -- consistent with
this being the architecturally-correct fix, not a coincidence that
happened to work for one system.

**NOT yet done (coordinator's items 2-4, continue from here)**: (2)
re-attempt the Lagrangian correction (Part 1/2) against this NOW-
validated baseline -- the four components moving in a coherent,
same-direction band (unlike the old, corrupted-pipeline coefficient
scan's contradictory per-component directions) means a single
coefficient may actually work uniformly this time; (3) UHF re-
verification -- `RCL_DRIVER_HBARSAFE`'s snapshot logic already covers
the UHF branch (lists 154/155/157/158/159), but `Scrnc_form_rcc_hbar.F`
was ONLY updated to read the snapshot for the RHF (`IUHF.EQ.0`) branch
-- UHF still reads directly from `LISTAA=53+ISPIN` (54/55), unchanged
from before this whole session, so UHF's results should be completely
unaffected, but this needs an actual rerun to confirm rather than
assume; (4) H2/DRCCD generalization beyond what's already covered by
the H2/STO-3G RCCD check above.

---

## SESSION 2026-08-20 (autonomous fork, independent-reference cross-check
per Ajith's direction, RESUME HERE) -- MAJOR REDIRECTION: reproduced the
2026-08-14 UHF Hbar-eigenvalue cross-check methodology for RHF and found
scrnc's own Hbar genuinely disagrees with ringcc/ringcl's independent
Hbar for RHF -- but traced WHY, and it is NOT the "missing Lagrangian
term" this whole thread has been chasing. It looks instead like a
**pipeline-ordering bug**: scrnc may be reading the Lambda-equation's
OWN dressed intermediate off lists 54/56, not the plain ground-state
Hbar `RCC_DWMBEJ`/`RCC_GENINT` built, because `xringcl`'s Lambda solver
(`rcl_driver.F`) runs between `xringcc` and `xscrnc` and overwrites
those same lists.

**Setup**: built the RHF counterpart of the existing UHF verify pair
(H2/STO-3G, `h2_scrnc_rccd_uhf_sto3g` + `h2_rpa_verify_rccd_uhf_sto3g`,
2026-08-14) -- `h2_scrnc_rccd_rhf_sto3g_verify` (scrnc's own production
pipeline, `CALC=RCCD,DENSITY=RESPONSE,REF=RHF`, ending in `xscrnc`) and
`h2_rpa_verify_rccd_rhf_sto3g[_noprops]` (pure `excite=EOM(Sf)-RCCD,
eomref=RCCD,REF=RHF`, ending in `xringcl`, no `xscrnc`/no
`DENSITY=RESPONSE`). Confirmed `xvprops` inclusion/exclusion makes zero
difference to the eigenvalues (ruled out as a confound).

**Read `aces2/ringcc/check_rhf_rcc.F` first, per Ajith's lead -- it is
NOT an Hbar/eigenvalue check.** It's a ground-state singlet/triplet
spin-adapted ENERGY consistency check (already-superseded machinery
from the `rcc_cmpeng.F` energy-formula fix in `d80f2742`), unrelated to
Hbar diagonalization. Did not pursue further.

**Key finding 1 -- `RCL_ADDONEH`/`SCRNC_ADDONEH` are byte-identical**
(diffed directly, only the subroutine name/line-wrap differs) -- ruled
out the Fock-diagonal-add step itself as a source of any RHF-specific
discrepancy.

**Key finding 2 -- `RCL_HEFFXC.F`'s own `LISTDUMP` is hardcoded to 56
regardless of `IMULT`** (unlike `scrnc_form_rcc_hbar.F`, which THIS
session's earlier work -- committed before this thread even started,
2026-08-19 -- already fixed to route `IMULT=2` to its own list 54). For
RHF RCCD, `Rcl_check_rcchbar` sets `Nmult=2` and calls `Rcl_heffxc`
twice (singlet then triplet), and the triplet call's `PUTALL` briefly
clobbers list 56 with the triplet's own Fock-added Hbar -- harmless
here only because `Rcl_drvrpa` diagonalizes immediately after each call
inside the same loop, before the next iteration's overwrite. Still a
real, unfixed bug in `RCL_HEFFXC.F` (should mirror scrnc's own
IMULT-aware `LISTDUMP` fix) -- flagged, not fixed this session (doesn't
explain the discrepancy found below, but is a latent risk for any other
consumer that reads list 56 between the two `Rcl_heffxc` calls).

**Key finding 3 (the real one) -- `rcl_dwmbej_r.F`'s own header comment
says exactly what it builds: `W(mbej) = 2 W'(mbej) - <mb||ej>`, "the
ring intermediates... appropriate for the CCD LAMBDA EQUATIONS,"
EXPLICITLY DIFFERENT from `W'(mbej)`, "the intermediate used in the
energy calculation."** `Rcl_dwmbej_r`/`_u` run unconditionally
whenever `.NOT. Solve_4lambda` inside `rcl_hbar.F` (i.e. for the
"verify" EOM(Sf) pathway too) and overwrite lists 54/56 with this
Lambda-specific combination BEFORE `Rcl_check_rcchbar`'s own diagnostic
runs -- meaning the 2026-08-14 "bit-identical" UHF result validated
that scrnc's and ringcl's DIAGONALIZATION ROUTINES agree when fed
IDENTICAL (already Lambda-transformed) input, not that the raw,
untransformed ground-state Hbar is itself what either program
diagonalizes. This alone isn't the bug (both UHF and RHF go through
the same transform) -- but combined with finding 4 below, it's the key
that unlocks the real mechanism.

**Key finding 4 -- `rcl_driver.F` (the ACTUAL Lambda-equation ITERATIVE
solver) runs AFTER `Rcl_hbar`/`Rcl_check_rcchbar`'s diagnostic, and
FURTHER modifies lists 54/56**, confirmed via a same-run checksum
comparison on H2O/STO-3G: `Rcl_dwmbej_r`'s own pre-iteration checksum
of list 56 (added this session, matches the diagnostic's own restored
value) vs. list 56's checksum printed later in the SAME run (after
`rcl_driver.F`'s "Ring coupled cluster doubles (rCCD) lambda equations
are solved" message) differ measurably, and scrnc's own EVENTUAL read
of list 56 (in the separate `xscrnc` process, after `xringcl` has fully
exited) differs again, roughly 2x in checksum norm from the
pre-iteration value in one specific comparison. `rcl_driver.F` only
runs a REAL iterative Lambda solve when `Solve_4lambda=.True.`, which
`DENSITY=RESPONSE` triggers regardless of the `excite=` keyword
(confirmed: even the "clean" EOM(Sf)-RCCD ZMAT, with `DENSITY=RESPONSE`
added just to let `xscrnc` run per its own hard requirement, printed
"lambda equations are solved").

**Decisive same-run comparison** (H2/STO-3G RHF RCCD,
`h2_rpa_verify_rccd_rhf_sto3g_withscrnc`: `excite=EOM(Sf)-RCCD,
eomref=RCCD,DENSITY=RESPONSE,REF=RHF`, full 10-step pipeline through
`xscrnc`, WITH this session's Part-1 Lagrangian correction active
(CORR_2P=-1.5, the currently-committed dabea903/aa22a7ea state) --
same ZMAT/keywords used for BOTH diagnostics, zero pipeline-config
confound):
- ringcl's own independent Hbar diagonalization (printed EARLIER in
  the SAME `xringcl` run, BEFORE `rcl_driver.F`'s Lambda solve):
  **26.22164379 eV (singlet), 16.66249500 eV (triplet)**
- scrnc's own Hbar diagonalization (`xscrnc`, reading lists 54/56
  AFTER `xringcl` has fully exited, i.e. AFTER the Lambda solve ran):
  **19.47037548 eV** (only one value -- scrnc's `Nmult=1` architecture
  only ever reports the singlet channel)

19.47 matches NEITHER of ringcl's own numbers -- and the relevant
comparison (scrnc reports singlet only) is 19.47 vs 26.22, a ~26% gap,
LARGER than any of the JJ-matrix percentage gaps this whole thread has
been chasing. A separate, cleaner run with NO `DENSITY=RESPONSE` (so no
Lambda solve at all, `h2_rpa_verify_rccd_rhf_sto3g[_noprops]`) gives
ringcl's truly pre-Lambda-solve reference: **26.69888286 eV (singlet),
17.44865788 eV (triplet)** -- close to but not identical to the
same-run 26.22/16.66 (confirming the Lambda solve DOES perturb list
54/56 measurably even before scrnc gets to read it, consistent with
finding 4).

**Working hypothesis, NOT yet proven, needs Ajith's physics read or
further tracing**: scrnc's Hbar-consuming machinery may be structurally
wrong to read lists 54/56 AFTER `xringcl` at all, in a pipeline where
`Solve_4lambda` triggers -- it may need the PRE-Lambda-solve
`RCC_GENINT`/`RCC_DWMBEJ` output specifically (the true ground-state
`W'(mbej)`), but by pipeline construction
(`...,xringcc,xringcl,xscrnc`) can only ever see whatever `xringcl`'s
Lambda solve leaves behind. **Why UHF apparently didn't show this
problem in the 2026-08-14 cross-check is unexplained** -- either UHF's
own Lambda-equation solve happens to leave list 54/55 unchanged/
equivalent to the pre-solve Hbar for this specific method (a real
possibility for some CC variants, not verified), or the 2026-08-14 UHF
"bit-identical" match was itself not as clean a test as assumed (worth
re-verifying with the SAME same-run-comparison technique used here for
RHF, not yet done for UHF this session).

**This REDIRECTS the whole Lagrangian-correction thread** -- confirmed,
not just hypothesized: implemented the fix (`Scrnc_driver.F`, guarded
`IF (IUHF.EQ.0 .AND. (RCCD.OR.DRCCD))`, right before `Spinad56`):
re-run `RCC_RNABIJ`+`RCC_GENINT` fresh, rebuilding lists 54/56 from the
converged T2 still on list 46 (untouched by the Lambda solve, which
only touches T with T already held fixed), bypassing whatever
`xringcl` left behind. Both routines already link into `xscrnc` fine
(this session's earlier `ringcc:vcc` addition to `ACES_LINK_LISTS`
covers it).

**H2O/STO-3G RHF-vs-UHF JJ-matrix result, pipeline fix ALONE (Part 1/2
Lagrangian correction fully off, CORR_2P=0):**

| | baseline (old, buggy pipeline) | pipeline-fix-only |
|---|---|---|
| JJ11 | +9.58% | -8.57% |
| JJ22 | -7.78% | -14.39% |
| JJ33 | +0.46% | -7.64% |
| JJ44 | **+2269.19%** (23.7x off) | **-20.71%** |

**JJ44's catastrophic 23.7x baseline error collapses to a normal,
sane -20.71%.** This is the headline result: the near-pole/negative-
eigenvalue instability chased via CORR_2P tuning in the PRIOR session
entry (below) was very likely primarily a symptom of scrnc reading the
WRONG (Lambda-solve-contaminated) Hbar all along, not a genuine
missing-Lagrangian-term physics issue. All 4 JJ components now
undershoot UHF by a moderate, coherent ~8-20% -- a plausible real
residual gap, nothing like baseline's wildly inconsistent sign pattern.

**Tested combining this fix with the OLD Part-1 correction
(CORR_2P=-1.5) on top of it: essentially no additional effect** (JJ11
0.79092 -> 0.79126, ~0.04% shift, similarly negligible on the other 3
components). Root cause: `RCC_LAGRANGIAN_FIX` (Part 1) runs inside
`xringcc`, modifying list 56 there -- but `Scrnc_driver.F`'s NEW fresh
`RCC_GENINT` rebuild runs AFTER that, in `xscrnc`, and OVERWRITES list
56 from scratch, discarding whatever Part 1 had written. **Part 1 is
now functionally inert with this fix in place.** Set `CORR_2P` back to
0 in `rcc_lagrangian_fix.F` (documented as superseded, not deleted --
the machinery may still be useful for a future correction applied
AFTER the pipeline-ordering fix rather than before/underneath it).

**SHIPPED**: the pipeline-ordering fix in `Scrnc_driver.F`, CORR_2P=0
(Part 1/2 Lagrangian correction off). This is a clear, unambiguous,
large improvement over the previous best state (aa22a7ea) and doesn't
carry aa22a7ea's own open problems (JJ(4,4) near-pole instability,
component-dependent CORR_2P tradeoff) at all.

**NOT yet done, next steps**: (1) the remaining ~8-20% undershoot is
now a much cleaner, single-signed residual -- worth a FRESH
investigation (not blind coefficient tuning on the old, now-obsolete
CORR_2P machinery) once this checkpoint is stable; (2) reproduce the
same-run ringcl-vs-scrnc comparison for UHF to confirm/rule out whether
UHF has the identical underlying pipeline issue, just coincidentally
invisible in the JJ-matrix metric so far (the 2026-08-14 UHF
cross-check predates this fix and was not re-verified with it); (3)
re-verify on H2/STO-3G (the originally-flagged smaller-gap case) and
DRCCD, not just H2O/RCCD, to confirm this generalizes -- not done this
session, deliberately, to get this large a result documented and
committed before spending more budget; (4) the still-real, unfixed
`RCL_HEFFXC.F` `LISTDUMP`-always-56 bug (finding 2 above) is worth
fixing on its own even though it didn't explain this particular
discrepancy.

---

## SESSION 2026-08-20 (autonomous fork, continued past dabea903, per
coordinator instruction) -- DIAGNOSED the Part-2/JJ(4,4) blowup
mechanism (near-pole crossing, not a bug), scanned Part 1's CORR_2P
coefficient (0 to -12), and found a decisive negative result: no
single scalar CORR_2P improves all 3 "well-behaved" JJ components at
once. Full detail:

**Item 1 (JJ(4,4) blowup, resolved -- not a bug):** `scrnc_form_
rcc_hbar_inv.F` already has a `_DEBUG_LVL0`-gated real-eigenvalue dump
of scrnc's own Hbar(ai,bj) matrix (10x10 for H2O/STO-3G's singlet
channel) right before it's explicitly inverted via `MINV` -- this ties
DIRECTLY into the older, already-documented "2 spurious negative Hbar
eigenvalues" thread (see "SESSION 2026-08-20 (continued)"/"(continued
2)" further below, from BEFORE this whole DCC-precedent investigation
started: RHF's baseline ring-CC Hbar has 2 negative eigenvalues,
-0.88/-0.63, that UHF's all-positive spectrum has no trace of). Traced
the SAME 2 eigenvalues across correction strength:
- baseline (no correction): -0.88, -0.63
- Part 1 only (CORR_2P=-1.5, the committed dabea903 state): -0.633,
  -0.540 -- still negative, but shrinking toward zero
- Part 1 + Part 2 (both, full strength): +0.122, +0.256 -- crossed
  zero, now positive, matching the SIGN UHF's reference spectrum says
  is correct, but landing very close to zero
- Part 1 + Part 2 (Part 2 at half strength, FACT2=0.5): -0.255, -0.142
  -- confirms the crossing is LINEAR in Part-2 strength, not a jump

A near-zero eigenvalue in a matrix being explicitly inverted for a
static response blows the inverse up in exactly the observed way (the
Part-1-only inverse's rows/cols 9-10, corresponding to these 2 small
eigenvalues, already show anomalous -43/-50 values against an O(0.01-1)
scale everywhere else) -- textbook near-pole response-theory
sensitivity (chi ~ 1/(Omega-Omega_pole) for a static Omega=0 property),
not evidence of a Part-2 implementation defect. JJ(4,4) itself tracks
this exactly: 1.34 (Part 1 only) -> 8.86 (Part 2 half strength) -> 33.2
(Part 2 full strength), growing monotonically worse as the eigenvalues
approach the crossing. **Conclusion: both corrections push these 2
eigenvalues in the physically-correct direction (toward UHF's all-
positive spectrum); the blowup is the expected numerical signature of
passing near/through a genuine near-degenerate excitation channel in
this specific system+basis, not proof Part 2 is wrong.** Whether a
near-zero static response is itself physically plausible here (i.e.
whether STO-3G/RHF/RCCD is close to an actual instability threshold for
this mode) is a genuine physics question, not resolved by code-reading
or further coefficient scanning -- flagged for Ajith if it becomes
relevant again, but NOT a reason to distrust Part 2's direction.

**Item 2 (tune CORR_2P, decisive NEGATIVE result):** scanned CORR_2P in
{0,-1,-1.5,-2,-3,-6,-10,-12} with Part 2 OFF (clear of the pole above).
JJ(1,1)'s RHF/UHF gap shrinks monotonically and crosses zero near
CORR_2P~-10.5 (0:+9.58% -> -1:+8.12% -> -2:+6.83% -> -3:+5.69% ->
-6:+2.91% -> -10:+0.26% -> -12:-0.75%). But JJ(2,2) and JJ(3,3)'s gaps
get monotonically WORSE over the exact same range (JJ(2,2): 0:-7.78% ->
-1:-11.26% -> -2:-13.09% -> -3:-14.14% -> -6:-15.43% -> -10:-15.85% ->
-12:-15.91%; JJ(3,3): 0:+0.46% -> -1:-1.47% -> ... -> -12:-6.33%) --
CORR_2P=0 (no correction at all) is their own individual best value.
**No single scalar CORR_2P is a strict improvement across all 3
components** -- this means Part 1's functional form (one uniformly-
scaled self-referential T*W contraction, same shape for every matrix
element) is structurally incomplete, not just wrongly-scaled. Further
blind coefficient scanning will not fix this; it needs either a real
re-derivation of rCCD's own Lagrangian H-bar (not a DCC-analogy guess)
or a genuinely different functional form (e.g. matrix-element- or
symmetry-block-dependent weighting, since DCC's own formula keeps a
disconnected/exchange term that rCCD's Lagrangian drops -- see the
original task framing). CORR_2P is LEFT AT -1.5 (unchanged from its
original guess) purely because no other single value is a strict
improvement -- not because -1.5 is validated.

**Item 3 (does UHF need the same correction? -- partially answered by
existing evidence, not newly tested):** UHF's own baseline (uncorrected)
ring-CC Hbar spectrum, already documented in the "(continued)" session
below, has NO negative eigenvalues at all (~20, all positive) --
unlike RHF's 2 negative ones. This is evidence that UHF's naive Hbar
does NOT show the same instability RHF's does, which is at least
consistent with (though doesn't strictly prove) UHF being a reasonable
comparison target as-is. Did NOT attempt a full UHF-side port of
RCC_LAGRANGIAN_FIX this session (UHF needs all 6 spin cases -- AAAA/
BBBB/ABAB/BABA for the two-particle piece, more list-number bookkeeping
than the RHF singlet/triplet pair -- a substantial separate effort, not
a quick extension of the existing RETURN-on-IUHF.NE.0 guard). Whether
UHF's ABSOLUTE JJ values are fully correct (as opposed to merely self-
consistent/non-unstable) is a genuinely open question this session did
not resolve -- would need an independent reference (finite-field
numerical check, literature value, or higher-level theory), not just
an eigenvalue-sign check.

**Working tree**: left in a state functionally IDENTICAL to the
committed dabea903 checkpoint (CORR_2P=-1.5, Part 2 disabled) -- only
comments/documentation changed, plus a small refactor (FACT2 scale
variable with an `IF(FACT2.NE.0)` guard, replacing hard-commented-out
Part-2 calls, so re-enabling/scaling Part 2 for future testing is a
one-line change instead of un-commenting code). See
`rcc_lagrangian_fix.F`'s own header for the fully updated STATUS block.

**NEXT STEPS, concretely, superseding the dabea903 list:**
1. Either re-derive rCCD's Lagrangian H-bar properly (the actual ask
   from the start of this thread) or find/construct a genuinely
   different functional form for the two-particle correction -- a
   single global scalar on the current self-referential contraction is
   now DISPROVEN as sufficient (item 2 above).
2. Once a better two-particle form exists, re-examine whether the same
   near-pole sensitivity (item 1) still shows up, and if so, whether
   that's an acceptable/expected feature of this specific system+basis
   or a sign that basis (STO-3G) or geometry needs revisiting for this
   particular test case.
3. Attempt the full UHF-side port (6 spin cases) IF the RHF-vs-UHF
   comparison target itself starts looking suspect -- not urgent given
   item 3's existing (partial) evidence UHF's spectrum is already
   well-behaved.
4. Re-verify on H2/STO-3G (the originally-flagged smaller ~1.2% gap
   case) and DRCCD once a functional-form fix exists worth testing
   broadly -- not done this session, deliberately, since H2O's richer
   4-component JJ matrix was more diagnostic for finding the component-
   dependent tradeoff in item 2.

---

## SESSION 2026-08-20 (autonomous fork, RCCD Lagrangian-Hbar port) --
IMPLEMENTED + PARTIALLY VALIDATED: ported the DCCSD Lagrangian-derived
Hbar correction pattern (see update 9/10/11 in the memory safety net,
`reference_scrnc_debug_log_location.md`, for the full derivation and
build-system/Fortran-implicit-typing debugging trail) into a new file,
`aces2/ringcc/rcc_lagrangian_fix.F` (`RCC_LAGRANGIAN_FIX`), called from
`rcc_driver.F` right after the amplitude iterations converge. Two
pieces: (1) two-particle I(mb,ej)=t(mn,ef)<nj||fb> correction subtracted
from list 56 (singlet)/54 (triplet) Hbar, coefficient CORR_2P=-1.5
(empirical guess, not derived); (2) one-particle Fmi/Fae quadratic
correction via the generic `PDCC_QUAD2`/`PDCC_QUAD3` routines, wired
into `scrnc_form_rcc_hbar.F`'s previously-always-zero one-particle
block (confirmed: RHF RCCD's Hbar(i,j)/Hbar(a,b) was PURE bare-Fock,
zero correlation contribution at all before this).

**RESULT (H2O/STO-3G RHF vs UHF "RPA JJ matrix" ground truth,
`/blue/ufhpc/perera/scrnc_test/h2o_scrnc_rccd_{rhf,uhf}_sto3g_regroup`
baseline, `..._lagrangianfix`/`..._lfix_part1only` corrected runs):**

| | JJ(1,1) (UHF=0.86509751290) | JJ(4,4) (UHF=0.03016663846) |
|---|---|---|
| baseline (uncorrected) | 0.94799229627 (9.58% off) | 0.71470364552 (23.7x off) |
| Part 1 only (two-particle) | 0.92958256128 (7.45% off) | 1.33723905424 (44x off) |
| Part 1 + Part 2 (both) | 0.89764317520 (3.76% off) | 33.22569814971 (1100x off) |

Part 1 alone is a real, if partial, improvement on the "mixed-symmetry"
components (1-3) and only mildly worsens the pre-existing JJ(4,4)
anomaly. Adding Part 2 improves the mixed components further but makes
JJ(4,4) blow up catastrophically -- Part 2 is the culprit, most likely
because JJ(4,4)'s AO-pair symmetry sector already sits close to some
near-degenerate/small-denominator structure (it's exactly symmetry-
decoupled from components 1-3 in both references, per the "SESSION
2026-08-20" entry below) and Part 2's one-particle shift pushes an
RPA-type response denominator through a near-zero crossing for that
sector specifically. NOT investigated further this session (deep
symmetry-sector-specific debugging, likely needs its own top-down
debug-print bisection per this project's established method).

**SHIPPED STATE (committed): Part 1 active, Part 2 DISABLED** (commented
out in both `rcc_lagrangian_fix.F`'s own PDCC_QUAD2/3 calls and
`scrnc_form_rcc_hbar.F`'s corresponding list-91/92 read, which reverts
to the original all-zero one-particle block) -- this is the verified,
non-catastrophic, net-improvement state. CORR_2P=-1.5 is UNCHANGED from
its original derivation guess and NOT tuned -- Part 1 alone only
partially closes the gap (9.58%->7.45%), so it is likely still not the
exactly correct coefficient. See `rcc_lagrangian_fix.F`'s own header
comment for the full derivation, build-system notes (had to add
`ringcc`+`vcc` to `Makefiles/GNUmakefile`'s `ACES_LINK_LISTS` scrnc
entry; had to `UPDMOI`-register list 91/92 sub-index 9 before Part 2's
PDCC_QUAD2/3 can PUTLST to it, since no pure ring-CC pipeline step does
this by default), and a real Fortran implicit-typing bug found+fixed
along the way (this whole class of ringcc routine needs `IMPLICIT
INTEGER (A-Z)`, not `IMPLICIT DOUBLE PRECISION` -- see
[[feedback_aces_fortran_implicit_typing]]).

**NEXT STEPS, concretely:**
1. Investigate JJ(4,4)'s Part-2-triggered blowup specifically -- likely
   a near-zero RPA response denominator in that symmetry sector, made
   worse by the one-particle Fmi/Fae shift. Top-down debug-print
   bisection (this project's established method) on that sector's own
   eigenvalue/denominator structure before Hbar is fed into the CI-like
   solve.
2. Re-derive or empirically scan CORR_2P (the two-particle coefficient)
   more carefully -- -1.5 was a structural-analogy guess from DCC's own
   bookkeeping, not verified against rCCD's own Lagrangian. Given Part 1
   alone only gets JJ(1,1) from 9.58% to 7.45% (not fully closed), try a
   few other coefficient values (e.g. -1.0, -2.0, -1.0 exactly matching
   DCC's OWN net formula without RCC's pre-existing +0.5 term folded in)
   and see which minimizes the RHF/UHF gap on components 1-3.
3. Once Part 2 is fixed/scoped correctly, re-verify on H2/STO-3G too
   (the smaller, originally-flagged ~1.2% gap case) and on DRCCD, not
   just RCCD.
4. Consider whether UHF's own Hbar (built via the SAME
   `RCC_DWMBEJ`/naive-CCSD-mirroring formula, just the UHF spin-case
   branches) has the SAME missing-Lagrangian-term problem -- if so,
   "matching UHF" was never fully correct as a validation target, and
   both references may need the correction (this was flagged as an open
   question, not yet tested empirically this session -- would need
   extending `RCC_LAGRANGIAN_FIX`'s UHF-branch early-return to also
   apply the analogous UHF spin-case correction and comparing corrected-
   RHF vs corrected-UHF instead).

## DEAD ENDS -- DO NOT REVISIT (RCCD thread specifically)

- **`DENSITY=RESPONSE` in RCCD's own pipeline**: proven bit-identical
  inert (2026-08-20) -- only matters for whether `xdens` runs, and
  RCCD's manual pipeline never calls `xdens`. Never add it back to test
  RCCD; `scrnc_set_vars.F`'s hard gate on it is stale gatekeeping.
- **`xvprops` in RCCD's own pipeline**: never needed -- scrnc responds
  to AO-pair density perturbations, not dipole integrals. Only the
  separate `PROPS=EOM_NLO`/`xvcceh` CCSD-reference pipeline needs it.
- **Singlet-only response channel (`Nmult=1`, killing Imult=2/triplet)**:
  ALREADY the shipped architecture since `498542f8`. Real and correct
  (Wigner-Eckart: spin-scalar density perturbation has zero triplet
  overlap), but does NOT close the gap -- H2 residual ~1%, H2O residual
  9.58%/23.7x. Don't re-derive or re-test this hypothesis again.
- **Old "~0.03% RHF/UHF agreement" milestone claim** (from
  `reference_scrnc_debug_log_location.md`, "update 5"): unreproducible
  after isolating every variable (2026-08-19 "continued 3") -- treat as
  false/unverified, not a real prior result to chase back to.
- **`ringcl/dxrcl_dwmbej.F` triplet-branch list-routing copy-paste fix**
  (`NMSYTA`/`LISTTB` mismatch, 2026-08-19 "continued 2"): a real,
  legitimate bug and correctly fixed, but fixing it alone made the
  RHF/UHF gap WORSE (0.933%->3.98%), not better -- it is not, by
  itself, a source of the response gap. Don't expect similar list-
  routing typos elsewhere to be the whole answer.
- **Apparent "3.98%->1.42%" improvement from fixing only
  `scrnc_rcc_hbar_read.F`'s hardcoded list-56 read** (2026-08-19
  "continued 5"): was a MEASUREMENT ARTIFACT -- `scrnc_form_rcc_hbar.F`
  had a second, canceling bug in its Hbar write-back (also hardcoded to
  list 56). Fixing both together restores the original 3.98%. Lesson:
  when a read-side and write-side of the same list pairing both look
  suspicious, fix and verify both together, not one at a time -- a
  lone-side fix can look like progress while actually just uncovering
  a second, compensating bug.
- **List34 as a "singlet-channel exchange partner" to T(list37) in
  `rcc_dwmbej.F`'s RHF Hbar construction** (2026-08-20 "continued 3"):
  tried mirroring `vcc/dwmbej.F`'s `2*T37-T39` pattern with list34 in
  place of the (trivially-equal-to-37) list39 -- made the Hbar's
  negative eigenvalues WORSE, not better. Traced via `RCC_RNABIJ`'s own
  source: list34 is built from a DIFFERENT reference list (44,
  `SPCAS2='AAAA'`, same-spin) via `RCC_SST003`, not a resort of
  list37/46 the way list39 genuinely is -- it's ring-CC's TRIPLET
  channel's own primary amplitude (confirmed: `rcc_dwmbej.F`'s own
  IMULT=2 branch reads it as such), not a singlet-channel correction
  term. Don't reuse list34 this way again without a real derivation
  behind it -- also casts doubt on whether `Scrnc_form_mubar_s.F`'s
  earlier, already-shipped `LISTT2B=34` fix (`498542f8`) rests on solid
  ground, though it did produce a real measured improvement at the
  time and has NOT been re-examined/reverted-and-retested.

---

## SESSION 2026-08-20 (continued 4) -- checked the Rishi/Perera/Bartlett JCP2020 paper directly (arXiv:2008.00336, same reference already used for the ground-state energy formula in d80f2742) for the RHF Hbar/W-intermediate equations -- paper does NOT give an explicit closed-shell rCCD form (explicitly says the naive CCSD H-bar formula "is not precisely correct... for (d)rCCD" without spelling out the alternative); every code-level candidate for the ~1.86-1.92x gap that could be DIRECTLY compared (CNTRCT itself, the C4 scale, the "missing" T2 combination) is now RULED OUT or shown to be a non-issue -- genuinely stuck at a point needing either Ajith's own derivation/recollection or a pilot-code cross-check, not more code tracing (RESUME HERE)

**Fetched the paper** (`https://arxiv.org/html/2008.00336`, via `WebFetch`,
two targeted passes after an initial pass was too lossy/summarized to
trust for exact coefficients) and read Section II.2 ("Towards excited
states through the formulation of Lambda equations and effective
Hamiltonian (H-bar) for ring-CCD") plus the Eq. 49 ground-state-energy
section already used as precedent.

**What the paper actually says, verbatim as transcribed**:
- Eq. 34 (general spin-orbital CCSD-style W-bar):
  `W-bar_iajb = <ia||jb> + Sum_{m,e} t(jm,eb) <mi||eb>` -- this is the
  SAME general ring/W(mbej)-type intermediate structure as
  `rcc_dwmbej.F`'s own header comment already states it's implementing.
- Eq. 30-31 (the naive H-bar = e^-T H e^T one-particle pieces) come
  with an explicit caveat immediately after: **"But as we have seen
  from our discussion, it is not precisely correct to assume that form
  for (d)rCCD. Instead, the expressions for H-bar elements have to be
  obtained from the corresponding Lagrangian equations."** The paper
  does NOT then spell out the corrected closed-shell/RHF-spin-adapted
  W-bar formula explicitly anywhere located -- only the general
  spin-orbital form (Eq. 34) and the separate ground-state-energy
  spin-adapted bare-integral B-matrices (Eq. 49, singlet
  `2<Ij||Ab>-<Ij||bA>` / triplet `-<Ij||bA>`, already confirmed
  correctly implemented in `rcc_dwmbej.F`'s existing bare-integral
  combination, lines 316-329, verified in "continued 6" 2026-08-19).
- Confirms rCCD needs the BARE (undressed) Fock matrix, not a
  similarity-transformed one (`F-bar = F`) -- matches
  `scrnc_form_rcc_hbar.F`'s existing construction (reads raw
  `SCFEVALA`/`SCFEVALB`, already verified correct via the exact-to-11-
  digits Fock-gap check in "continued 2" above).

**Fallback (per instruction): derive from the spin-orbital equation,
compare to code, rule things out systematically.** Since the paper's
Eq. 34 W-bar formula is structurally the same ring/W(mbej) intermediate
`vcc/dwmbej.F` (CCSD, proven-correct) already implements, re-examined
whether ring-CC's version is missing something CCSD's has:
1. **`CNTRCT` itself** (the actual routine performing the T*W
   contraction): `ringcc/cntrct.f` vs `vcc/cntrct.f` diffed directly,
   byte-for-byte -- ONLY whitespace and one commented-out debug
   checksum line differ. The contraction logic is IDENTICAL. Ruled
   out as a source of any scale/coefficient difference.
2. **`C4` scale argument**: both `rcc_genint.F` and `vcc/genint.F` set
   `C4 = -0.50D0` for the standard Hbar-building call path. Identical.
   Ruled out.
3. **The "missing T2 exchange combination"** (this session's earlier
   attempted fix, already reverted): re-examined now that list34 is
   known to be the wrong quantity (previous entry). The STRUCTURALLY
   correct analog of CCSD's list39 -- built the same way, via
   `RCC_RNABIJ`'s own `TYPE='T'` resort of the SAME reference list
   (46) -- is list39 itself, and list39 is ALREADY proven to equal
   list37 exactly for RCCD (both from the 2026-08-19 empirical
   checksum finding and from tracing `RCC_RNABIJ` directly this
   session). This means "no T2 combination" (what the code currently
   does) and "2*T37-T39" (the CCSD-analogous formula, IF you use the
   structurally-correct list39) are the SAME OPERATION for RCCD --
   there never was a missing combination at this step. My earlier
   "continued 3" fix attempt conflated superficial code-pattern
   similarity (both CCSD and RCC do a "2T-T_X" pattern SOMEWHERE) with
   an actual missing physics term, and picked the wrong X (list34,
   which is a genuinely different quantity -- the triplet channel's
   own amplitude, see previous entry) because the RIGHT X (list39)
   trivially reduces to a no-op and so isn't "missing" at all.

**Where this leaves things**: every piece of `rcc_dwmbej.F`'s IMULT=1
construction that could be DIRECTLY compared against `vcc/dwmbej.F`'s
proven-correct analog -- the bare-integral B-matrix (Eq. 49, matches),
the CNTRCT routine (byte-identical), the C4 scale (identical), and the
T2-combination step (provably a no-op either way for RCCD) -- checks
out clean. The paper doesn't supply the missing closed-shell rCCD-
specific H-bar formula to compare against for anything deeper (its own
explicit caveat is that the naive CCSD-style formula, which is what
all of the above comparisons have been implicitly assuming transfers
directly, "is not precisely correct" for (d)rCCD -- the paper just
doesn't say what to replace it WITH). This suggests the actual missing
piece may not be inside `rcc_dwmbej.F` at all -- candidates not yet
examined: whether RCCD's own converged T2 amplitudes (upstream of
`rcc_dwmbej.F` entirely, from whatever solves the ground-state RCCD
amplitude equations) are themselves the right ones for this
excited-state/Hbar construction specifically (the paper's own
caveat about needing "Lagrangian equations... instead" suggests the
one-particle Hbar pieces, or possibly the amplitudes feeding them,
need Lambda/Lagrangian-derived quantities that plain ground-state T2
doesn't supply), or something in how `scrnc_rcc_hbar_read.F`/
`Scrnc_form_rcc_hbar_inv.F` consume the written Hbar downstream
(already partially checked in earlier sessions but not re-verified
against this specific lead).

**FLAGGED for Ajith, now at the limit of profitable code-only
tracing**: this needs either (a) Ajith's own derivation/recollection
of what the Lagrangian-based, non-naive rCCD H-bar formula actually
is (the paper explicitly says one is needed but doesn't give it), or
(b) building an independent, minimal pilot-code cross-check of the
RHF RCCD Hbar diagonal for this exact H2O/STO-3G system (matching this
thread's own established "Pilot-Code Debugging" precedent -- see
memory) against a hand-coded or literature-verified reference value,
rather than continuing to compare `rcc_dwmbej.F` against CCSD's
`dwmbej.F` term-by-term -- that comparison is now exhausted for
everything directly comparable and hasn't found the gap.

**Housekeeping**: no code changes this entry (pure paper-reading +
source comparison, no rebuilds needed). Working tree unchanged from
the entry below (still just the pre-Fock debug print in
`scrnc_form_rcc_hbar.F` and the RHF-RCCD-trap test override in
`scrnc_set_vars.F`, both uncommitted, both still needed, neither
shipped).

---

## SESSION 2026-08-20 (continued 3) -- Ajith's call ("it's (b), the ratio IS the bug") acted on: implemented the vcc/dwmbej.F-analogous `2*T(37)-T(34)` combination in rcc_dwmbej.F's RHF/IMULT=1 branch -- made the negative eigenvalues WORSE, not better; traced WHY via RCC_RNABIJ's own source and RULED OUT list34 as any kind of "singlet exchange partner" -- it's ring-CC's TRIPLET channel's own primary amplitude (built from a genuinely different reference list, 44/AAAA-type), not a permuted view of list37/46 the way CCSD's list39 is. FIX REVERTED. Real question is now sharper but still open: does RHF RCCD's Hbar ring-term need an exchange-type correction to T(37) AT ALL, or was applying CCSD's "2*T-T" pattern to ring-CC's non-antisymmetrized amplitude never right in the first place -- FLAGGED for Ajith, needs the actual RCCD spin-adapted equations, not another code-side substitution guess

**Attempted fix, per Ajith's direct call that the ~1.86-1.92x ratio
(not clean 2.0) IS the bug**: mirrored `vcc/dwmbej.F`'s RHF ABAB block
(lines 434-435, read-only reference, NOT modified) in
`rcc_dwmbej.F`'s RHF/IMULT=1 branch -- added `CALL GETLST(...,LISTTB)`
+ `SSCAL(...,TWO,...)` + `SAXPY(...,ONEM,...)` to build `2*T(list37) -
T(list34)` before the `CNTRCT` T*W ring contraction, exactly matching
CCSD's `2*T(list37)-T(list39)` pattern structurally (same SSCAL/SAXPY
idiom, same I010/I020 buffer reuse timing, verified safe against the
routine's own later I010 reuse for the bare-J integral read).

**Result: WORSE, not better.** Rebuilt `xringcc`
(`f6acd6e76e1a5a9cf30d72946f8a12e7`), ran H2O/STO-3G RHF fresh
(`/blue/ufhpc/perera/scrnc_test/h2o_scrnc_rccd_rhf_sto3g_dwmbejfix`).
Hbar eigenvalues: still 2 negative, and BOTH became MORE negative than
the unfixed baseline (`-0.88460752`->`-1.06894620`,
`-0.62703242`->`-0.68825671`). This is decisive: the fix moves the
result in the WRONG direction, ruling out "subtract list34 from
2*list37" as the correct term, not just leaving it unverified.

**Traced WHY, via `RCC_RNABIJ` (ring-CC's own analog of CCSD's
`RNABIJ`, called from `rcc_driver.F`) -- this settles list34's actual
identity definitively, not just by inference from `check_t2.F`'s
generic label as before.** For `TYPE='T'` (the T2-amplitude resort
mode): list37/39 (`LIST(1)`/`LIST(3)`) ARE built together, for RHF too,
via `SST002`+`SSTRNG` resorts of the SAME reference list (`REFLST=46`,
the primary `T(Ij,Ab)` amplitude) -- confirming the earlier session's
"`T(list37)==T(list39)` exactly" finding is a real mathematical fact
(RCCD's non-antisymmetrized `T(Ij,Ab)` is permutationally invariant
under the specific index-swap `SST002`/`SSTRNG` implement), not a code
bug to work around. List34 (`LIST(5)`), by contrast, comes from a
COMPLETELY SEPARATE second loop in the SAME routine, built via
`RCC_SST003` from **`REFLST-3+ISPIN` = 44, NOT 46**, with
`SPCAS2='AAAA'` explicitly labeled same-spin. List34 is genuinely,
structurally NOT a permuted view of the list37/39/46 ABAB-tensor
family at all -- it is built from an entirely different source list
(44). This matches (and now explains, rather than just observes)
`rcc_dwmbej.F`'s OWN IMULT=2 branch, which reads list34 as the
TRIPLET channel's PRIMARY amplitude (`CALL GETLST(...,LISTTB)` where
`LISTTB=34`) -- not as a correction term folded into the singlet
formula. Reusing the triplet channel's own primary amplitude inside
the singlet (IMULT=1) construction, the way this attempted fix and the
earlier `Scrnc_form_mubar_s.F` fix (2026-08-19, "continued 6" below)
both did, has no structural justification once you trace where list34
actually comes from -- it's not "the RCCD-specific stand-in for CCSD's
antisymmetrization partner," it's a different physical quantity
(same-spin T2) with no general reason to appear in a singlet-channel,
opposite-spin (ABAB-type) formula at all.

**Reverted the fix cleanly**: `git checkout -- aces2/ringcc/rcc_dwmbej.F`,
confirmed `git diff` empty against HEAD, rebuilt `xringcc`
(`f7d67117453febd3ae3d596a9cc89823`, matches the pre-session-start
binary). Note for future sessions, NOT a new issue: this rebuild
surfaced a PRE-EXISTING, unrelated inconsistency -- the git-committed
`aces2/ringcc/rcc_dwmbej.f` (the auto-generated preprocessed file,
also tracked in this repo) is stale relative to the committed `.F`
(missing a `MERGE(...)` checksum-size expression that IS present in
`.F` and gets correctly regenerated by any fresh build) -- cosmetic,
diagnostic-print-only, doesn't affect any computed result, not
introduced this session, not fixed here (out of scope, would need a
`.f`-regeneration-and-recommit housekeeping pass, unrelated to RCCD
physics).

**This also puts a question mark over the EARLIER, already-shipped
`Scrnc_form_mubar_s.F` fix** (2026-08-19 "continued 6" below, `LISTT2B
= 34` for RCCD, committed in `498542f8`) -- it rests on the SAME
premise (substitute list34 for the trivial list39) that today's
attempt just showed doesn't hold up once you trace list34's actual
construction. That fix DID produce a real, reproducible numerical
improvement at the time (H2/STO-3G: 2.45%->1.41%, H2/DZ: 11.6%->1.02%)
-- so it is NOT necessarily wrong in the same direction as today's
attempt (property-dressing and Hbar-construction are different
formulas, and "improves the number" is real, if inconclusive, evidence
that today's Hbar attempt didn't have). But given list34 is now known
to be the triplet channel's own primary amplitude rather than any kind
of singlet exchange partner, that earlier fix's own theoretical
grounding is now suspect too, and the ~1-1.4% residual it left behind
unresolved (never fully closed even at the time) may be related. NOT
re-investigated this session (would need to re-verify H2/STO-3G and
H2/DZ RCCD again with `Scrnc_form_mubar_s.F` reverted to its pre-fix
state to see whether removing it changes anything now, given today's
Nmult=1 architecture postdates that fix) -- flagged for whoever picks
this up next, not acted on, to avoid another unverified guess-and-test
cycle without Ajith's input first.

**FLAGGED for Ajith, sharper than before but still not code-resolvable
alone**: does RHF RCCD's Hbar ring-term (`rcc_dwmbej.F`'s IMULT=1
`T*W` contraction) need ANY exchange-type correction to bare `T(37)`
at all? The routine's own header comment ("There are no antisymmetrized
T2 amplitudes") may have been right all along for THIS specific
formula -- CCSD's `2*T37-T39` pattern exists because CCSD's `T2` is
genuinely antisymmetrized (so its own list39 is a real, distinct
quantity); RCCD's `T2` is NOT antisymmetrized by construction, so the
straightforward analogy may simply not transfer, and the ~1.86-1.92x
deviation from UHF's per-spin value may come from something else
entirely (a missing FACTOR/prefactor somewhere in the T*W contraction
call itself -- e.g. the `X`/`C4` scale arguments passed to `CNTRCT`,
worth comparing literally line-by-line against CCSD's own call
signature at vcc/dwmbej.F line 446-448, not yet done -- or from the raw
`2J-K` bare-integral piece, or from somewhere entirely outside
`rcc_dwmbej.F`). This needs either the actual closed-form RHF-RCCD
spin-adapted Hbar equation (ideally cross-checked against the
Rishi/Perera/Bartlett JCP153,234101(2020) paper already used as this
thread's precedent for the energy formula) or Ajith's own recollection
of the correct construction -- not another blind list-substitution
attempt. Do not try list34, list39, or any other existing list as a
"missing term" again without first confirming from the derivation
which quantity (if any) is actually supposed to appear there.

**Housekeeping**: no fix committed (the attempted one was verified
WRONG and reverted). Working tree: `rcc_dwmbej.F` back to exact HEAD
state (confirmed via empty `git diff`); `rcc_dwmbej.f` now shows a
small, pre-existing, unrelated diff (see above, not this session's
doing, harmless); `scrnc_form_rcc_hbar.F`'s pre-Fock debug print from
the previous entry, and `scrnc_set_vars.F`'s RHF-RCCD-trap test
override, both still in place as before, still not shipped. Test
artifact (not committed): `h2o_scrnc_rccd_rhf_sto3g_dwmbejfix`.

---

## SESSION 2026-08-20 (continued 2) -- traced the 2 negative Hbar eigenvalues down to a specific term: definitively RULED OUT "missing Fock-diagonal" (it's added correctly, verified against SCF orbital gaps to 8 digits); narrowed to the raw 2J-K correlation-kernel diagonal for the HOMO(MO5)->virt(MO6,MO7) pairs being large enough (~-1.0 to -1.15 Ha) to nearly cancel the small HOMO-LUMO-type gap (~0.99-1.13 Ha) and flip the net sign -- whether that raw kernel magnitude is itself correct is a physics judgment call, FLAGGED for Ajith, not resolved by code-reading alone

**Added ONE minimal, targeted debug print** (per Ajith's "lean on
existing infra, not a lot of new one-off prints" -- the existing
instrumentation didn't cover this specific pre-Fock intermediate, so
added exactly one, matching the existing `_DEBUG_LVL0`/`call output`
style): in `scrnc_form_rcc_hbar.F`, right after the `IMULT.EQ.1` raw
`GETALL(...,56)` read and before the `(ai,bj)->(ab,ij)` reorder /
`SCRNC_ADDONEH` Fock-diagonal-add call, dumps the pristine list-56
matrix in the SAME `(ai,bj)` index convention as the final "drCCD Hbar
matrix" dump -- so row/col N here is directly comparable to row/col N
in the already-known post-Fock dump, no reindexing needed. RHF-only,
`IMULT.EQ.1` only, `_DEBUG_LVL0`-gated (kept in the working tree
uncommitted, matching this week's `scrnc_set_vars.F` test-override
precedent -- needed for continued testing, not yet shipped).

**"Missing Fock-diagonal" hypothesis: definitively RULED OUT.** Ran
H2O/STO-3G RHF fresh with the new print
(`/blue/ufhpc/perera/scrnc_test/h2o_scrnc_rccd_rhf_sto3g_prefock_test`).
Diffed the new pre-Fock dump against the already-known post-Fock "drCCD
Hbar matrix" dump, row by row, and checked the difference against the
actual SCF orbital-energy gaps (`ORBITAL EIGENVALUES` block, same run):
```
Row  MO pair     Fock added (post-pre)   Computed gap eps(virt)-eps(occ)
1    MO1->MO6    20.84501690838           20.845 (exact match, deep core)
3    MO2->MO6     1.86991641999            1.870 (exact match)
9    MO5->MO6     0.99399238313            0.99399238313 (EXACT, 11 digits)
10   MO5->MO7     1.13009122331            1.13009122331 (EXACT, 11 digits)
```
The Fock-diagonal term is being added correctly and completely for
EVERY row checked, including rows 9/10 (the anomalous ones) -- to the
full 11 printed digits, no truncation, no sign error, no missing term.
This rules out the entire "SCRNC_ADDONEH"/list-routing/irrep-blocking
class of hypothesis that motivated adding the print in the first
place.

**Narrowed the anomaly to the raw correlation kernel itself.** Row
9/10's PRE-Fock diagonal (the raw `2J-K`-type kernel from
`rcc_dwmbej.F`, list 56) is `-1.02650109829` / `-1.15313423144` --
large enough in magnitude to nearly cancel the (comparatively small)
`0.994`/`1.130` Ha HOMO(MO5,B1)->LUMO(MO6/MO7) gap, flipping the net
sign. This ISN'T obviously anomalous in isolation -- EVERY row's raw
kernel is O(0.05-1.65) Ha in magnitude (row1: -1.48, row2: -1.65, row8:
-0.95, etc.) -- rows 9/10 just happen to be the ones where the gap
itself is small enough for the kernel to dominate and flip the sign,
rather than being a small correction on top of a large gap (rows 1-4,
deep-ish orbitals) or staying subdominant anyway (rows 5-8).

**Cross-checked against UHF's analogous per-spin value**: UHF's alpha-
only post-Fock diagonal at the same composite indices (9,10 of its
20x20 matrix, confirmed same `(i-1)*NVRT+a` convention, alpha block =
rows/cols 1-10) is `0.44207134278` / `0.52913151138` -- POSITIVE,
unlike RHF. Backing out UHF's own Fock gap (same molecule, same
orbitals, must match RHF's `0.994`/`1.130` exactly) gives UHF's raw
per-spin correlation kernel as `-0.552`/`-0.601`. Ratio RHF-singlet-
kernel : UHF-alpha-kernel = **1.86 and 1.92** -- close to but not
exactly 2.0. This is right on the boundary of "looks like the expected
`2J-K` (RHF, roughly double a per-spin `J-K`-type UHF quantity)
qualitative relationship" vs. "could still be hiding a real, smaller-
scale error" -- NOT distinguishable from code-reading or checksums
alone; would need either (a) the actual closed-form CC-theory relation
between RHF's spin-adapted singlet ring-CC Hbar diagonal and UHF's
spin-orbital one (a derivation, not a code trace), or (b) Ajith's own
read on whether that ~1.86-1.92x ratio (vs an exact 2.0) is within
expected non-trivial-spin-adaptation territory or is itself the
smoking gun.

**FLAGGED for Ajith, not guessed at**: is the RHF ring-CC singlet
Hbar's raw correlation-kernel diagonal, for an occ-virt pair where the
kernel and the bare orbital gap are comparable in magnitude, supposed
to relate to UHF's per-spin analogue by a clean factor of 2 (in which
case today's ~1.86-1.92x is the bug, and worth chasing directly in
`rcc_dwmbej.F`'s IMULT=1 branch -- candidates there:
the `2J-K`-vs-`-K` bare-integral construction identified in "continued
6" above, or the `T*W` contraction step whose checksums are
already-existing instrumentation, not yet compared row-by-row) -- or
is some deviation from exactly 2x physically expected once genuine
spin-adaptation (not just a naive per-spin doubling) is done correctly,
in which case the bug is elsewhere entirely (e.g. genuinely in how
`Scrnc_form_mubar_s.F`/`scrnc_jj_rcc_assemble.F` CONSUME this Hbar
downstream, not in the Hbar itself) and today's whole Hbar-eigenvalue
lead would be a red herring dead-end for RCCD. Do not proceed further
on the Hbar-construction line without this call -- it determines
whether the next step is "trace `rcc_dwmbej.F` term-by-term" or "go
back to comparing consumption of a CORRECT Hbar."

**Housekeeping**: no fix committed (none verified). Working tree now
also has the new debug print in `scrnc_form_rcc_hbar.F`/`.f`
(uncommitted, `_DEBUG_LVL0`-gated, RHF/IMULT=1-only, kept for
continuity -- harmless, adds output only, changes no computed values).
New test artifact (not committed):
`/blue/ufhpc/perera/scrnc_test/h2o_scrnc_rccd_rhf_sto3g_prefock_test/`
(first attempt in this dir crashed with `Density=Response is not set`
-- reused the wrong ZMAT template missing the keyword; refixed by
copying `_cleanbaseline`'s ZMAT, which has it; the failed first attempt
is not a real finding, just a self-inflicted test-setup mistake, noted
here so it isn't mistaken for something else later).

---

## SESSION 2026-08-20 (continued) -- top-down debug-print bisection (Ajith's method: use EXISTING instrumentation, start high, walk backward, stop at first divergence): FOUND IT -- RHF's ring-CC Hbar has 2 spurious NEGATIVE eigenvalues (-0.88, -0.63) that UHF's equivalent spectrum (all 20 positive) has no trace of; this is upstream of and likely explains everything chased in every prior RCCD session

**Method, per Ajith's explicit instruction**: don't write new debug
prints -- the codebase already has extensive `_DEBUG_LVL0`-gated
instrumentation throughout the RCCD path (confirmed: this macro is
ALREADY on by default in this build,
`aces2/Makefiles/GNUmakefile.src` line 10 uncommented `DEFINES_EXTRA
= ... -D_DEBUG_LVL0 ...` -- nothing needed to "turn on", it was already
firing in every run this week, including today's `_cleanbaseline` RHF/
UHF logs from earlier this session). Walked the existing H2O/STO-3G
RHF and UHF `_cleanbaseline` run logs top-down: final JJ matrix (known,
9.58%/23.7x gap) -> response-only JJ (before bare-v add) -> per-AO-pair
R-side/L-side checksums -> Hbar eigenvalues. Stopped as soon as
something looked wrong, per instruction, rather than continuing deeper.

**Housekeeping note first**: the UHF `_cleanbaseline` log
(`/blue/ufhpc/perera/scrnc_test/h2o_scrnc_rccd_uhf_sto3g_
cleanbaseline/run_scrnc.log`) contains 7 stray embedded NUL bytes
(harmless -- looks like a raw/unformatted debug print of an integer
buffer somewhere upstream, likely in `xintprc`'s `@INTRN1` integral-
transform debug block, not scrnc-specific) that make plain `grep`
silently treat the file as binary and stop matching partway through --
this is why last entry's "UHF log has no RPA JJ matrix" note was a
FALSE ALARM, not a real bug: the full pipeline ran correctly and DOES
print the JJ matrix (`0.86509751290`/`0.03016663846`, matches the
REGROUP session exactly), just invisible to naive `grep -n`. Retract
that flagged concern. Used Python (`.decode('latin-1')`, byte-safe) for
all log analysis below instead of `grep`.

**Level 1 -- Response-only JJ matrix (before bare-v add),
`scrnc_jj_rcc_assemble.F`'s own `_DEBUG_LVL0` dump**: RHF JJ(4,4) =
`0.69148053436`, UHF JJ(4,4) = `0.00694352740` -- ratio **~99.6x**,
even bigger than the final 23.7x. Checked the bare-v piece separately
(full JJ(4,4) minus response-only JJ(4,4)): RHF `0.02322`, UHF
`0.02322` -- IDENTICAL. So the bare (uncorrelated) Coulomb piece is
completely fine and reference-independent as expected; the entire
23.7x final discrepancy is 100% inside the response construction,
and is actually LARGER (99.6x) before the shared bare-v addition
partially dilutes it. Confirms: don't waste time on `scrnc_prep_2ints.F`
or the bare-v list-220 machinery, it's clean.

**Level 2 -- per-AO-pair R-side checksums, `scrnc_lrspn_rcc_main.F`**:
extracted the `Npert_count=4` block (AO-pair index 4, same index as the
bad JJ(4,4) component) for both refs. The DRESSED bare property
checksum (`Mubar_s(Iside=1)`, computed BEFORE the Hbar^(-1) contraction)
is RHF `0.117846531745E+00 / 0.707396779915E-02` vs UHF
`0.238630500016E+00 / 0.142786183332E-01` -- a clean, consistent ~2.02x
ratio on both checksum components. Not obviously wrong (plausibly just
an ISPIN-loop bookkeeping artifact, UHF's routine summing over 2 spins
where RHF's sums over 1) -- didn't chase further given what's below.
Immediately AFTER the Hbar^(-1) contraction (the `1st order amps`
checksum, i.e. `R = Hbar^(-1) * Mbar`): RHF
`-0.375279167468E+00 / 0.347799818980E+00` vs UHF
`0.665798062707E-01 / 0.354861392438E-02` -- NOT a clean ratio anymore,
DIFFERENT SIGN on top of a ~5-10x magnitude mismatch. This pins the
first clearly-anomalous step to the Hbar^(-1) contraction itself (or
whatever it's built from) -- not the bare-property dressing, which
still looks internally consistent going in.

**Level 3 -- Hbar eigenvalues, `scrnc_form_rcc_hbar_inv.F`'s own
`_DEBUG_LVL0` `DGEEV` diagnostic (Ajith's own 2026-08-14 request,
"independent Hbar diagonalization... compare against ringcc's own
eigenvalues") -- FOUND IT.** (Note: text-matching "The real
eigenvalues" alone is ambiguous, it also matches an unrelated SCF-
orbital-energy print elsewhere in the pipeline -- anchor on the unique
"scrnc's own Hbar(ai,bj) is fully diagonalized" header instead, which
prints exactly once per reference, immediately followed by the REAL
Hbar spectrum.)
```
RHF Hbar eigenvalues (MATDIM=10, singlet-channel, list 56):
  526.96238440   3.84383016  19.76892375  29.26769779  525.99926507
   27.48199452   5.49605640  10.63082475  -0.88460752   -0.62703242

UHF Hbar eigenvalues (MATDIM=20, alpha+beta combined):
  547.09011687 545.39016197 548.45785580 547.29146058  39.68624555
   34.63841383  28.39704398  40.98300845  37.87447405  19.73887279
   16.60471345  12.80987637  21.88606748  19.03727989  17.88331701
   10.97710036  14.60126343  13.08165647  13.73465055  15.06216374
```
**RHF has 2 NEGATIVE eigenvalues (`-0.88460752`, `-0.62703242`). UHF's
full 20-eigenvalue spectrum has ZERO negative or near-zero values --
smallest is `10.97710036`.** This is not a subtle numerical-precision
issue (RHF's own two negative values are ~0.6-0.9 in magnitude, nowhere
near the ~1e-8 convergence floor this codebase shows everywhere else)
-- it's a qualitative, structural difference. `rcc_dwmbej.F`'s RHF
IMULT=1 construction is the standard closed-shell Casida singlet
`A=2J-K` matrix (confirmed in an earlier session, "continued 6" above)
-- for a genuine, stable, non-degenerate closed-shell ground state
(H2O at equilibrium, no reason to expect an RPA/electronic instability)
this matrix should be positive-definite, like an orbital Hessian.
UHF's matching spectrum being cleanly all-positive is exactly the
physically-expected behavior; RHF's two negative values are the
anomaly. High-energy eigenvalues correspond loosely across references
(RHF's `526-527` pair vs UHF's four `545-548` values -- roughly the
expected 1:2 ratio if UHF's alpha/beta are near-degenerate copies of
what RHF represents once), so the CONSTRUCTION isn't wholesale broken
-- something specific and localized (a missing term, sign error, or
incomplete contraction affecting only the lowest 1-2 excitation(s) in
the singlet manifold) is corrupting just the low end of RHF's spectrum.

**Per Ajith's bisection instruction ("stop descending as soon as you
find the first stage where something looks wrong") -- STOPPING HERE.**
This is upstream of everything chased in every prior RCCD session
(mubar/T2-dressing, Imult/Loffset routing, singlet/triplet
recombination -- all of those operate on OR consume this same Hbar,
so a corrupted Hbar spectrum would explain the response-level symptoms
at every level checked above and in every previous entry in this file,
without needing a second, independent bug anywhere else). NOT yet
traced into `rcc_dwmbej.F` itself to find which specific term produces
the 2 bad eigenvalues -- that's the concrete next step, not done this
session (this entry is the "first divergence found" checkpoint per
Ajith's own stated method, not the final root-cause trace).

**Next step, concretely**: open `ringcc/rcc_dwmbej.F`'s RHF (`IUHF.EQ.0`,
`IMULT.EQ.1`) branch and either (a) print the raw 10x10 Hbar matrix
itself (already available via the existing `_DEBUG_LVL0` "The drCCD
Hbar matrix" `call output(...)` dump right before the DGEEV call in
`scrnc_form_rcc_hbar_inv.F`, not yet inspected element-by-element this
session) and look for an obviously wrong row/column/asymmetry, or (b)
compare it term-by-term against `vcc/dwmbej.F`'s RHF construction
(CCSD's, known-correct, never decouples) the same way the "continued 6"
entry above already did for the OVERALL `2J-K`/`-K` structure -- but
this time at the level of individual contraction terms, not just the
top-level formula, to find which specific piece is missing or sign-
flipped for whichever occ-virt pair(s) are producing the negative
eigenvalues. A finite-difference or symmetry cross-check (e.g. is the
matrix supposed to be symmetric? `rcc_dwmbej.F` builds a non-symmetric
Hbar in general, so don't assume symmetry should hold, but DO check
whether the pristine "drCCD Hbar matrix" dump itself already looks
obviously wrong in a specific row/column before assuming the bug is in
the diagonalization/eigenvalue step rather than the matrix construction
itself).

**Housekeeping**: no code changes this entry (pure log-based
investigation of the two already-completed `_cleanbaseline` runs, no
new binary rebuilds or submissions needed since `_DEBUG_LVL0` output
was already captured). Working tree unchanged from the entry below.

---

## SESSION 2026-08-20 -- RCCD-only session (per Ajith's mandate): rediscovered that the "singlet-only" response hypothesis is ALREADY committed (498542f8) and does NOT close H2O's gap; proved `DENSITY=RESPONSE` is bit-identical-irrelevant to RCCD's own numbers (Ajith's correction, confirmed empirically); real bug (H2O/STO-3G 9.58%/23.7x gap) remains OPEN, no new fix

**Two live corrections from Ajith mid-session, both now reflected in
how this thread should run going forward:**
1. `xvprops` must NEVER be in scrnc's own manual pipeline (dipole/
   property integrals scrnc's AO-pair-density construction never
   consumes) -- only the `PROPS=EOM_NLO`/`xvcceh` reference pipeline
   needs it.
2. **Supersedes item 1's `DENSITY=RESPONSE` guidance and the "ruled
   out" conclusion in the 2026-08-19 entry below** ("`DENSITY=RESPONSE`
   hypothesis tested and RULED OUT"): Ajith says do NOT use
   `DENSITY=RESPONSE` in EITHER pipeline at all -- it only matters for
   whether `xdens` runs, and neither pipeline (scrnc's own, or the
   `xvcceh` reference) ever calls `xdens`. The earlier finding that
   removing it crashes `xlambda` on list 144 was about the CCSD/
   `xvcceh`-comparison thread specifically (out of scope today) --
   not re-verified this session, flagged for whoever picks that thread
   back up. **For RCCD specifically, verified directly this session**
   (see below): the keyword is provably inert.

**RCCD's own pipeline does NOT need `DENSITY=RESPONSE` -- proven, not
assumed.** Built a fresh H2O/STO-3G RHF RCCD test
(`/blue/ufhpc/perera/scrnc_test/h2o_scrnc_rccd_rhf_sto3g_nodensitytest`,
manual pipeline `xjoda,xvmol,xvmol2ja,xvscf,xvtran,xintprc,xringcc,
xringcl,xscrnc`, NO `xvprops`, NO `DENSITY=RESPONSE` in the ZMAT) against
an otherwise-identical clean baseline WITH the keyword
(`..._cleanbaseline`) -- required temporarily bypassing
`scrnc_set_vars.F`'s own hard gate on `Respn_d` to even test this
(`.FALSE. .AND. .NOT. Respn_d`, test-only, reverted before this entry
was written, NOT shipped). `diff` of the two full run logs: the ONLY
differences are ASLR-randomized `@ACES_INIT-DEBUG` memory addresses and
the `DENSITY` flag echo itself -- every actual numeric result, including
the full RPA JJ matrix to all 11 printed digits, is bit-identical.
RCCD's pipeline never calls `xdens` (confirmed: only consumer of the
`Respn_d`/`IFLAGS(h_IOPPAR_density)` flag we could find that's actually
in RCCD's call path is `scrnc_set_vars.F`'s own gate -- traced
`Scrnc_set_dens_vars.F`, the file that actually USES this flag
meaningfully via `COMMON/DERIV/`, and confirmed it's `xdens`-specific
machinery RCCD's manual pipeline never invokes) and never calls
`xlambda` either (ring-CC builds its own Hbar entirely in-module). So
for RCCD, `scrnc_set_vars.F`'s `Respn_d` gate is pure legacy
gatekeeping with zero effect on the actual computation -- worth
permanently removing, but NOT done this session (the gate is shared
with the CCSD path in the same file, which is out of scope today and
was NOT re-verified against this same claim -- don't blanket-delete the
gate without checking CCSD's own dependency on it first, in case that
thread's `xlambda`/list-144 crash is genuinely tied to this flag rather
than to `PROPS`/`IHBAR2` as Ajith suspects).

**Rediscovered, the hard way, that the "singlet-only" response
hypothesis is not new -- it's already shipped, in `498542f8`, and
already known to be insufficient.** Independently derived (before
looking closely at `scrnc_driver.F`) that the density-density response
function feeding `W=v+v*Chi*v` should, by the Wigner-Eckart theorem, get
ZERO contribution from the triplet (Imult=2) channel -- a spin-scalar
one-particle density operator has no matrix element between a singlet
ground state and a triplet excited state, so `scrnc_jj_rcc_assemble.F`'s
0.25(singlet)/0.75(triplet) combination looked like a plausible bug
(reusing the ground-state-ENERGY recombination formula somewhere it
doesn't belong). Wrote a test edit (`Wt=1.0` for Imult=1, `Wt=0.0` for
Imult=2, `rcc_cmpeng.F`'s own energy formula untouched) and reran H2O/
STO-3G RHF RCCD fresh -- result was **bit-identical to the untouched
baseline** (JJ(1,1)=`0.94799229627` both ways, JJ(4,4) differed only in
the 10th digit, clearly just Davidson-iteration noise). Root cause of
the null result: `scrnc_driver.F` ALREADY hardcodes `Nmult = 1` (see its
own header comment, dated 2026-08-19, committed in `498542f8`) --
"confirmed (Ajith) the singlet/triplet channels must NOT be decoupled
and weight-recombined for a response quantity... A first-principles
check (density-type perturbations couple only to the singlet manifold,
zero overlap with the Sz=0 triplet excitation operator) supports Nmult=1
(singlet-only)" -- i.e. this EXACT hypothesis, with nearly identical
reasoning, was already tested and shipped a session ago. Because Nmult
is now always 1, `scrnc_jj_rcc_assemble.F`'s `If (Nmult .EQ. 2)` branch
(where the 0.25/0.75 weights live) is dead code -- my test edit there
could never have had any effect, hence the null result. Reverted the
test edit (back to matching HEAD exactly, verified via `git diff`).
**Important for future sessions: don't re-derive or re-test the
singlet-only-channel hypothesis again -- it's already the shipped
architecture, and `498542f8`'s own commit message says the result was
"gap 3.98%->~1%" on H2/STO-3G, with an explicitly acknowledged
"residual normalization-scale mismatch against UHF... not yet
resolved."** The REGROUP session's H2O/STO-3G 9.58%/23.7x gap (below)
was measured AFTER this Nmult=1 fix was already in place -- it is that
same unresolved residual, just far more exposed on H2O than it was on
H2, not a new or different bug.

**The real open bug, unchanged, not newly investigated this session
beyond the above**: H2O/STO-3G RHF RCCD's RPA JJ matrix vs UHF --
JJ(1,1) 9.58% off, JJ(4,4) 23.7x off (`0.71470364552` vs
`0.03016663846`), several off-diagonal elements (JJ(1,2),JJ(1,3),
JJ(2,2)) also off by several percent, NOT just component 4. Sharpest
lead, not yet chased: component 4 is EXACTLY symmetry-decoupled from
components 1-3 in BOTH references (RHF off-diagonals are exact zero to
11 digits; UHF's are ~1e-7, consistent with Davidson-convergence noise
around a true zero) -- so index-4's AO pair belongs to a symmetry
sector that doesn't mix with 1-3's sector in either reference, yet its
OWN diagonal magnitude is wrong by 23.7x while the mixed-sector
diagonals (1,1)/(2,2)/(3,3) are only off by single-digit percent. A
uniform normalization/scale bug would hit every component by the same
factor -- this pattern (isolated component wildly wrong, mixed
components mildly wrong) instead suggests something specific to
whichever symmetry sector AO-pair-4 belongs to, possibly a missing or
mis-scaled contribution that gets partially diluted/cancelled in the
mixed components but shows up at full strength in the pure one. Have
NOT yet identified what that sector is (which MOs/irrep AO-pair-index 4
maps to) or traced a specific routine -- this is the concrete next
step, not the singlet/triplet-channel question (already closed, see
above) and not the `DENSITY=RESPONSE`/`xvprops` pipeline-composition
question (already closed, see above).

**Housekeeping**: nothing new committed this session (pure
investigation + two reverted test edits, verified via `git diff`
matching session-start state exactly). Working tree still has the
same pre-existing uncommitted state as every session this week:
`scrnc_set_vars.F`'s RHF-RCCD-trap `.FALSE. .AND.` override (still
needed for testing, still must NOT ship) and the whitespace-only
`scrnc_form_rcc_hbar.F` diff. New test artifacts (not committed, data
only) under `/blue/ufhpc/perera/scrnc_test/`:
`h2o_scrnc_rccd_{rhf,uhf}_sto3g_cleanbaseline` (RHF/UHF, manual
pipeline, no `xvprops`, WITH `DENSITY=RESPONSE`, matches the REGROUP
session's numbers exactly), `h2o_scrnc_rccd_rhf_sto3g_nodensitytest`
(same, `DENSITY=RESPONSE` removed, bit-identical to `_cleanbaseline`,
proves the keyword inert for RCCD), `h2o_scrnc_rccd_rhf_sto3g_
singletonly_clean` (the now-understood-to-be-inert Wt=1/0 test, kept
for reference). Note: `_cleanbaseline`'s UHF run's log shows no "RPA JJ
matrix" print at all despite exit code 0 -- not chased down this
session (lower priority than the main open bug, and the RHF-side
numbers needed for today's tests were already available and consistent
with the REGROUP session's own numbers), worth a quick look next
session in case it's a real (separate, UHF-specific) issue rather than
a `_DEBUG_LVL0` build-flag artifact.

**Next step, concretely**: identify AO-pair-index-4's symmetry/MO
origin (probably via the `_DEBUG_LVL0` `Isympert`/`Iaopop` diagnostic
prints already built into `scrnc_jj_rcc_assemble.F`, or a rebuild with
that ifdef enabled), then trace which specific routine builds that
sector's response contribution differently between RHF and UHF. Given
today's finding that this is NOT a singlet/triplet-combination issue
(that's already fixed/closed), the remaining candidates are back to
`rcc_dwmbej.F`'s own list-56 Hbar construction, `Scrnc_form_mubar_s.F`'s
list-34 T2-dressing fix (verified on H2/STO-3G's single-irrep case
only -- never re-verified on a genuinely multi-irrep system like H2O),
or `rcc_cmpeng.F`-adjacent list-routing now that Nmult=1 changed which
lists are actually populated vs `scrnc_jj_rcc_assemble.F`'s Loffset=0
assumption.

---

## SESSION 2026-08-19 (continued 7) -- "REGROUP" session: proved the CI-like framework itself is NOT the problem (RHF=UHF to ~1e-8 via the correct ACES2 pathway); the bug is entirely within scrnc's own construction (or whatever DENSITY=RESPONSE triggers, even for plain CCSD). RESUME POINT for 2026-08-20: check scrnc's CCSD path first, then ringcc/RCCD (RESUME HERE)

**Context**: after the H2/STO-3G RCCD list-34 fix (see "continued 6"
below) got the gap to ~1-1.4%, Ajith asked to regroup: stop chasing H2
(too trivial, 1-orbital-pair, possibly masking/distorting the real
picture) and re-derive everything on H2O/STO-3G instead, proving CCSD
RHF=UHF works FIRST before trusting anything about RCCD.

**H2O/STO-3G RCCD, fresh, via scrnc (`xjoda...xringcc,xringcl,xscrnc`,
`DENSITY=RESPONSE`)**: RHF vs UHF RPA JJ(1,1) gap **9.58%**, and
`JJ(4,4)` differs by ~23.7x (`0.71470364552` vs `0.03016663846`) --
MUCH worse than H2 ever showed, while the ground-state energy still
matches to ~1e-10 (`-74.994555231233` vs `-74.994555231114`). Confirms
H2's ~1-1.4% residual was NOT representative -- H2O's richer symmetry
structure exposes a much bigger problem in the response construction
specifically.

**Sanity-checked the H2O/STO-3G CCSD baseline via scrnc's own
`DENSITY=RESPONSE` pipeline first**: only ~0.0145% RHF/UHF agreement
(`0.86106633095` vs `0.86094156406`) -- suspiciously loose given CCSD-LR
was verified to ~1e-8 on H2 earlier this same session. Ajith flagged
this directly: "should match 10^-6 or better."

**Verified independently with ACES2's own, decades-old, dedicated
polarizability code (`xvcceh`, J.F. Stanton 1993 -- nothing scrnc-
related)**, run in a genuinely fresh directory (a first, contaminated
attempt where `xscrnc` had already run in the SAME directory gave
wildly wrong numbers -- classic non-idempotent-list-reuse trap, redo
in a clean dir): `DENSITY=RESPONSE`'s own polarizability came out
~0.2-0.5% RHF/UHF, NOT the expected precision, and threw an
intermittent (non-deterministic, doesn't correlate with the actual
numeric result -- confirmed by rerunning and getting the identical
number with no error message) `Intel oneMKL ERROR: Parameter 13 was
incorrect on entry to DGEMM` in the RHF run only. The DGEMM message
looks like a red herring (uninitialized-memory read warning), not the
real cause.

**Found the actual correct methodology, at Ajith's direction**: the
real ACES2 test suite (`aces2/test/script.075b`/`075c`, H2O/DZP,
`OH=0.957201, ANG=104.5200`) uses `PROPS=EOM_NLO,EOMPROP=QUADRATIC`
via the standard `xaces2` driver -- NOT `DENSITY=RESPONSE`. Ran it: RHF
and UHF "Molecular Polarizabilities" match to all 6 printed digits
(`0.048619, 4.915536, 2.108859`), energy to ~1.8e-9. But Ajith
clarified QUADRATIC doesn't apply to scrnc -- "we have not added the
quadratic terms, we are in the CI-LIKE framework" (confirmed by an
existing prior-session pilot script,
`aces2/scrnc/independent_pilot/ci_like_vs_lr_polarizability.py`,
which independently verified via PySCF that CI-like [`IEOMPROP=0`,
frozen-Lambda expectation value] and LR/QUADRATIC [`IEOMPROP=2`,
fully-relaxed via Hellmann-Feynman/Schwarz] are genuinely different
constructions with different guaranteed-symmetry properties -- that
pilot was about xy/yx Onsager reciprocity specifically, not RHF/UHF,
so it didn't directly answer this question, but it's the right
precedent/reference for how "CI-like" is defined in this codebase).

**Reran the SAME test-suite input with `EOMPROP` simply omitted**
(defaults to `IEOMPROP=0`, CI-like) -- this is the DECISIVE result:
```
CI-like Pol. matrix
RHF: 2.92266723692  7.56686877462  5.89953761093
UHF: 2.92266716621  7.56686860801  5.89953760962
```
Agreement ~2e-8 across all three components -- **CI-like itself has NO
bug and gives RHF=UHF to the same precision as QUADRATIC.** This
proves the ~0.5% mismatch seen earlier was never a CI-like theory
limitation -- it's specific to whatever code path `DENSITY=RESPONSE`
actually exercises (even for plain CCSD, independent of scrnc/RCCD
entirely), or to scrnc's own construction, or both.

**Conclusion, and the plan for tomorrow (2026-08-20)**: "everything is
in scrnc" (Ajith). Since `PROPS=EOM_NLO`'s CI-like route (in `vcceh`,
via `xvcceh`/the standard `xaces2` CCSD pipeline) is now PROVEN correct
to ~1e-8, and CCSD's OWN `DENSITY=RESPONSE`-triggered path is NOT
(only ~0.0145%-0.5%), the highest-value next step is a **direct,
side-by-side comparison of the two CCSD code paths** -- both are
notionally "CI-like," both live in the shared `vcceh` module family
(scrnc's own `Scrnc_form_mubar_00.F`/`Scrnc_form_mubar_s.F` are
scrnc-local near-duplicates of `vcceh`'s own routines, per earlier
session's file-comparison work) -- so there is a real, findable
divergence between them worth tracing to its exact origin. **Order of
operations, per Ajith**: (1) scrnc's CCSD path first (compare directly
against the proven-correct `PROPS=EOM_NLO` CI-like construction,
same H2O system, find exactly where/why `DENSITY=RESPONSE`'s numbers
diverge), THEN (2) once CCSD is understood and fixed, move to
ringcc/RCCD with that understanding in hand -- RCCD's much larger
H2O gap (9.58%, `JJ(4,4)` off by 23.7x) is very likely the SAME root
cause, just more exposed by ring-CC's structure.

**Added hypothesis (Ajith, right after this entry was first written):**
"I do not think density flag is needed. We don't run the density
code at all." -- i.e. `DENSITY=RESPONSE` may not actually be necessary
for scrnc's own construction at all; it may just be an assumed/legacy
requirement in `scrnc_set_vars.F`'s hard gate (see "Practical test-
running notes" earlier in this file: "scrnc_set_vars.F hard-requires
REF=UHF... and DENSITY=RESPONSE"). `DENSITY=RESPONSE` likely also
triggers a SEPARATE relaxed-one-particle-density computation (for
things like gradients/population analysis) that scrnc's own
`Scrnc_form_mubar_*` routines never actually consume -- if that extra,
unneeded machinery is where the imprecision/DGEMM issue lives, the
fix could be as simple as confirming scrnc only needs Lambda solved
(xlambda run, T/L amplitudes converged) and NOT the `DENSITY=RESPONSE`
keyword itself. **Cheap first check for tomorrow, before the deeper
line-by-line trace below**: does scrnc (or a plain CI-like `xvcceh`-
equivalent check) still work correctly with Lambda solved but
`DENSITY=RESPONSE` OMITTED from the ZMAT? If yes, that may directly
explain and fix everything chased this session with far less tracing.

**Checked immediately (before stopping): `xdens` hypothesis ruled out**
-- confirmed the WORKING `PROPS=EOM_NLO` pipeline (`xjoda,xvmol,
xvmol2ja,xvscf,xvtran,xintprc,xvcc,xlambda,xvprops,xvcceh`) never
calls `xdens` either (matches Ajith's "we solve lambda but we do not
call dens"), so `xdens` was never the differentiator. Also tried
moving `xvprops` to match the working pipeline's LATE position
(after `xlambda`, right before the property step) instead of scrnc's
usual early position (right after `xvmol2ja`) -- **no change at all**,
identical numbers (`0.04896437544`/`0.04871500302`) either way. Ruled
out.

**`DENSITY=RESPONSE` hypothesis tested and RULED OUT** (Ajith: "Just
eliminate density. Even that is harmless. They do nothing in scrnc" --
meant as "not the source of the bug," not literally removable, per
what actually happened): removing it from the ZMAT does NOT reveal
anything about the precision issue -- it just breaks the pipeline
outright. `xlambda` itself crashes first (`@GETLST: Assertion failed.
List (1,144) does not exist.` -- list 144 apparently only gets
registered when `DENSITY=RESPONSE` is set, needed by plain CCSD
Lambda regardless of scrnc), and `scrnc_set_vars.F`'s own explicit
gate refuses to proceed ("The Density=Response is not set in the ZMAT
and calculation can not proceed.") even if `xlambda` didn't crash
first. So `DENSITY=RESPONSE` is structural plumbing, required for the
pipeline to run at all -- it doesn't itself introduce the imprecision,
but it can't be tested-by-removal either. Confirmed harmless in the
sense Ajith meant (not the bug's source), ruled out as an active line
of investigation.

**Ajith's follow-up, key clarification, not yet acted on**: "xvprops
is for dipole and other integrals. We do not need that for scrncc" --
Ajith also confirmed `xvprops` itself is harmless/unnecessary
(separately from `DENSITY=RESPONSE`).
i.e. `xvprops` builds DIPOLE-type property integrals specifically --
scrnc doesn't use dipole integrals at all (it responds to AO-PAIR
DENSITY perturbations, a completely different, larger operator set) --
so `xvprops`'s presence/absence/position in scrnc's pipeline may be
entirely irrelevant busywork, inherited by convention from CCSD-LR/
dipole-property test setups rather than being something scrnc
actually needs. **Not yet tested**: does scrnc run correctly (and does
precision change) with `xvprops` REMOVED from the pipeline entirely?
This is the natural next cheap check, alongside re-examining whether
`DENSITY=RESPONSE` itself is even needed (Ajith's still-open point
from the previous entry) -- both are "maybe scrnc's assumed pipeline
carries unnecessary/wrong legacy steps" hypotheses, worth testing
together first thing next session, before the deeper line-by-line
`Scrnc_form_mubar_*` vs `vcceh` trace.

**ROOT CAUSE FOUND for the CCSD/CI-like path.** Traced why `xlambda`
crashes on list 144 without `DENSITY=RESPONSE`, and why even WITH it
scrnc's own numbers were only ~0.0145% precise (not ~1e-8): the
`IHBAR,IHBAR2,IHBAR2` control-variable line printed by `vlamcc.F`
differs between the working `PROPS=EOM_NLO` pipeline (`F T F`) and
scrnc's own `DENSITY=RESPONSE`-only pipeline (`F F F`, same with or
without `DENSITY=RESPONSE` present). Traced `IHBAR2`'s assignment in
`aces2/lambda/setmet.F`:
```fortran
IHBAR  = IFLAGS(43).NE.0                    ! IFLAGS(43) = HBAR keyword
IHBAR2 = IHBAR .AND. (IFLAGS(91).NE.0 .OR. IFLAGS(18).GT.100
     &   .OR. IFLAGS(18).EQ.8 .OR. .EQ.9 .OR. .EQ.10 .OR. .EQ.13
     &   .OR. IFLAGS(18).EQ.11)              ! IFLAGS(18) = PROPS keyword
     &   .AND. IFLAGS2(geom_opt).EQ.0 .AND. IFLAGS(54).EQ.0
     &   .AND. IFLAGS2(3).EQ.0
IF (IFLAGS2(151).EQ.1) IHBAR2 = .TRUE.       ! Raman override
IF (IHBAR2) IHBAR = .FALSE.                  ! avoid doing it twice
```
`PROPS=EOM_NLO` sets `IFLAGS(18)=11`, one of the qualifying values --
`DENSITY=RESPONSE` alone never touches `IFLAGS(18)` at all, so
`IHBAR2` never turns on in scrnc's own pipeline. `IHBAR2` controls
whether the converged L1/L2 (and the full effective Hamiltonian/Hbar)
get written to lists 143/144 at all -- without it, scrnc's `Scrnc_
form_mubar_*` routines were silently working from incomplete/missing
data (a real number, not a crash, hence the precision loss rather than
an error).

**Fix verified, decisively, H2O/STO-3G CCSD**: added `PROPS=EOM_NLO`
alongside `DENSITY=RESPONSE` in the ZMAT (same manual pipeline
otherwise) -- confirmed `IHBAR2` flips to `T` (`F T F`, matching the
working case exactly), and scrnc's own LR JJ matrix now agrees
RHF/UHF to **~9e-9** (was ~1.45e-4) -- a ~15,000x precision
improvement, matching the same floor seen everywhere else in this
codebase that's known to be correct.

**Ajith's guidance on how to use/generalize this finding**:
- "We need all HBAR pieces" -- confirms `IHBAR2` (build the COMPLETE
  effective Hamiltonian, not just solve Lambda) is genuinely required,
  not incidental.
- "We will find a different way to set this flag for scrnc" -- don't
  permanently require users to remember `PROPS=EOM_NLO` in their
  ZMAT; the real fix is for `scrnc_set_vars.F` (or wherever scrnc's
  own flag setup runs) to force the equivalent of `IFLAGS(18)` into a
  qualifying value (or set `IHBAR2` directly) INTERNALLY, regardless
  of what the user's ZMAT says -- not yet implemented, `PROPS=EOM_NLO`
  is a working stand-in for now.
- "EOM_NLO is for polarizabilities. For the time being we can work
  with EOM_NLO." -- acknowledged this keyword is semantically about a
  different job type (polarizability property calcs specifically),
  just a convenient way to force the needed flag for now, not the
  intended long-term mechanism.
- **"I do not think we need this for ringcc path though."** -- RCCD/
  ring-CC likely does NOT need this fix at all. Ring-CC builds and
  consumes its OWN Hbar directly (`Scrnc_form_rcc_hbar.F`, reading
  lists 56/54 -- see much earlier entries this session, already
  independently verified correct via eigenvalue checks), entirely
  separate from CCSD's Lambda-based lists-143/144 mechanism that
  `IHBAR2` controls. So RCCD's much larger H2O gap (9.58%, `JJ(4,4)`
  off 23.7x) is very likely a SEPARATE, still-unexplained bug, NOT
  fixed by this `IHBAR2`/`PROPS=EOM_NLO` discovery -- do not assume
  applying the same fix to RCCD will help without testing it
  specifically, and don't be surprised if it doesn't move RCCD's
  numbers much (or at all) if ring-CC's pipeline never reaches this
  `IHBAR2`-gated code at all.

**Checked immediately, confirmed**: RCCD's pipeline (`xringcc`/
`xringcl` run.log) never prints the `Hbar control variables` line at
all, and `grep`ing the ENTIRE `ringcc`/`ringcl` source trees for
`IHBAR2` returns zero hits -- ring-CC's Lambda-equivalent step
(`dxrcl_dwmbej.F` etc.) is completely separate code that never touches
`vlamcc.F`'s `IHBAR2` mechanism at all. **Definitively confirms**:
this `IHBAR2`/`PROPS=EOM_NLO` fix is CCSD-specific and structurally
irrelevant to RCCD -- ring-CC's much bigger H2O gap (9.58%, `JJ(4,4)`
off 23.7x) remains a genuinely separate, still fully open bug,
untouched by anything found in this entry.

**Session end, plan for 2026-08-20 (Ajith)**: tomorrow's session
focuses EXCLUSIVELY on scrnc's ring-CC (RCCD) path -- CCSD is now
understood and CONFIRMED WORKING (via the `PROPS=EOM_NLO` stand-in,
~9e-9 precision), to be used only as a reference/guidance pattern for
HOW a correct Hbar-consuming construction looks, not as something
needing further debugging itself. Concretely: resume the RCCD
response-vector trace from "continued 6" below (H2O/STO-3G, the
9.58%/23.7x gap), informed by everything learned today about IHBAR2/
Hbar completeness as a general CLASS of bug (missing/incomplete
upstream data silently causing wrong-but-plausible numbers, not
crashes) even though the SPECIFIC IHBAR2 mechanism doesn't apply to
ring-CC. Test artifacts (not committed): `h2o_ccsd_{rhf,uhf}_
sto3g_ihbar2test`, `/blue/ufhpc/perera/scrnc_test/`.

**Concretely, fallback if a permanent scrnc-internal flag fix isn't
found quickly**: pick ONE property component
(e.g. the xx one) on H2O/STO-3G (not DZP -- match scrnc's own basis
choice, or redo the EOM_NLO CI-like check on STO-3G first to have an
apples-to-apples number), and trace `Scrnc_form_mubar_00.F`/
`Scrnc_form_mubar_s.F`/`Scrnc_form_mubar_d_rhf.F` line-by-line against
`vcceh`'s own equivalent routines (`secprp.F` and whatever it calls)
to find the FIRST point where a checksum diverges between the
DENSITY=RESPONSE-triggered path and the PROPS=EOM_NLO-triggered path,
for the SAME molecule/reference/method. That divergence point is very
likely the actual, single root-cause bug behind everything chased
this entire multi-day thread (H2's ~1.4%, H2O's 9.58%, DRCCD's missing
channel may or may not be related -- still out of scope per Ajith's
standing instruction).

**Housekeeping / state as of session end**: nothing new committed this
entry (this was pure diagnosis, no code changes). Working tree still
has the local, uncommitted, test-only `scrnc_set_vars.F` RHF-RCCD-trap
override (`.FALSE. .AND. ...`) -- still needed for continued testing,
still must NOT ship. All other files match the last commit
(`498542f8`). Test artifacts from today's regroup session (not
committed, data only): `h2o_scrnc_{ccsd,rccd}_{rhf,uhf}_sto3g_regroup`,
`h2o_ccsd_{rhf,uhf}_sto3g_vcceh_clean`, `h2o_eomnlo_{rhf,uhf}_sto3g`,
`h2o_eomnlo_cilike_{rhf,uhf}_dzp` (failed, bad keyword, harmless),
`h2o_eomnlo_default_{rhf,uhf}_dzp` (the decisive CI-like-matches proof),
all under `/blue/ufhpc/perera/scrnc_test/`.

---

## SESSION 2026-08-19 (continued 2) -- terminal froze mid-investigation, recovered next session; found+fixed a real copy/paste bug in ringcl's dxrcl_dwmbej.F (Lambda-side twin of the already-fixed rcc_dwmbej.F bug), rebuilt+reran clean RHF-vs-UHF RCCD -- agreement got WORSE (0.933->3.98%), not better; root cause NOT the list-routing bug alone, still OPEN (RESUME HERE)

Picked up cold after a frozen terminal. Found substantial uncommitted
work already in place from before the freeze (not written up in a log
entry yet -- reconstructed from the working-tree diff, see below) that
had gone well beyond the "continued" entry above's stated next step.

**What was already done, uncommitted, before the freeze:**
1. `Ajith's physics call`: RHF RCCD's ring-CC=RPA equivalence needs
   BOTH the singlet (Imult=1, Hbar list 56) and triplet (Imult=2, Hbar
   list 54) channels, each with its own full Hbar/Hbar^(-1)/L/R
   response solve, combined only at the very end with the standard
   0.25 (singlet) / 0.75 (triplet) weight -- matching
   `check_rhf_rcc.F`'s `E_total=0.75*E_t+0.25*E_s` reference formula.
   Implemented across `scrnc_driver.F` (new `Imult`/`Nmult` loop around
   the existing Icount dynamic-frequency loop, with a Freq-sign-flip
   fix for the Itop=1 case), `scrnc_lrspn_rcc_main.F` and
   `scrnc_jj_rcc_lrspn.F` (both take a new `Imult` arg, route list
   numbers via `Loffset=4*(Imult-1)` so each channel gets its own
   373-376/377-380 list block), `scrnc_form_lrspn_rcc_lists.F`
   (registers 377-380 for `Iuhf=0`), and `scrnc_jj_rcc_assemble.F`
   (new `Imult`/`Wt` loop applying the 0.25/0.75 weights, and
   `Iuhf.EQ.0 .AND. Nmult.EQ.1` guard on the `Spntsing` call since the
   Nmult=2 channels are already spin-adapted at the Hbar source and
   Spntsing expects genuinely-UHF-shaped input).
2. A real bug fix in `ringcc/rcc_dwmbej.F`'s `IMULT.EQ.2` (triplet)
   branch: was reading list `LISTTB` but sizing/contracting with the
   singlet `NMSYTA`/`DSSYTA`/`LISTTA` variables instead of the TB-
   suffixed ones. Fixed to use `NMSYTB`/`DSSYTB`/`LISTTB` throughout
   that branch (both the `GETLST` call and the `RCC_TRP_CNTRCT` call).
3. `vcc.F`: `METHOD.NE.48.AND.METHOD.NE.49` guard added to the
   singles-contribution block (RCCD/DRCCD have no singles, but the old
   `METHOD.GE.9` catch-all let them through once RHF RCCD's upstream
   crashes were cleared, hitting a real divide-by-zero in
   `t2t1aa1.f:86`). Also registers list 190 explicitly before the
   `DRMOVE` call in both branches (was only registered later, inside
   `xscrnc`, too late for this call under RHF+DENSITY=RESPONSE where
   `SING1` can come out true even for singles-free ring-CC).
4. `checkhbar.F`: disabled the Hbar(AB,CI)/Hbar(IA,JK) diagnostic
   checksum blocks entirely (both UHF and RHF branches) -- ring-CC
   never builds those Hbar pieces at all (lists 107-110, 127-130), so
   the diagnostic's `IDSYMSZ` call was asserting (`ITYPL=0 ITYPR=0`)
   the first time ring-CC ever reached this routine. Diagnostic-only,
   zero effect on actual results.
5. `scrnc_set_vars.F`: the `REF=UHF` trap for RCCD/DRCCD temporarily
   forced open (`.FALSE. .AND. ...`) for testing -- NOT a permanent
   change, needs reverting to the selective CCSD-only-allow state once
   this thread resolves.

**New this session: found the same bug's twin, unfixed.** `ringcl`
(the Lambda-equation solver) has its own copy of this exact routine,
`ringcl/dxrcl_dwmbej.F`, and its `IMULT.EQ.2` branch (lines ~267-293)
had the identical copy/paste bug: `CALL GETLST(...,NMSYTA,...,LISTTB)`
and `RCC_TRP_CNTRCT(...,DSSYTA,NMSYTA,...,LISTTA,...)` where the
already-computed `DSSYTB`/`NMSYTB`/`LISTTB` should have been used.
`xringcl` had NOT been rebuilt since 2026-08-14 -- this bug was live in
every RHF RCCD run done earlier today. Fixed identically to the
`rcc_dwmbej.F` fix, rebuilt `xringcl` (`make` in `aces2/ringcl`),
installed to `aces2/bin`, md5-verified.

**Verification setup**: confirmed all four relevant binaries
(`xvcc` 13:46, `xringcc` 12:52, `xscrnc` 13:48, `xringcl` now 14:14)
are freshly built and consistent. Discarded the prior UHF comparison
run (`h2_scrnc_rccd_uhf_sto3g_fresh2`, ran at 13:47) as invalid --
it predates the 13:48 `xscrnc` rebuild by a minute, so its `@IDSYMSZ`
crash reflects stale code, not real signal.

**Ran two genuinely fresh, clean, full manual pipelines**
(`xjoda,xvmol,xvmol2ja,xvprops,xvscf,xvtran,xintprc,xringcc,xringcl,
xscrnc`, per this file's own "Practical test-running notes") in new
directories seeded from the same H2/STO-3G `ZMAT`s used earlier today,
`/blue/ufhpc/perera/scrnc_test/h2_scrnc_rccd_{rhf,uhf}_sto3g_finalcheck/`:

- Both completed with exit 0, no crashes, no assertions -- the
  `checkhbar.F`/`vcc.F` fixes hold up cleanly on a from-scratch run.
- Bare-v-independent **response-only** JJ(1,1): RHF `0.09340450647` vs
  UHF `0.06020595989` -- ratio 1.55, i.e. RHF's response piece is
  ~55% too large, not a small correction.
- Full **RPA JJ matrix** JJ(1,1): RHF `0.86801033153` vs UHF
  `0.83481178495` -- a **3.98% relative gap**, i.e. WORSE than the
  ~1.2% gap recorded in the "continued" entry above (which predates
  both this session's fixes AND the singlet/triplet-combination
  rewrite). Off-diagonal JJ(1,2)/JJ(2,2) etc. remain bit-identical
  between RHF and UHF (`0.44345130488`/`0.29626627913`), confirming
  the bare-v piece and the doubly-off-diagonal response elements are
  fine -- the gap is entirely localized to the same diagonal/(1,4)
  response elements as before.

**Conclusion: the dxrcl_dwmbej.F fix was a real bug, worth having
fixed, but it is NOT the (sole) source of the RHF/UHF gap** -- fixing
it and adding the singlet+triplet combination together made the
disagreement bigger, not smaller. For H2/STO-3G's single irrep (sym=
off), `NMSYTA` may simply equal `NMSYTB` numerically, making that fix
a correctness improvement with near-zero effect on this particular
test case -- consistent with the gap barely moving from the fix itself
and instead being dominated by something in the new singlet/triplet
combination machinery (items 1 above): the 0.25/0.75 weight assignment
itself, a sign convention, or the `Loffset` list-routing scheme are
now the prime suspects, NOT yet individually isolated.

**Next step (needs Ajith's physics judgment, not just more code-
reading)**: isolate which piece of the Imult=1/2 combination is wrong
by comparing the Imult=1-only (singlet, weight forced to 1.0, Nmult
temporarily forced to 1) response-only JJ against a hand-derivable
expectation, then add the triplet channel back in and see at which
step the 3.98% appears. Also worth double-checking the 0.25/0.75
assignment isn't swapped (0.75/0.25) -- an accidental swap would still
combine to a real number, just the wrong one, and wouldn't show up as
a crash. Not yet done this session.

**Housekeeping still needed, not done this session**: (1) revert
`scrnc_set_vars.F`'s temporary `.FALSE. .AND.` override once resolved
-- currently RCCD is wide open for RHF with a known-wrong answer,
should not ship like this. (2) nothing committed this session --
all of items 1-5 above plus the `dxrcl_dwmbej.F` fix and this log
entry remain uncommitted in the working tree.

---

## SESSION 2026-08-19 (continued 5) -- the "3.98%->1.42%" fix in the entry below was a MEASUREMENT ARTIFACT, not real progress: found a SECOND bug (`scrnc_form_rcc_hbar.F`'s Hbar write-back, also hardcoded to list 56) that exactly CANCELS the first one when both are present -- properly fixing BOTH restores the ORIGINAL 3.98% number exactly. The real RHF/UHF response gap is still 100% unexplained (RESUME HERE)

**What happened**: after the "continued 4" entry below reported
3.98%->1.42% from fixing `scrnc_rcc_hbar_read.F`'s hardcoded list-56
read, kept digging (per Ajith's "continue, goal is to fix this") and
found `scrnc_form_rcc_hbar.F` ALSO hardcodes its Hbar *write-back*
target to list 56 regardless of `Imult` (`LISTDUMP=56` unconditionally,
before this fix) -- the completed, Fock-diagonal-added triplet Hbar
was never actually being stored to list 54 at all. Fixed it
(`LISTDUMP=54` for `Imult.EQ.2`), also fixed a matching hardcoded-to-
type-56 `ISIZE` calculation in the same routine for consistency.

**Then the numbers went BACKWARD**: rebuilt, reran fresh (properly
this time, not reusing a directory -- caught myself reusing one once,
which corrupted a UHF number that should never have changed; redid it
clean) -- RHF landed EXACTLY back on `0.86801033153`, the untouched
ORIGINAL pre-any-fix number, to all 11 printed digits. Too precise to
be coincidence. Traced it with round-trip checksums (write value vs.
immediately-following read value for list 54) -- the write/read pair
IS internally consistent, no corruption there. The real explanation:
**the original code had TWO bugs, both hardcoded to list 56, that
canceled each other out for this exact call pattern.** `Scrnc_driver.F`
calls `Scrnc_form_rcc_hbar` then immediately `Scrnc_rcc_hbar_read`,
same `Imult`, back to back, every pass. In the untouched original
code: Imult=2's `Scrnc_form_rcc_hbar` correctly READ raw W from list
54, built the complete Hbar, then WROTE it to list 56 (the write bug)
-- and the very next call, `Scrnc_rcc_hbar_read`'s own list-56-only
read bug, picked up EXACTLY that freshly-written value, because both
bugs pointed at the same wrong list at the same moment. Two wrongs
made a right, by accident, only because of this specific call
adjacency. Fixing ONLY the read bug (as "continued 4" did) broke that
accidental cancellation -- Imult=2's read then legitimately hit list
54, but list 54 still only held the STALE, Fock-diagonal-incomplete
raw W (from `rcc_dwmbej.F`'s much earlier write during `xringcc`,
before `xscrnc` ever ran) -- giving a real but WRONG number
(`0.82296269100`), not a partial fix. Fixing BOTH bugs together
restores self-consistency and reproduces the untouched original
result exactly, because both paths now agree on which list holds the
real answer.

**Correct, sobering conclusion**: `scrnc_rcc_hbar_read.F`'s and
`scrnc_form_rcc_hbar.F`'s list-56 hardcoding were both real bugs
(fragile, list-56-was-never-supposed-to-be-hardcoded), and both are
now fixed correctly and are worth keeping -- but they were NEVER the
cause of the RHF/UHF response mismatch. The "continued 4" entry's
3.98%->1.42% narrative is **wrong** and should not be trusted or cited
-- it measured an inconsistent intermediate state (read fixed, write
not yet fixed), not real progress. **The actual RHF/UHF gap is back to
its original, still fully unexplained 3.98%** (confirmed via fresh,
non-reused directories:
`h2_scrnc_rccd_{rhf,uhf}_sto3g_isizefix`/`callcount`/`roundtrip`).

**Lesson for next time, worth remembering**: when two bugs are both
"hardcoded to list N" in a read/write pair called back-to-back, fixing
only one can look like huge progress (a real, large, correctly-
observed effect) while being completely wrong -- because breaking an
accidental cancellation looks identical to a real partial fix from a
single before/after comparison. Should have suspected this the moment
solving for the empirical combination weight kept giving non-physical
fractions (0.406, then 0.043) instead of converging toward something
clean -- that was the tell that the "raw" numbers being combined
weren't actually trustworthy yet, not that the weight itself was hard
to pin down.

**Follow-up, same entry's continuation**: with the two Hbar bugs
properly, self-consistently fixed (confirmed reproducible across two
independent fresh runs -- identical raw values both times, ruling out
non-determinism), pulled trustworthy raw per-channel response numbers
one more time: `R_S(Imult=1)=0.080635953018`,
`R_T(Imult=2)=0.097660690958`, UHF target `0.06020595989`. **No linear
combination `Wt_S*R_S + Wt_T*R_T` reaches the target at all** -- tried
0.25/0.75 (current, 55% too high), 0.25/0.25 (energy-style, 26% too
low), 0.5/0.5 (48% too high), and solving directly for whatever weight
would fit gives `Wt_S~2.2, Wt_T~-1.2` -- not a physically meaningful
weight (outside [0,1], wrong sign). This is qualitatively different
from earlier attempts (which at least gave weights in a plausible
range, just not a clean fraction) -- here NO weight works.

**Working hypothesis, needs Ajith's physics judgment, not more code
tracing**: UHF's own Hbar (`scrnc_rcc_hbar_read.F`'s `ELSE` branch,
`IUHF.NE.0`) is NOT block-diagonal in a "singlet channel + triplet
channel" sense -- it explicitly builds AABB/BBAA OFF-DIAGONAL coupling
blocks (`IOFFAABB`, `IOFFBBAA`) alongside the AAAA/BBBB diagonal ones,
a genuine 2x2-coupled system. A coupled matrix's EIGENVALUES can still
exactly match the eigenvalues of its symmetry-decoupled diagonal
blocks (basic linear algebra for a block matrix that's actually
block-diagonalizable in a rotated basis) -- which is exactly why the
Hbar eigenvalue check passed. But the RESPONSE FUNCTION (built from
how a perturbation vector Y projects onto and gets inverted against
that structure) does NOT have to reduce to a simple weighted sum of
two INDEPENDENT block inversions, even when the eigenvalues do decouple
cleanly. The current RHF two-channel design (`Scrnc_driver.F`'s
`Imult=1,2` loop + `Scrnc_jj_rcc_assemble.F`'s final `Wt`-weighted sum)
assumes exactly that reduction. If UHF's AABB/BBAA cross-coupling
genuinely contributes to the response (not just formally present but
numerically negligible), no amount of weight-tuning on two
INDEPENDENT channel solves can ever reproduce it -- the architecture
itself, not any single line, would need to change to include an
analogous cross-term.

**RESOLUTION of the architecture question (2026-08-19, same day,
"continued 6")**: Ajith confirmed directly -- "NO, they should not
decouple. We do not do that in CCSD route." Compared `rcc_dwmbej.F`'s
RHF two-channel construction against `vcc/dwmbej.F`'s (CCSD's, RHF,
never decouples) side by side. Found the two channels are NOT an
arbitrary decoupling -- `rcc_dwmbej.F` IMULT=1's bare-integral term is
literally `2*J-K` and IMULT=2's is `-K` (list 56, list 54
respectively) -- this is exactly the standard closed-shell Casida
singlet/triplet RPA matrix construction
(`A^singlet=2K_direct-K_exchange`, `A^triplet=-K_exchange`). So list 56
alone (`2J-K`) IS the complete, genuine singlet RPA matrix, not half
of something needing external recombination. First-principles check:
a spin-independent (density-type) perturbation operator
`V=Sum_pq v_pq(a+_pa a_qa + a+_pb a_qb)` has ZERO overlap with the
Sz=0 triplet excitation manifold (only the symmetric/singlet
combination) -- so the RESPONSE to such a perturbation should come
from the singlet channel ALONE, weight 0 on triplet, not 0.25/0.75.
Searched literature (3 angles: spin-adapted RPA, closed-shell
polarizability, GW screened Coulomb) -- confirms the general
singlet/triplet decoupling principle for eigenvalues but no exact
published formula for this specific CC-dressed ring-CC response
quantity; the derivation above is first-principles, not a citation.

**Implemented**: `Scrnc_driver.F`'s `Nmult` forced to 1 always
(dropped the `Iuhf.EQ.0.AND.Rccd`-triggered `Nmult=2` -- RCCD now
matches DRCCD's pre-existing Imult=1-only treatment). Removed the
`Spntsing` call from `Scrnc_jj_rcc_assemble.F` entirely at first (it
was only ever valid for CCSD-LR's stacked T1+T2(ABAB) shape, never
ring-CC's pure-doubles vector) -- but this made DRCCD WORSE (1.44%->
2.59%) while helping RCCD (3.98%->2.45%). Investigated why: DRCCD
(`drcc_dwmbej.F`) drops the exchange term entirely (`2J` only, direct-
RPA by definition) -- it never had a second channel to begin with, and
its own Hbar-eigenvalue check (an earlier entry today) already showed
UHF's direct-RPA Hbar still has a genuine SECOND eigenvalue (`33.927`)
that RHF's single list-56 channel doesn't reproduce -- DRCCD is
missing a real, unbuilt piece of Hbar, a bigger separate problem;
`Spntsing` was accidentally, if wrongly, compensating for some of that
gap. Restored `Spntsing` for DRCCD specifically (added the
`/Eominfo/` `Rccd`/`Drccd` common block to `Scrnc_jj_rcc_assemble.F`,
guarded `Spntsing` to `Iuhf.EQ.0 .AND. Drccd` only) -- confirmed fresh:
DRCCD reproduces its exact original baseline (`0.79307663018`)
unaffected, RCCD keeps its improvement.

**Current state, RCCD, fresh verified**: RHF `JJ(1,1)=0.85524177807`
vs UHF `0.83481178495` -- **2.45% gap**, down from the original 3.98%,
real and reproducible. Response-only piece: RHF `0.08063595302` vs UHF
`0.06020595989`, ratio ~1.339 -- not a clean/recognizable scale factor
(not 4/3, not sqrt(2), not 2). Traced further into
`Scrnc_form_mubar_s.F` (the T2-contraction dressing that builds the
480/482 perturbation lists `scrnc_lrspn_rcc_main.F` reads) -- its RHF
branch (`IUHF.EQ.0`) already correctly does the same `2*T(AB)-T(list
39)` CCSD-style combination, so this piece looks correctly
implemented on inspection; the remaining ~2.45% source not yet
isolated further. Not a clean weight or scale-factor problem this
time (unlike the earlier "1.42%" false lead) -- genuinely still open.

**FOUND IT, and it's the real one.** Instrumented `Scrnc_form_mubar_s.F`
directly (temp debug prints, since removed) to print `T(list37)` and
`T(list39)` side by side, for both H2/STO-3G and H2/DZ. Result:
**`T(list37) == T(list39)`, element for element, exactly**, in both
systems (not a coincidence of the trivial STO-3G case -- confirmed at
DZ's multi-element MATDIM=3 too). Root cause: `rcc_dwmbej.F`'s own
header comment says it outright -- "There are no antisymmetrized T2
amplitudes" for ring-CC. List 39 is built generically (via
`RNABIJ(TYPE='T')`, resorting list 46 into an `Ab<->bA`-swapped
version) for CCSD's genuinely-antisymmetrized T2, where that swap
produces a real, different quantity. For RCCD's un-antisymmetrized
T2, that same swap is a no-op -- list 39 degenerates to an exact copy
of list 37. `Scrnc_form_mubar_s.F`'s RHF dressing (shared between
CCSD and RCCD, called from both `scrnc_lrspn_cc_main.F` and
`scrnc_lrspn_rcc_main.F`/`scrnc_jj_rcc_lrspn.F`) does
`2*T(37)-T(39)`, exactly matching CCSD's own working formula -- but
for RCCD, `2*T(37)-T(39) = 2*T(37)-T(37) = T(37)`, silently DROPPING
the entire second (exchange-type) contribution the formula is meant
to subtract. This is the actual bug: CCSD's dressing formula was
reused as-is for RCCD's response without accounting for RCCD's
different (non-antisymmetrized) amplitude structure.

**Fix**: confirmed `T(list34)` (the SAME exchange-type partner
`rcc_dwmbej.F`'s own, already-verified-correct Hbar construction uses
for list 37) holds genuinely distinct, nonzero data (`0.159334637447`
vs list37/39's `-0.096723866173` at STO-3G) -- reused it.
`Scrnc_form_mubar_s.F` now branches on the `/Eominfo/` `Rccd` flag
(added the common block): `LISTT2B = 39` normally, `= 34` when
`Rccd` -- CCSD and DRCCD (deferred, separate, unaffected) keep list
39 exactly as before.

**VERIFIED, fresh, both systems, both references**:
```
                RHF/UHF gap, before -> after this fix
H2/STO-3G RCCD:  2.45%  ->  1.41%
H2/DZ     RCCD:  11.6%  ->  1.02%
```
DZ -- the system with real off-diagonal Hbar coupling, previously the
WORST case -- improved the most, strong evidence this was the right
mechanism (a bug that only bites once there's genuine multi-element
T2 data to lose, exactly matching the earlier "gap grows with system
size" observation). DRCCD reconfirmed unaffected (exact original
baseline `0.79307663018`), ground-state energy reconfirmed unaffected
and still ~1e-10 RHF/UHF agreement (this fix only touches the
response-dressing path, never `rcc_cmpeng.F`).

**Still not fully closed**: ~1-1.4% residual remains, both systems,
notably more consistent in size across systems now than before (was
2.45% vs 11.6%, wildly different scaling; now 1.41% vs 1.02%, close)
-- possibly a smaller remaining issue, possibly approaching the
`cc_conv=9`/`scf_conv=10` convergence floor (though CCSD-LR's own
~1e-8 agreement suggests floor should be tighter than 1%). Checked
`Scrnc_form_mubar_00.F` (the OTHER RHF-specific dressing piece,
computes the `MBAR00` scalar) -- looks clean, its own T1-list (90)
read is naturally zero for RCCD (no singles), no analogous bug found
there. `scrnc_jj_rcc_lrspn.F` calls the SAME (now-fixed)
`Scrnc_form_mubar_s.F`, so the L-side is automatically covered too,
no separate fix needed. Not yet root-caused further.

**Next step**: decide whether the remaining ~1-1.4% is worth chasing
further (diminishing returns territory) or is close enough to accept
as a convergence-level residual -- Ajith's call. `scrnc_set_vars.F`'s
RCCD/DRCCD RHF trap is STILL force-disabled for testing -- revert
before shipping. Nothing from today's scrnc-level response-side
changes committed yet (only the earlier, unrelated verified fixes
from `d80f2742` are committed) -- this fix (and the earlier Hbar-read/
write Imult fixes, and the Nmult=1/Spntsing architecture change) are
all real, verified improvements worth committing once the residual is
resolved or accepted.

**Superseded decision, kept for the narrative** (below was written
before the "NO, they should not decouple" resolution above; the actual
outcome differs from what this predicted):
Next session should start by thinking through the correct RHF response
formula on paper/in discussion -- specifically whether the singlet/
triplet channels genuinely decouple for THIS response quantity (they
do for the ground-state energy, confirmed today, and for the bare
Hbar eigenvalues) or whether a cross term is structurally required,
before writing any more code. Not yet resolved. `scrnc_set_vars.F`'s
RCCD/DRCCD RHF trap is STILL force-disabled for testing -- revert this
before anything ships regardless of how the architecture question
resolves.

---

## SESSION 2026-08-19 (continued 4) -- ROOT CAUSE FOUND AND FIXED: `rcc_cmpeng.F` was using the WRONG singlet/triplet recombination formula for RHF RCCD ground-state energy; corrected + verified against UHF to ~1e-10; confirmed independently by the Rishi/Perera/Bartlett JCP paper's own Eq. 49 (RESUME HERE)

Ajith's own suggested debugging method, stepwise from the ground up:
"ringcc energy must match, ringcl energies must match. If they are,
then problem is in scrnc." Checked step 1 first (skipping the
scrnc-level two-channel design entirely) -- and it immediately failed:
fresh RHF vs UHF RCCD ground-state energies (H2/STO-3G,
`h2_scrnc_rccd_{rhf,uhf}_sto3g_finalcheck`, both from "continued 3")
differed by **15 mHartree** (`-1.143309875506` vs `-1.128488116293`),
not the sub-percent response-only gap this whole thread had been
chasing -- the real problem was one level upstream of everything
looked at so far.

**Root cause, found in `aces2/ringcc/rcc_cmpeng.F`'s `RHF.AND.RCCD`
branch (~line 117)**: the routine already computes TWO candidate total
energies from the same singlet/triplet channels (`E_S`, `E_T` from
`RCC_RHF_CMPENEG`) and PRINTS both, correctly labeled:
```
E_flip     = 0.25*(E_S + E_T)        <- printed "with spin-flip excitations"
E_no_flip  = 0.75*E_T + 0.25*E_S     <- printed "without spin-flip excitations"
```
with an adjacent hardcoded comment stating `E_flip` "is identical to
the energy computed with the spin-orbital formulation" -- but a prior,
dated `2026-08-13` change set `ETOT=E_no_flip` anyway, based solely on
matching an archived reference file (`Testdir_rcc/rccd_rhf_adz.out`,
H2O/cc-pVDZ) that was apparently never itself cross-checked against an
independent UHF run.

**Two independent confirmations that `E_flip` (not `E_no_flip`) is
correct:**
1. Fresh UHF RCCD (H2/STO-3G) converges to correlation energy
   `-0.011836838064` -- matches `E_flip`'s `-0.011836838068` to 8
   decimal places, does NOT match `E_no_flip`'s `-0.026658597092`.
2. Ajith recalled a co-authored paper covering exactly this: Rishi,
   Perera & Bartlett, *J. Chem. Phys.* 153, 234101 (2020) (arXiv
   2008.00336). Its Eq. 49 gives the plain rCCD ground-state
   correlation energy as the EQUAL-weight combination
   `E_RrCCD=(1/4)(3B.3T+1B.1T)` -- exactly `E_flip`'s shape (no
   triplet-degeneracy factor). The paper's OWN 3x-triplet-weighted
   formula (Eq. 60, `E_no_flip`'s shape) is for a DIFFERENT quantity --
   spin-flip-inclusive excitation/response energies -- not the
   ground-state energy. This also means the scrnc two-channel
   0.25/0.75 response design (matching `E_no_flip`'s shape) may be the
   physically CORRECT choice for the response function specifically,
   even though that same shape is wrong for the energy -- narrows the
   still-open scrnc response gap to an implementation bug (list
   routing/sign/etc.), not a wrong-formula problem.

**Fix applied**: `rcc_cmpeng.F`'s `ETOT`/`ESPIN` switched from
`E_no_flip` to `E_flip`, old comment kept (superseded, not deleted) for
the record. Rebuilt via `make` in `aces2/ringcc` (also tried the
top-level `./build.sh aces2` first for proper multi-target
lib-propagation safety, but it fails early in the unrelated `pdens`
module -- pre-existing, not investigated, worked around by rebuilding
`aces2/ringcc` directly, which turned out to already auto-install to
`aces2/bin/xringcc`, md5-verified).

**VERIFIED, fresh directories** (`h2_scrnc_rccd_{rhf,uhf}_sto3g_
energyfix`, `/blue/ufhpc/perera/scrnc_test/`): RHF `-1.128258565382`
a.u. vs UHF `-1.128258565299` a.u. -- **agree to ~8e-11**, i.e. exactly
at the convergence floor, essentially machine-precision agreement.
This is the ground-truth check Ajith's stepwise method called for, and
it now passes cleanly.

**Ringcl checked too, per the same stepwise method**: `xringcl` links
`libringcc.a` (contains `rcc_cmpeng.o`) but had NOT been relinked since
the energy fix -- relinked (`make` in `aces2/ringcl`), reinstalled,
reran fresh. RHF `-1.128980631497` vs UHF `-1.128980631422` a.u. --
agree to ~7.5e-11, same convergence-floor quality as ringcc alone.
Both stages of Ajith's stepwise check now pass.

**Then reran the full pipeline through scrnc with everything above in
place**: the RPA JJ matrix is **UNCHANGED** --
RHF `JJ(1,1)=0.86801033153` vs UHF `0.83481178495`, still the same
3.98% gap as before any of today's fixes. Makes sense on reflection:
`rcc_cmpeng.F`'s `ETOT`/`ESPIN` only control what gets PRINTED as the
total energy (and the Jacobi convergence check) -- they never feed
back into the T2 amplitudes or Hbar that scrnc actually reads from
`JOBARC`. So this is a clean, meaningful result under Ajith's own
method, not a null result: **ringcc and ringcl are now independently
verified correct and RHF/UHF-consistent; the remaining gap is entirely
inside scrnc's own response construction**, not inherited from
upstream. This rules out an entire class of explanation (bad
amplitudes/Hbar feeding scrnc) that was never actually eliminated
before today.

**Committed** (`d80f2742`): the 5 verified fixes above (`rcc_cmpeng.F`
energy formula, `rcc_dwmbej.F`/`dxrcl_dwmbej.F` triplet dimension bugs,
`vcc.F`, `checkhbar.F`) + this log. The scrnc two-channel design and
`scrnc_set_vars.F`'s temporary test override deliberately left
uncommitted -- still open, see below.

**Tried to empirically pin down the Imult=1/Imult=2 combination
weight** (Ajith's own suggested method: find the right answer from the
UHF/RPA route, then work out the RHF adjustment) by instrumenting
`scrnc_jj_rcc_assemble.F` to print each channel's raw, unweighted
response contribution to JJ(1,1) before combining (temporary debug
print, marked for removal, not committed). Result:
```
singlet (Imult=1) raw: 0.005411120836
triplet (Imult=2) raw: 0.097660690958
UHF target:             0.06020595989
```
Neither the current 0.25/0.75 weight (gives 0.0746, 24% too high) nor
the energy's 0.25/0.25 weight (gives 0.0258, 57% too low) reproduces
the UHF target. Solving for whatever weight *would* fit
(`Wt_S~0.406`) gives a non-physical fraction -- not a match to any
textbook singlet/triplet weight. Ajith's own correction, worth
recording as method: **a single data point can't tell you whether the
weight is wrong or one of the raw channel numbers themselves is
wrong -- both are live possibilities, don't assume it's just a weight
problem.**

**Decisive follow-up test, exploiting DRCCD's simpler structure**:
DRCCD (direct RPA, no exchange) never enters the two-channel
combination at all (`Scrnc_driver.F`'s `Nmult` stays 1 for DRCCD even
under RHF -- `drcc_dwmbej.F` has no triplet list) -- so a fresh RHF vs
UHF DRCCD comparison removes the weight question entirely and tests
the shared single-channel machinery on its own. Result
(`h2_scrnc_drccd_{rhf,uhf}_sto3g_isolationcheck`, fresh dirs, all of
today's other fixes in place):
- Ground-state energy: RHF `-1.137618408411` a.u. == UHF
  `-1.137618408411` a.u. -- EXACT match, amplitudes provably fine here.
- Response-only JJ(1,1): RHF `0.01847080512` vs UHF `0.03008849944` --
  still disagree (ratio 0.61), reproducing the original ~1.2%
  total-JJ gap from the very first "continued" thread entry, now
  isolated with zero ambiguity.

**Conclusion: the bug is NOT a weight-combination problem at all.**
It's a real defect in scrnc's single-channel response-construction
code itself (`scrnc_lrspn_rcc_main.F`/`scrnc_jj_rcc_lrspn.F`, the
machinery both DRCCD and RCCD's Imult=1 pass share) -- upstream
amplitudes/Hbar are fine (energy matches exactly), the two-channel
weighting question is now moot until this more fundamental bug is
found, since even the "cleanest" single-channel case is already wrong.

**Ajith's redirect**: don't assume the code diff will show the bug --
use H2's known physics instead. (Not "CCSD=RCCD for H2", which is
false -- ring vs full CCD genuinely differ even here -- but: finish
the Aug-14 planned-but-never-run Hbar-eigenvalue cross-check.)

**Hbar eigenvalue check, finally run (fresh dirs,
`h2_rpa_hbar_eig_{rhf,uhf}_sto3g`, `excite=EOM(Sf)-RCCD`,
`eomref=RCCD`, manual pipeline through `xringcl` which does the actual
`Heff(ai,bj)` diagonalization and prints `The real eigenvalues`):**
```
UHF (one diagonalization, both channels): 25.26716566  15.09016923
RHF (two diagonalizations, singlet then triplet): 25.26716566   15.09016922
```
**Match to 7 decimal places.** This is decisive: the Hbar itself
(`rcc_dwmbej.F`, `rcc_trp_cntrct.F`, `Scrnc_form_rcc_hbar.F`) is
CORRECT for RHF -- the physical excitation spectrum agrees with UHF
exactly, RHF just gets it from two separate singlet/triplet
diagonalizations instead of UHF's one combined one, exactly as
expected. Rules out the Hbar construction entirely.

**IMPORTANT CORRECTION, found right after this entry was first
written (kept below for the narrative, but its "3.98%->1.42%" claim is
WRONG -- see the follow-up entry after this one for what's actually
true).**

**FOUND, and it's a real one: `scrnc_rcc_hbar_read.F` always read list
56 (singlet Hbar), regardless of `Imult`.** This routine (called once
per `Imult` pass inside `Scrnc_driver.F`'s `Do Imult = 1, Nmult` loop)
never had an `Imult` parameter at all -- `LISTDIAG` was hardcoded to
56. `Scrnc_form_rcc_hbar.F` correctly writes the triplet Hbar to list
54 for `Imult=2` (confirmed correct by the eigenvalue check above),
but nothing ever read it back -- the "triplet channel" response was
actually built from a second, redundant re-read of the SINGLET Hbar,
mislabeled as triplet. `Scrnc_form_rcc_hbar_inv.F`'s `MATDIM`
calculation had the identical hardcoded-to-56 gap (harmless for H2/
STO-3G since list 54 and 56 happen to share the same dimension there,
but wrong in general).

**Fixed**: added `Imult` to both `Scrnc_rcc_hbar_read.F` and
`Scrnc_form_rcc_hbar_inv.F` (`LISTDIAG=54` when `Imult.EQ.2`, else 56 --
matching `Scrnc_form_rcc_hbar.F`'s own convention), updated
`Scrnc_driver.F`'s two call sites to pass it through. Also fixed the
same gap in `scrnc_rcc_hbar_inv.F` (a thin wrapper, found to be
orphaned/uncalled dead code with its own pre-existing arg-count
mismatch predating the 2026-08-15 WHBAR/W buffer split -- left mostly
as-is, not in the active path, just added the missing `Imult` pass-
through for whoever finds it next).

**Rebuilt, reinstalled, reran fresh RHF RCCD (temporary debug print
still in `scrnc_jj_rcc_assemble.F` from the earlier empirical-weight
check, not yet removed)**: the raw Imult=2 (triplet) response
contribution changed completely --
`0.097660690958` (stale singlet reread) -> `0.062672114317` (genuinely
different data now). Full JJ(1,1): **3.98% gap -> 1.42% gap**
(`0.86801033153`->`0.82296269100` vs UHF's `0.83481178495`). Response-
only JJ(1,1): `0.04835686595` vs UHF `0.06020595989` (still ~20% off
on this piece alone, though the bare-v-dominated total is much
closer). Real, major bug, confirmed by a 3x reduction in the gap --
but NOT fully closed yet. Remaining possibilities: the 0.25/0.75
weight itself still isn't quite right even with correct data now, or
a smaller residual bug remains somewhere in the list-routing/`Loffset`
machinery. Not yet root-caused further -- next step, if continued:
check whether a different weight (solve for it the same way as
before, now with TRUSTED raw numbers) reproduces UHF exactly, which
would mean the weight was the last piece; if not, keep tracing.
`scrnc_set_vars.F`'s RCCD/DRCCD RHF trap is STILL force-disabled for
testing -- still needs reverting before this ships. Nothing from this
entry committed yet (temp debug print still in place).

---

## SESSION 2026-08-19 (continued 3) -- found a conflicting, unverifiable claim in the higher-frequency memory safety net (Claude's own `reference_scrnc_debug_log_location.md`, "update 5": "MAJOR MILESTONE... ~0.03% agreement"); could NOT reproduce it after isolating every variable; current best evidence is the gap is real, ~4-12%, and NOT explained by any single uncommitted change (RESUME HERE)

Before doing anything else this entry, re-read
`reference_scrnc_debug_log_location.md` (the memory-system safety net
this project uses specifically because this terminal keeps freezing).
It contains a LATER, never-committed "update 5" claiming all of this
was already resolved: 4 bugs fixed (the 3 in "continued 2" above plus a
`vcc.F` E4S/T2T1AA1 divide-by-zero guard), and RHF/UHF RCCD JJ(1,1)
agreeing to **~0.03%** (`0.78952789957` vs `0.78977444812`) -- a
"MAJOR MILESTONE". That update also explicitly flagged its own
uncertainty: the two-channel `Imult`/`Nmult` design in
`scrnc_driver.F` etc. was "NOT yet revisited... may or may not still
be needed" and the result was "not independently confirmed but matches
expectation" -- i.e. even the session that wrote it wasn't fully sure
before the next freeze hit.

**Tried to reproduce it and could not, despite isolating every
variable found:**
- Current full working tree (all 5 known bugs fixed, including this
  session's new `dxrcl_dwmbej.F` fix, two-channel `Nmult=2` design
  active): RHF `JJ(1,1)=0.86801033153` vs UHF `0.83481178495` --
  **3.98%**, not 0.03%.
- Same tree with the `dxrcl_dwmbej.F` fix reverted (genuinely
  reproducing the pre-this-session bug, md5-verified against the
  original Aug-14 binary): IDENTICAL numbers, `0.86801033153` /
  `0.83481178495` -- confirms the fix is numerically inert on this
  H2/STO-3G/sym=off system (single irrep, `NMSYTA` apparently equals
  `NMSYTB` here even though the code was reading the wrong list before)
  and is NOT what's blocking agreement.
- Same tree with `Nmult` forced to 1 (single-channel, i.e. reverting
  just the two-channel design update-5 itself doubted): RHF
  `JJ(1,1)=0.93587773109` -- **WORSE** (12.1% vs UHF), not better.
  Rules out "just drop the two-channel design" as the fix update-5's
  own doubt implied it might be.

**Leading hypothesis, not yet confirmed**: the original "listfix" test
directory that produced `0.78952789957` (both in this session's very
first check AND apparently behind update 5's claim) shows suspicious
timing -- its `ZMAT`/`GENBAS` are timestamped 13:14 but
`DERGAM`/`MOABCD`/`JOBARC`/etc. are all 13:43, a 29-minute gap for a
2-atom STO-3G job that runs in ~1-3 seconds end-to-end. That strongly
suggests the directory was reused for multiple partial reruns (e.g.
rerunning just `xscrnc` or `xringcl` on top of leftover `JOBARC` state
rather than a single clean `xjoda`-to-`xscrnc` pass) -- exactly the
non-idempotent-reuse pitfall this file's own "Practical test-running
notes" section warns about, and the exact failure mode update 5 itself
caught happening to an OLDER "reference" UHF directory ("the old
'successful' UHF reference run really was relying on stale directory
state the whole time"). If `listfix`'s RHF number is ALSO contaminated
by stale state, then `0.78952789957` was never a trustworthy number to
begin with, in either this session's first check or update 5's.

**Every number in THIS entry, by contrast, comes from directories
created fresh this session** (`ZMAT`+`GENBAS` copied in, immediately
run `xjoda` through `xscrnc` once, no reuse):
`h2_scrnc_rccd_{rhf,uhf}_sto3g_finalcheck`,
`h2_scrnc_rccd_{rhf,uhf}_sto3g_origringcl`,
`h2_scrnc_rccd_rhf_sto3g_nmult1test` (all under
`/blue/ufhpc/perera/scrnc_test/`). All exit 0, no crashes -- so the
disagreement is not a stale-state artifact on THIS side, at least.

**Current best-evidence conclusion**: with everything currently
believed fixed, RHF and UHF RCCD genuinely disagree by several percent
on their response-only contribution, and neither of the two design
choices tried (two-channel-weighted vs single-channel) reproduces the
close agreement update 5 claimed. That claimed result is more likely
explained by stale-directory contamination than by a real prior state
of the code that has since regressed -- but this is not proven, just
the best available explanation. **Genuinely open, needs Ajith's
steer**: worth deciding whether to (a) trust this session's
fresh-directory numbers as ground truth and treat update 5's claim as
a false positive, not chase it further, or (b) try to reconstruct
update 5's exact test conditions more carefully before giving up on
reproducing it.

**Working tree restored to a clean, understood state** before writing
this entry: `scrnc_driver.F`'s `Nmult` logic is back to its intended
`Nmult=2 If (Iuhf.EQ.0 .AND. Rccd)` form (the temporary `Nmult=1` test
edit was reverted), `xscrnc`/`xringcl` rebuilt+reinstalled+md5-verified
to match. Nothing committed. `scrnc_set_vars.F`'s RCCD/DRCCD RHF trap
is STILL force-disabled (`.FALSE. .AND. ...`) for testing -- still
needs reverting before this ships, per "continued 2"'s note.

---

## SESSION 2026-08-19 -- NEW THREAD: RHF-support question for scrnc raised, probed, and scoped; then RESOLVED for CCSD-LR (2 small bugs, not a 28-file rewrite) and verified against UHF to ~1e-8; RCCD/DRCCD left open with a real, separate ~1.2% gap (RESUME HERE, see "continued" below for the resolution)

Ajith, mid-typing when the terminal froze (recovered next session):
questioned whether scrnc needs `REF=UHF` at all, given ringcc/ringcl and
the CCSD-LR route (what scrnc actually uses) don't inherently need it --
his original reason for going UHF-throughout was that a plain RHF
diagonalization only gets you the singlet ph-manifold, and
ringcc=RPA-equivalence needs both singlet AND triplet channels summed
into the response function.

**Physics resolved**: he's right that this isn't fundamental. ringcc
already has a live spin-adapted `IUHF=0` path --
`ringcc/rcc_rng_4rhf_rccd.F` (explicit "SPIN-ADAPTED RHF CODE" branch,
`IUHF.EQ.0`) and `ringcc/rcc_trp_cntrct.F` (triplet-channel contraction,
header-commented "SPIN ADAPTED CODE", called from both
`rcc_dwmbej.F` and `ringcl/dxrcl_dwmbej.F` -- i.e. wired into both the
amplitude and Lambda sides, not dead code). So RHF-restricted
ringcc/ringcl already produces both singlet and triplet contributions
without full alpha/beta spin-orbitals.

**Engineering status, already probed this morning (before the
freeze)**: found an uncommitted comment in `scrnc_set_vars.F` dated
2026-08-19 08:32 recording that the `REF=UHF` trap there had already
been experimentally disabled and reverted. Result: NOT an artificial
holdover -- removing it segfaults in scrnc's own spin-adaptation
assembly layer (`spinad3.f`, via `spntsing.f`/
`scrnc_jj_rcc_assemble.f`), which is `Iuhf`-branched throughout (28
source files) and expects UHF-shaped spin-orbital lists. Real RHF
support needs actual rewriting of that assembly layer to consume
RHF-shaped (spin-adapted) lists, not a flag flip.

**Decision (Ajith)**: worth doing. Two payoffs -- efficiency (RHF-only
avoids carrying redundant alpha/beta copies) and, more importantly, an
independent consistency check: run a closed-shell case through both the
existing UHF path and a new RHF path and confirm the JJ-matrix/QP
results agree, per his standard exhaust-every-consistency-check method
([[feedback_consistency_check_debugging]]). Housekeeping committed
first (`6ec22620`: keeps the `REF=UHF` trap in place with the
2026-08-19 explanatory comment, plus an unrelated already-verified
change disabling `scrnc_driver.F`'s infeasible brute-force AO2MO
cross-check call). **Not started**: the actual RHF-support rewrite in
`spinad3.f`/`spntsing.f`/`scrnc_jj_rcc_assemble.f`.

**Next step**: scope exactly what shape `spinad3.f`/`spntsing.f`/
`scrnc_jj_rcc_assemble.f` currently assume (probably paired AAAA/ABAB-
style UHF list layouts) versus what ringcc's own `IUHF=0` routines
already produce/consume, to find the smallest real change that lets
scrnc's spin-adaptation layer read RHF-shaped lists instead of
UHF-shaped ones. Not yet investigated which of the 28 Iuhf-branched
files actually need changes vs. already degrade correctly to `Iuhf=0`
(scrnc's own `scrnc_built_respn_den.F`, for instance, looks like it
already sizes its beta blocks as `Iuhf*ND1BB` etc. -- may already be
RHF-clean; the segfault is specifically in the spin-adaptation trio
above, not necessarily everywhere).

**2026-08-19, later same session -- RESOLVED for CCSD-LR, NOT a 28-file
rewrite after all.** Rather than trust the morning's "real code changes
needed, not a flag flip" diagnosis, actually traced it. Two narrow bugs,
both found by reading the exact call path the segfault backtrace named:

1. **`scrnc_jj_cc_lrspn.F:289`** (the CCSD-LR path): 
   `Call Spntsing(Isympert(Iocc_pair),Work(I000),Work(Iend),Memleft)` --
   `Iocc_pair` appears NOWHERE else in the entire file (never declared,
   never assigned), so under Fortran implicit typing it's an
   uninitialized garbage integer. Two lines earlier the same routine
   correctly computes `Irrepx = Isympert(Iao_pair)`, and the sibling
   ring-CC file `scrnc_jj_rcc_assemble.F:95` uses the correct
   `Isympert(Iao_pair)` in the identical call -- classic copy/paste
   typo. Feeding garbage as an irrep index into `Spntsing`->`Spinad3`
   (which indexes symmetry COMMON blocks by that irrep) is exactly the
   kind of thing that segfaults unpredictably. Fixed: `Iocc_pair` ->
   `Iao_pair`.

2. **`scrnc_newsize.F`'s `_DEBUG_LVL0` diagnostic print** (a SECOND,
   independent bug, hit immediately after fixing #1): the real `Nsizec`
   computation (lines 16-24) already correctly branches on `Iuhf` and
   never touches the BB-specific list types (Idsymsz types 1,3/2,4)
   when `Iuhf=0` -- but the diagnostic print two lines below
   unconditionally evaluated `Idsymsz(Irrepx,2,4)` (and a second
   `ISYTYP(list44/45/46)` block) regardless of `Iuhf`, and those list
   types aren't registered under RHF, so `Idsymsz` itself asserts
   (`@IDSYMSZ: Assertion failed, ITYPL=0 ITYPR=0`). Debug-print-only bug,
   zero effect on actual results -- fixed by mirroring the same `Iuhf`
   branch the real computation already uses.

**Verified, not just "doesn't crash"**: with both fixed, H2/STO-3G CCSD
ran cleanly end-to-end under `REF=RHF` through the full manual pipeline
(`xjoda,xvmol,xvmol2ja,xvprops,xvscf,xvtran,xintprc,xvcc,xlambda,
xscrnc`). Its final `LR JJ matrix` (`JJ(1,1)=JJ(4,4)=0.83435407686`,
`JJ(1,4)=JJ(4,1)=0.50949592047`) agrees with a fresh same-build UHF
rerun (`JJ(1,1)=0.83435406798`, `JJ(1,4)=0.50949591159`) to ~1e-8 --
right at the `cc_conv=9`/`scf_conv=10` convergence floor. This is a
genuine, positive UHF-vs-RHF consistency check on the CCSD-LR route --
exactly what Ajith asked for, not just an absence-of-crash result.

**RCCD/DRCCD (ring-CC) is a separate, still-open story.** Also stopped
crashing once the trap was lifted, but its numbers do NOT agree with
UHF: H2/STO-3G DRCCD gives RHF `JJ(1,1)=0.79307663018` vs a fresh UHF
rerun's `JJ(1,1)=0.80469432450` -- ~1.2% off, far outside convergence
noise. The bare-v contribution is bit-identical between the two runs
(confirmed via the `ifile,I,J,K,L` checksum lines), which rules out a
trivial cause and localizes the gap to the response (`v*chi*v`) piece
specifically. NOT root-caused this session -- natural suspects are
`ringcc/rcc_rng_4rhf_rccd.F`'s `IUHF=0` branch or
`ringcc/rcc_trp_cntrct.F`'s triplet contraction, but neither actually
inspected yet.

**Action taken**: `scrnc_set_vars.F`'s `REF=UHF` trap is now selective
-- allows RHF for CCSD (verified), still blocks it for RCCD/DRCCD
(`Iuhf.EQ.0 .AND. (Rccd.OR.Drccd)`) with a specific error message,
rather than either a blanket allow or blanket block. Test artifacts on
disk (not in git): `/blue/ufhpc/perera/scrnc_test/
h2_scrnc_ccsd_rhf_sto3g/`, `h2_scrnc_ccsd_uhf_sto3g_retest/`,
`h2_scrnc_drccd_rhf_sto3g/`, `h2_scrnc_drccd_uhf_sto3g_retest/`.
Committed: `scrnc_jj_cc_lrspn.F`, `scrnc_newsize.F`, `scrnc_set_vars.F`
(+ their `.f` twins), this log entry.

**Next step**: root-cause the ring-CC RHF response-piece gap (start by
comparing `Response-only JJ matrix` -- i.e. before the bare-v add --
between the RHF and UHF DRCCD runs above, isolating exactly which
contraction diverges) -- not started, no urgency since CCSD-LR (the
route actually used) is already solid.

---

## SESSION 2026-08-18 (continued 27) -- NEW THREAD, separate from the Nsel/Pade numerics arc below: scrnc_sigma_driver restructured to a static/per-frequency split, verified bit-identical on 3 cases

This is a new architectural thread (see
`project_scrnc_pipeline_restructure.md` memory), not a continuation of
"continued 26"'s Nsel escalation work -- that arc is still exactly where
it left off (see below), just no longer the RESUME HERE.

**Problem**: `scrnc_sigma_driver.F` (`xscrncsigma`) reran the ENTIRE
pipeline (`xjoda...xscrnc` or `xjoda...xringcc,xringcl,xscrnc`) fresh
from `ZMAT`/`GENBAS` for every one of the `N` Gauss-Legendre quadrature
nodes, even though most of that chain is frequency-independent. Traced
which steps actually touch the `frequency` file to confirm a split
would be correct, not just faster (getting this wrong would be a
silent-correctness bug): `xjoda, xvmol, xvmol2ja, xvprops, xvscf,
xvtran, xintprc` and the ground-state amplitude step (`xvcc,xlambda`
for CCSD, `xringcc,xringcl` for ring-CC) have zero frequency/omega
references anywhere in source -- only `xscrnc` reads `frequency`
(`scrnc_set_vars.F`/`vcceh.F`). A second copy of `scrnc_set_vars.F` in
`aces2/rpa/` (gated by rpa's own `-D_RCC` switch) looked like it might
mean `xringcc` is frequency-aware too -- ruled out: `rpa/` builds a
wholly separate binary (`xrpa`), and `ringcc/GNUmakefile.src` /
`ringcl/GNUmakefile.src` link neither `librpa.a` nor any `rpa/` object,
confirmed via both `grep -il "freq\|omega" ringcc/*.F ringcl/*.F`
(zero hits) and reading the link lines directly.

**Design**: `scrnc_sigma_driver.F` now reads TWO control files instead
of one -- `static_cmds` (the frequency-independent prefix, run ONCE
before the loop) and `pipeline_cmds` (now just the frequency-dependent
step, e.g. `"xscrnc"`, run once per node, same file name as before but
repurposed). A naive "just loop `xscrnc` alone on the same directory"
implementation would NOT be safe -- this driver's own existing header
comment already documents that `xscrnc`'s CCSD linear-response path
writes dynamic list records (into `MOABCD`) that collide on a second
same-directory invocation (the `@GETLST: Assertion failed` failure
mode that originally forced the old always-rerun-everything design).
Fix: after `static_cmds` completes once, snapshot every generated file
into `.static_snapshot/`; each loop iteration clears the directory back
to that pristine snapshot (never to bare `ZMAT`/`GENBAS`, never to the
previous node's leftover state) before running `pipeline_cmds`. One
implementation snag hit and fixed along the way: `find -exec cp {}
DIR/ \;` needs the `;` shell-escaped, which is fragile through
Fortran's `EXECUTE_COMMAND_LINE` string-literal path (no backslash-
escape processing without `-assume bscc`) -- switched to `cp -p $(find
...) .static_snapshot/` instead, simpler and avoids the escaping
question entirely.

**Verified bit-identical to the pre-restructure baseline** on 3 cases
(fresh directories, e.g. `h2o_scrnc_sigma_debug_ccsd_restructure_test`,
seeded from each canonical case's own `ZMAT`/`GENBAS`/`sigma_input`,
compared via plain `diff` against `qp_energies.dat`/`sigma_c_output.dat`
copied out as baselines before this change):
- H2O/STO-3G/DRCCD (ring-CC path, sym=off): `diff` clean on both files,
  including the exact same MO 1 "needed escalation to 10 points"
  behavior as the pre-restructure run.
- H2O/STO-3G/CCSD (the path the `MOABCD`-collision comment specifically
  warns about): `diff` clean on both files, no `@GETLST` assertion,
  confirming the snapshot/restore design actually does what its own
  design comment claims.
- H2O/DZ/DRCCD (larger basis): `diff` clean on both files.
- sym=on (`h2o_scrnc_sigma_debug_symon`) NOT separately re-verified this
  session -- its `ZMAT` predates the `SAVE_INTS=ON` fix from "continued
  20" (missing that keyword), a pre-existing gap unrelated to this
  restructure; only its `static_cmds`/`pipeline_cmds` were updated in
  place (see below), not run.

**Rough timing signal** (H2O/DZ, N=4, tiny toy system so this
understates the real-world win): static prefix alone ~1.6s wall, one
`xscrnc` alone ~1.2s wall -- old design ~4x(1.6+1.2)=11.2s, new design
~1.6+4x1.2=6.4s, ~43% even here. The whole point of the split is that
the static-prefix cost stays FIXED regardless of `N` while the
per-node cost is the only one actually justified to repeat, so the win
grows with both `N` and system size -- not separately benchmarked on a
larger/more realistic system this session.

**Canonical test directories updated in place** (matching this
project's established convention, e.g. "continued 26"): the 4 existing
`h2o_scrnc_sigma_debug*` directories' `pipeline_cmds` (previously the
full chain) backed up alongside as `pipeline_cmds.bak_fullchain_<epoch>`,
replaced with the new split `static_cmds`/`pipeline_cmds` pair -- so
they stay directly runnable with the newly-installed `xscrncsigma`
rather than immediately failing its new `static_cmds not found` check.
Verification itself was done in separate `*_restructure_test` sibling
directories, left on disk as evidence (scratch, not archived, not part
of the commit).

**Build**: `ifx -i8 -o scrnc_sigma_driver scrnc_sigma_driver.F`, same
command as before; installed to `aces2/bin/xscrncsigma` (`md5sum`
matched between local and installed copies).

**Not yet done / next steps**: (1) fold this into `xaces2` as a proper
keyword-driven job type instead of an external orchestration program --
this session only did the static/per-frequency split, not the
`xaces2` integration, which was the second half of the original design
goal. (2) sym=on re-verification with the correct `SAVE_INTS=ON` ZMAT.
(3) a real-system timing benchmark now that the split is verified
correct, to get an honest performance number instead of the toy-system
estimate above.

**Committed**: `scrnc_sigma_driver.F`, this log entry. Test-directory
`pipeline_cmds`/`static_cmds` updates and `*_restructure_test` scratch
dirs left on `/blue/ufhpc/perera/scrnc_test/`, not committed (data, not
source, per project convention).

---

## SESSION 2026-08-18 (continued 26) -- Nsel escalation IMPLEMENTED + verified: closes the STO-3G MO 1 gap from "continued 25" cleanly (10 points, pilot-confirmed pole-free); full comprehensive re-verification across STO-3G sym=off/sym=on/CCSD + DZ, all MOs converge, none unconverged -- COMMITTED (459a1571)

Picked up "continued 25"'s open next step directly.

**Implemented**: `Scrnc_sigma_qp_solve`'s per-MO selection is now in
a `Do While` escalation loop -- starts at `Nsel0=6` (static + 5
nearest to `|eps_HF|`, as before); if `Scrnc_sigma_qp_bracket` can't
find a clean bracket (now using the "continued 25" data-scaled pole
guard), retries with `+Nselstep=2` more (still all selected by
nearest-distance to the SAME anchor -- local, not the huge-shared-grid
regime that caused "continued 23"'s instability), up to `Nselmax=20`
or until `Nomega` is exhausted. Added a one-line diagnostic
(`"MO n needed escalation to N points"`) so the actual count used is
always visible in the log, not just inferred. Rebuilt, reinstalled,
md5-verified.

**Result on H2O/STO-3G**: MO 1 now converges at `Nuse=10`
(`eps_QP=-20.36516`, `Im[Sigma_c]=0.0714`). Independently pilot-
verified the actual 10-point selection the driver used (confirmed via
a driver diagnostic print, since a hand-derived guess at the selection
turned out wrong the first time -- worth remembering: don't guess
which points an escalating/adaptive selection picked, have the code
say so): single sign change, `|B(N)|` stays `2.9-606` throughout (no
near-pole anywhere), tiny slope (`g'=0.067`) -- genuinely clean, not
just guard-satisfying.

**Comprehensive re-verification, all four existing test
directories, fresh `sigma_input` (2-line auto-grid format) in each**:
- STO-3G sym=off: all 7 MOs converged (MO 1 via escalation to 10, as
  above).
- STO-3G sym=on (`SAVE_INTS=ON`): all 7 converged, MO 1 needed the
  SAME escalation to 10 points and lands at `eps_QP=-20.36513` --
  matches sym=off's `-20.36516` to `3D-5`, as tight as every other MO
  pair in this comparison. Auto-grid + escalation logic is symmetric-
  vs-non-symmetric-run deterministic/consistent.
- STO-3G CCSD: all 7 converged, MO 1 converged directly at the base
  `Nsel0=6` (no escalation needed) -- a different CC level apparently
  gives a smoother `Sigma_c(z)` for this same MO, consistent with
  escalation being a genuine per-fit adaptive response, not a fixed
  requirement baked into any one MO index.
- H2O/DZ: all 14 MOs converged, none needed escalation beyond 6 this
  time (this session's grid differs slightly in extension-point
  values from "continued 24"'s manual union grid, since it's now
  auto-derived from `Emax=43.3`, so exact point placement shifted) --
  same qualitative results as "continued 24"'s per-orbital design
  (core and high-virtual both previously pilot-verified clean).

**No unconverged MOs anywhere across all four cases.** This closes the
open gap from "continued 25" -- the per-orbital-scaled auto-grid
design (universal, no manual per-system tuning) plus data-scaled pole
rejection plus adaptive point-count escalation together give a
mechanism that: (a) reaches whatever energy range a given molecule/
basis actually needs, (b) never silently accepts a spurious root
(fails loudly with a WARNING instead), and (c) recovers automatically
when the base 6-point count isn't enough, without needing to guess a
bigger number for everyone up front.

**Not yet done**: no molecule/basis LARGER or more chemically
different than this project's existing H2O/STO-3G and H2O/DZ test
cases has been tried -- "universal" is a design intent verified
against two related test systems, not proven across a wide chemical
range. `Nsel0=6`, `Nselstep=2`, `Nselmax=20`, the `20x`/`0.01` pole-
guard constants, and the `1.3x`/`x1.6` auto-grid constants are all
reasonable-looking but essentially engineering judgment calls, not
independently derived or swept.

**Committed**: `459a1571` ("scrnc: escalating per-MO Nsel selection
closes STO-3G MO 1 convergence gap"), 2026-08-18. Changed:
`scrnc_sigma_driver.F`/`.f` (auto-grid, per-orbital escalating
selection, scaled pole guard -- the full "continued 24"-"continued
26" arc lives in this one file's diff) and this log. `sigma_input` in
all four test directories reduced to the 2-line auto format (old
versions backed up alongside with timestamped names, e.g.
`sigma_input.bak_<epoch>`) -- scratch, not archived, not part of the
commit.

**Next**: consider archiving updated `qp_energies.dat`/logs for all
four cases to `sigma_debug_archive/` as the new baseline; then try a
molecule/basis larger or more chemically different than H2O/STO-3G
and H2O/DZ to actually test the "universal" grid/escalation design
claim (see "Not yet done" above).

---

## SESSION 2026-08-18 (continued 25) -- auto-scaled universal Omega_ext grid IMPLEMENTED (no more per-system hand-tuning); bracket-search pole guard tightened from a fixed absolute cutoff to a data-scaled one (catches a real, previously-missed spurious-root case); net effect: STO-3G MO 1 now HONESTLY reports unconverged instead of silently returning a wrong answer -- new open gap, not yet resolved (RESUME HERE)

Ajith's question after "continued 24": a hand-typed `sigma_input`
Omega_ext range is a real limitation for other molecules/basis sets --
"we need a grid that works more or less universally."

**Implemented auto-grid generation** in `scrnc_sigma_driver.F`:
`sigma_input` now only carries `N`/`Omega_c` (2 lines, the
Nomega/Omega_ext block removed). Right after node 1's `wcreduc.dat`
gives the real `Eps(1..Nmo)` (before that node's own fold step runs),
the driver builds its own grid: the proven `{0,0.3,0.6,0.9,1.2,1.5}`
valence sub-grid, then geometric extension (x1.6/step) until exceeding
`1.3*Max(|Eps(p)|)` for THIS run's own orbitals -- zero extra pipeline
cost either way. Verified on DZ (Emax=43.3 from MO 14, the highest
virtual, not the core -- auto-reached to 64.4, Nomega=14) and STO-3G
(Emax=20.2, Nomega=13): MOs whose own energy was already within the
fixed valence sub-grid are bit-identical to "continued 24"'s
hand-picked-grid results (provably a no-op for them); DZ's core (p=1)
and highest virtual (p=14) both independently pilot-verified clean
(no pole, in-domain, small slope).

**But STO-3G's own core orbital (p=1, `eps_HF=-20.24`) came back
wrong**: `eps_QP=-20.470`, `Im[Sigma_c]=+0.681` -- two orders of
magnitude out of scale vs. every other converged MO (`~0.005-0.15`)
and vs. every actual SAMPLED value for this same MO
(`Im` maxes out at `0.0036` across its 6 fit points). Pilot-code check
confirmed a real, if moderate, nearby pole (`min|B(N)|~0.26` over a
wide window, one lone sign change, right where the driver's root
landed) -- the existing bracket-search's `|Sigma_c|<1000` pole guard
(fixed absolute cutoff, unchanged since "continued 21") completely
missed this since `0.68 << 1000`, even though it's obviously spurious
relative to the actual data scale.

**Fixed the guard, not the symptom**: `Scrnc_sigma_qp_bracket` now
takes an `Fscale` argument (`Max(|Fpt|)` over the MO's own 6 selected
sample points, computed by the caller) and rejects any candidate
bracket where `|Sigma_c|` at either endpoint exceeds
`Max(20*Fscale,0.01)` instead of the old blanket `1000` -- scaled to
what's actually physically plausible for THIS fit, not an arbitrary
global number. Also switched the check from `|Re[Sigma_c]|` to the
full `|Sigma_c|` (catches Im blowing up too, which is what actually
happened here). Rebuilt, reinstalled, md5-verified.

**Result: STO-3G MO 1 now correctly reports `converged=0` ("no
bracketed root found for MO 1") instead of silently returning the
wrong answer** -- the intended behavior (fail loudly, per this
project's standing practice) but NOT a fix for the underlying gap:
this particular MO's 6-point selection genuinely doesn't have a
trustworthy nearby root with the current design. All other STO-3G MOs
(2-7) unaffected -- confirms the tighter guard doesn't false-positive
on well-behaved fits.

**Open question, not yet resolved**: DZ's core orbital (very similar
`eps_HF~-20.6`, similar-shaped 6-point selection) worked fine;
STO-3G's core orbital (`eps_HF~-20.2`) doesn't. `Nsel=6` is
apparently not uniformly sufficient -- sometimes it is, sometimes it
isn't, and it isn't obvious in advance which. Natural next step
(not yet tried): automatically escalate `Nsel` (e.g. retry with 8,
then 10 local points, still all selected near that MO's own energy --
NOT the same failure mode as "continued 23"'s huge-shared-grid
instability, since these stay local) when the 6-point bracket search
comes back empty, before giving up.

**Committed**: nothing yet. Changed: `scrnc_sigma_driver.F`/`.f`
(auto-grid generation, `Fscale`-relative pole guard). `sigma_input`
in both `h2o_scrnc_sigma_debug{,_dz}` reduced to the new 2-line
format (scratch, not archived/committed).

**Next**: try automatic Nsel escalation on bracket-search failure;
re-verify sym=on and CCSD paths with the new auto-grid once STO-3G's
core orbital is resolved; then get sign-off to commit the whole
"continued 24"/"continued 25" arc together.

---

## SESSION 2026-08-18 (continued 24) -- per-orbital-scaled small grid IMPLEMENTED + verified: fixes DZ p=1 (core) cleanly, one bug found+fixed along the way (must always keep the omega=0 anchor point), all 14 DZ MOs now converge to sane, pole-free, in-domain roots (RESUME HERE)

Ajith's call after "continued 23"'s diagnosis (regularized/robust fit
vs. per-orbital small grid): per-orbital grid -- smaller lift, and the
data isn't actually noisy (CC-response Sigma_c is accurate to
1e-9-1e-14; a robust/least-squares fit solves a problem we don't have)
-- and "might as well do it now" rather than defer.

**Implemented in `Scrnc_sigma_qp_solve`**: `sigma_input` now carries
one modest UNION grid (13 points, `0..25`, same as "continued
22"/"23", cheap since it's still just algebraic folds of the existing
`Wc_reduced`). For each MO, instead of fitting all Nomega points, the
routine now selects only `Nsel=6` of them -- the omega=0 (static)
point ALWAYS included, plus the 5 nearest to that MO's own
`|eps_HF|` by simple insertion selection -- sorts them back into
ascending-omega order, and runs the exact same (unchanged)
`Scrnc_sigma_pade_coeff`/`Scrnc_sigma_qp_bracket` on that 6-point
subset. This automatically recovers the original `[0,1.5]` 6-point
fit for every valence/virtual MO (already its own 6 nearest points)
while giving deep orbitals their own locally-relevant subset.

**Bug found and fixed during verification, not glossed over**: first
version (distance-only selection, no forced inclusion of omega=0)
made MO 2's fit badly unphysical -- `Im[Sigma_c(eps_QP)]=+2.94`,
two orders of magnitude out of scale vs. every other MO (`~0.01-0.16`).
Root cause: MO 2's `Anchor=1.362` is roughly equidistant from `0.0`
(dist 1.362) and `2.0` (dist 0.638) but `2.0` is numerically closer,
so pure nearest-6 selection silently swapped `2.0` in for `0.0`
-- and losing the physically-meaningful static point destabilized the
fit even though every selected point was individually "close" to the
target. Confirmed via the same independent pilot-code check
(`data_selected` vs `data_orig` in scratch): min`|B(N)|` over a wide
window dropped from `0.58` (with 0.0) to `0.24` (without) -- not
technically AT a pole, but close enough that the fitted value at the
root went wildly out of scale. Fixed by always forcing the omega=0
(or whichever `Omega_ext` value is smallest) index into the selection
first, then filling the remaining `Nsel-1` slots by nearest-distance
among the rest. Rebuilt, reinstalled, md5-verified both times.

**Verified on H2O/DZ (all 14 MOs)**:
- MOs 3,4,5,6,7,8,9,10 -- IDENTICAL to the original 6-point-grid
  rootfind baseline (their nearest-6 already was `[0,1.5]`, so this
  change is provably a no-op for them). No regression.
- MOs 2,11,12,13 -- small, sane shifts (their nearest-6 swaps one
  distant `[0,1.5]` point for a slightly-closer wider-grid point,
  e.g. MO 13 trades `0.3` for `2.0`); MO 2 specifically now gives
  `Im[Sigma_c]=-0.159`, back in normal range after the static-anchor
  fix.
- **MO 1 (core) -- the original target**: `eps_QP=-20.74686`,
  `Im[Sigma_c]=-0.113`. Independent pilot-code check on its own
  6-point selection (`0,5,8,12,18,25`): single clean root, no nearby
  pole (`min|B(N)|~0.92`, nowhere near zero over `[-30,10]`), root
  INSIDE the sampled domain, tiny slope (`g'=0.046`) -- a genuinely
  well-conditioned, trustworthy fit, unlike every prior attempt at
  this orbital (6-point extrapolation: untrustworthy by construction;
  13/36-point all-points fit: right ballpark but itself part of the
  now-discredited "more points" family). Lands within ~0.03-0.06 Ha
  of the interpolated-but-unstable wide-grid answers from "continued
  22"/"23" (`-20.69` to `-20.82` across those), which is reassuring
  convergence given none of those were individually trustworthy on
  their own.

**Not yet done**: did not re-verify STO-3G (only DZ tested this
session) or re-run sym=on with the new grid; `Nsel=6` and the 13-point
union grid's specific values (`0..25`) are reasonable defaults but not
rigorously optimized/swept themselves (would need the same kind of
small-N stability study flagged as unexplored in "continued 23").

**Committed**: nothing yet -- pending Ajith's go-ahead now that this
is implemented and verified (per "wait to fix convergence issues").
Changed: `scrnc_sigma_driver.F`/`.f` (the per-orbital selection).
Scratch (not archived): `h2o_scrnc_sigma_debug_dz/sigma_input` (holds
the 13-point union grid), `qp_energies_perorbital{,_v2}.dat`.

**Next**: get sign-off to commit; re-verify STO-3G + sym=on with the
new grid/selection logic before calling this fully done; consider
whether `Nsel`/grid-point choices need a small stability sweep of
their own.

---

## SESSION 2026-08-18 (continued 23) -- N-convergence study DONE: the eps_QP scatter across grid sizes is NOT double-precision roundoff (ruled out via 100-digit arbitrary-precision arithmetic reproducing the double-precision answer exactly) -- it's the well-known intrinsic ill-posedness of high-order exact-interpolation Thiele/Pade continuation; the fix is a small per-orbital-scaled grid, not one large shared grid (RESUME HERE)

Picked up "continued 22"'s open next step: an actual N-convergence
study, per Ajith's "wait to fix convergence issues" instruction.

**Swept N = 8, 12, 18, 26, 36** for the H2O/DZ case, using a graded
grid `omega_i = Omega_c*(1+x_i)/(1-x_i)` with `x_i` evenly spaced in
`(-1,1)` (same functional form as the codebase's own
`Scrnc_sigma_freq_grid`, so N naturally grows both denser AND wider at
once -- covers up to `|omega|~71` by N=36, comfortably past even MO
14's `eps_HF=43.3`). (Process note: had to generate the grid files and
copy results OUT of the working directory each iteration -- the
driver's own per-node `find ... -delete` cleanup, added in "continued
14", wipes anything not named `sigma_input`/`ZMAT`/`GENBAS`/
`pipeline_cmds`/`*.log` at the start of its NEXT invocation, so
same-named scratch files silently got deleted mid-loop on the first
attempt; fixed by keeping generated grids and per-N output copies
outside the pipeline directory entirely.)

**Result: eps_QP does NOT converge smoothly as N grows -- it
scatters, non-monotonically, by amounts up to ~0.25 Ha** (e.g. p=12:
1.294->1.233->1.487->1.323->1.258 across N=8..36; p=2: -1.457->
-1.121->-1.572->-1.341->-1.477). This looked at first like classic
numerical instability, but which kind mattered: double-precision
roundoff (fixable with more precision) vs. genuine ill-posedness of
exact high-order rational interpolation (not fixable by precision
alone).

**Diagnosed which, using mpmath arbitrary-precision arithmetic**
(installed this session, `pip install --user mpmath`) -- reimplemented
the identical Thiele algorithm a third, independent way, and evaluated
p=1's and p=2's N=36 fit at 15, 30, 50, and 100 decimal digits.
**Result: bit-identical root at every precision level, matching the
double-precision Fortran/earlier-Python result almost exactly.**
This RULES OUT double-precision roundoff as the cause of the N-scatter
-- each individual N's fit is computed to full accuracy, precision is
not the problem. **The scatter is therefore the well-documented
intrinsic behavior of exact-interpolation Thiele/Pade continued
fractions at high order**: an exact rational interpolant through many
points is only weakly constrained between/beyond its nodes (a
Runge's-phenomenon-like effect for rational rather than polynomial
fits), so its near-real-axis pole/zero structure -- and hence the
extracted on-shell root -- can shift substantially with each added
point, even in exact arithmetic. This matches long-standing cautions
in the numerical-analytic-continuation literature (Vidberg-Serene
1977's own original method; robustness studies going back to
Beach-Gooding-Marsiglio 2000) that naive many-point Thiele
continuation becomes increasingly non-unique/ill-posed well before
double-precision roundoff would ever matter -- this is NOT a bug in
this codebase, and NOT fixable by adding more sample points, more
precision, or a better root-finder (the bracket/bisection solver from
"continued 21" is doing its job correctly on whatever fit it's handed
-- the fit itself is the moving target).

**Conclusion, reframing "continued 22"'s open question**: the right
design is NOT "one large shared Omega_ext grid, more points = more
trustworthy" (this session disproves that) and NOT "the original
6-point grid was basically fine, just check convergence" (also
disproved -- p=2/p=12/etc's 6-point-ish low-N answers scatter just as
much as the high-N ones do; there is no evidence ANY of the N tested
is "the converged answer"). The literature's usual fix applies
directly: keep N small (order ~6-10, where this scatter is presumably
tamer, though NOT yet checked systematically at low N the same way)
and place those few points to BRACKET each orbital's own energy scale
-- i.e. "continued 15" item 3's original per-orbital-scaled-grid
proposal, now with a concrete reason (not just "wider domain is nice
to have") for why it beats a single big shared grid. An unexplored
alternative also worth naming: a regularized/robust rational fit
(least-squares over more points than free parameters, or an
SVD-truncated Thiele variant) instead of exact interpolation --
standard in the GW/DMFT analytic-continuation literature for exactly
this reason, but a bigger implementation lift than a per-orbital grid.

**Not done**: did not yet test whether SMALL N (6-10) is itself
actually stable/convergent when scanned the same way (e.g. N=4,6,8,10
all bracketing a given orbital) -- assumed favorable by literature
precedent but not verified in this codebase; did not implement any
per-orbital grid logic; did not touch production code this session
(all experimentation stayed in scratch `sigma_input`/copied-out
results, per Ajith's "wait to fix convergence issues" -- nothing
committed).

**Committed**: nothing yet. Changed (scratch only, not archived/
committed): `h2o_scrnc_sigma_debug_dz/sigma_input` (currently holds
the last-swept N=36 grid; `sigma_input.orig_6pt` backup from
"continued 22" untouched). Scratch files (not in repo):
`sigma_input_N{8,12,18,26,36}`, `qp_energies_N{8,12,18,26,36}.dat`,
`p1_N36_raw.txt`, `p2_N36_raw.txt`, `pade_highprec_probe.py`.

**Next**: get Ajith's read on which fix direction to pursue
(per-orbital-scaled small grid vs. a regularized/robust rational-fit
algorithm) before implementing anything in production; if per-orbital
grid, first verify small-N (4-10) stability systematically (the one
open empirical gap above) before trusting it as the final design.

---

## SESSION 2026-08-18 (continued 22) -- widened Omega_ext grid genuinely fixes DZ p=1's extrapolation-domain problem, but reveals a DEEPER, more general issue: even previously "clean" valence/virtual MOs' eps_QP are NOT converged w.r.t. Pade sample count -- NOT resolved, needs Ajith's numerics call (RESUME HERE)

Picked up "continued 21"'s open next step (widen `Omega_ext` for the DZ
core orbital). Per the driver's own architecture (`Wc_reduced` reused
across any `Omega_ext(J)`, no extra pipeline cost), extended
`h2o_scrnc_sigma_debug_dz/sigma_input` from 6 points (`0..1.5`) to 13
(`0,0.3,0.6,0.9,1.2,1.5,2,3,5,8,12,18,25` -- original backed up as
`sigma_input.orig_6pt`). Reran `xscrncsigma`: all 14/14 MOs still
converge, no bracket-search failures.

**p=1 (core) -- genuinely fixed this time**, confirmed via the same
independent pilot-code check: single clean root at `z=-20.692`, now
INSIDE the sampled `[0,25]` domain (a real interpolation, not
extrapolation), no nearby pole (`|B(N)|~5-7D5`, nowhere near zero),
tiny slope (`g'=0.036`). But the value moved substantially from the
6-point result: `eps_QP` went from `-20.552` (correction `+0.007` Ha)
to `-20.692` (correction `-0.133` Ha) -- a SIGN FLIP in the QP shift.
Confirms the original "continued 13"/"continued 20" suspicion: the
6-point extrapolated answer was not just formally out-of-domain, it
was actually wrong.

**But p=2 and p=13 (valence/virtual, NOT domain-violation cases --
both roots were already inside the old `[0,1.5]` sampled range even
with only 6 points) ALSO shifted substantially** when the same wider
grid was used to build their Pade fit: p=2's `eps_QP` moved
`-1.335->-1.377`, p=13's moved `1.692->1.745` -- shifts of the same
order as the corrections themselves, not small refinements. Checked
both with the pilot code: both remain single clean roots, no nearby
pole, well-conditioned slopes (`g'=-0.25`, `g'=0.12`, both safely
`|g'|<1`) -- so this isn't a new instability being introduced by the
wider grid. **The likely explanation: the original 6-point/N=6 Thiele
fit was never actually verified to be CONVERGED in N** -- i.e. this
project has been implicitly assuming 6 sample points is enough to
pin down the analytic continuation accurately, but never tested that
assumption by comparing against a higher-N fit, for ANY MO. Adding
more points changing the valence-orbital numbers by O(0.01-0.05) Ha,
comparable to the physics differences reported as significant in
"continued 14"'s CC-vs-dRPA comparison, means **every eps_QP number
produced so far (STO-3G included) has an unquantified, potentially
non-negligible Pade-order truncation error that has never been
checked** -- this is a more fundamental gap than the core-orbital
domain issue that motivated this session's experiment.

**Not resolved, deliberately not patched further**: this needs an
actual N-convergence study (rerun with several different N/Omega_ext
choices -- e.g. 6, 10, 13, 20 points -- and check whether eps_QP
stabilizes as N grows, the standard way to establish a numerical
result is converged) before trusting ANY previously-reported eps_QP
number, including the already-published-in-this-log STO-3G MO-5 case
from "continued 21". Whether a fixed shared grid (this session's
approach) or a per-orbital-scaled grid (closer to "continued 15" item
3's original wording) is the right design also still needs Ajith's
numerics call, now sharpened by this finding.

**Committed**: nothing yet this session. Changed (not yet archived to
`sigma_debug_archive/` or committed): `h2o_scrnc_sigma_debug_dz/
sigma_input` (widened, original backed up alongside as
`sigma_input.orig_6pt`), `xscrncsigma_widegrid.log`/`qp_energies.dat`
in that directory. Pilot script (scratch, not in repo):
`pade_dz_widegrid_probe.py`.

**Next**: get Ajith's read on the Pade-order convergence question
before doing anything else with Sigma_c/QP numbers -- this is now the
central open numerics question for the whole pipeline, superseding
the narrower "core orbital domain" framing this session started with.

---

## SESSION 2026-08-18 (continued 21) -- Pade/QP-solve robustness: root cause was the naive fixed-point iteration, NOT a bad Pade fit; replaced with bracketed bisection; H2O/STO-3G MO 5 now converges (matches independent pilot-code prediction exactly), sym=on/off still consistent, DZ p=7 now genuinely clean, DZ p=1 (core) "converges" but is STILL an unresolved extrapolation-domain artifact -- not silently treated as fixed (RESUME HERE)

Picked up "continued 20"'s open next step, per Ajith's own framing:
diagnose the Pade robustness gap (H2O/STO-3G MO 5 hitting the 51-iter
cap post-fold-fix) "reasoning from math, not tuning to what looks
right." Per [[feedback_pilot_code_debugging]]/
[[feedback_consistency_check_debugging]].

**Pilot-code diagnosis (independent Python reimplementation of the
exact Thiele/Vidberg-Serene algorithm, reusing the driver's own 6
`(z_i,f_i)` sample points verbatim)**: scanned `F(z)=z-eps_HF-
Re[Sigma_c(z)]` and the continued fraction's denominator `B(N)(z)`
over real z for MO 5 (`eps_HF=-0.3912734`). Found exactly ONE clean
real root, `z=-0.49680`, `Sigma_c(root)=-0.1055-0.1518i`, with no
nearby Pade pole (`|B(N)|` stays order-0.1-1 throughout). But the
local slope `g'(root)=dRe[Sigma_c]/dz` at that root is **-2.14**.
Since fixed-point iteration `eps_new=eps_HF+Re[Sigma_c(eps_qp)]` only
converges locally when `|g'(root)|<1`, this root is mathematically
unreachable by the original iteration scheme -- reproducing the exact
production algorithm in the pilot confirmed it settles into a stable
period-2 limit cycle (-0.408 <-> -0.562, averaging near but never
landing on -0.497), exactly matching the real 51-iteration
non-convergence. **Conclusion: MO 5's Pade continuation was fine all
along; only the iteration scheme was inadequate** -- a different
failure mode than the DZ p=1/p=7 concerns from "continued 13"
(extrapolation-domain violation / sparse-grid kink), and exactly what
the external review's point 4 ("continued 15") already anticipated
(replace fixed-point with a bracketed root-find).

**Fix implemented**: `scrnc_sigma_driver.F`, `Scrnc_sigma_qp_solve` --
replaced the fixed-point loop with a new `Scrnc_sigma_qp_bracket`
subroutine: scans an expanding window centered on `eps_HF` (starting
+-2 Ha, doubling up to 8x) for the sign change of `F(z)` closest to
`eps_HF`, rejecting any interval where `Re[Sigma_c]` itself blows up
(`>1000`, the signature of straddling a pole rather than crossing a
genuine root), then bisects to `Tol=1D-8`. Cost is negligible (Pade
eval is already just an O(N) continued-fraction recurrence). Rebuilt
(`ifx -i8`), installed to `aces2/bin/xscrncsigma`, md5-verified.

**Verified, not just plausible**:
- H2O/STO-3G sym=off: MO 5 now converges to `eps_QP=-0.4967992,
  Im=-0.1517794` -- matches the independent pilot prediction to 4
  decimal places. All other 6 MOs reproduce their prior fixed-point
  results bit-for-bit (no regression). 21 bisection iterations for
  every MO (expected -- fixed iteration count from the doubling+
  bisection design, not a red flag).
- H2O/STO-3G sym=on (`SAVE_INTS=ON`): all 7 MOs now converge too, and
  still match sym=off by `eps_HF` to ~1e-4-1e-6 (one MO, the same
  MO-5/4 pair, is ~6e-4 -- an order of magnitude looser than the rest
  but still small; plausibly amplification of the usual
  `cc_conv=9`-scale truncation noise near this particular root's
  steep slope, not investigated further).
- H2O/STO-3G CCSD path: all 7 MOs converge cleanly too.
- H2O/DZ (larger basis, the original "continued 13" case): 14/14 MOs
  now report `converged=1` (previously 12/14, with p=1 and p=7
  failing). **Did NOT stop at "it converged" -- re-ran the same
  pilot-code pole/slope analysis on both previously-failing MOs**:
  - **p=7 (`eps_HF=0.311`)**: genuinely clean now -- single isolated
    root at `z=0.3182`, `|B(N)|` stays ~0.39 throughout (no nearby
    pole), root sits INSIDE the sampled `|omega| in [0,1.5]` domain,
    tiny slope (`g'=0.0085`). The original "sharp kink" concern
    doesn't reproduce under the fold-fixed Sigma_c samples -- this MO
    is resolved for real, not just numerically quieted.
  - **p=1 (`eps_HF=-20.559`, O 1s core)**: also finds a single,
    pole-free root (`z=-20.552`, small `|g'|=0.007`) -- BUT that root
    sits at `|z|=20.55`, ~14x past the sampled `|omega| in [0,1.5]`
    domain edge, i.e. still a severe Pade EXTRAPOLATION, not
    interpolation. The extrapolated `Sigma_c` there
    (`Re=+0.0069,Im=-0.053`) is qualitatively different from the flat,
    small, same-sign trend seen across the actual sampled points
    (`Re` stays `~-0.0187` to `-0.0188` the whole way) -- exactly the
    kind of sign/magnitude discontinuity a rational-function
    extrapolation produces when asked far outside its fitted range,
    with nothing in the fit's own internal structure (no pole, small
    slope) able to flag it as untrustworthy. **This MO's
    "convergence" does NOT mean the extrapolation is trustworthy --
    the original "continued 13" domain-of-validity concern for core
    orbitals is UNRESOLVED, just no longer visible as an iteration
    failure.** Deliberately not claimed as fixed.

**Not done / left open, per standing "not silently patched" practice**:
the core-orbital extrapolation-domain gap needs the already-scoped fix
from "continued 15" item 3 (widen/rescale the external `Omega_ext`
evaluation grid per orbital energy scale -- e.g. sample out to ~20+ Ha
for deep-core states -- reusing the SAME already-computed inner
omega'-quadrature `Wc_reduced` per the driver's own N-node
architecture, no expensive pipeline rerun needed) rather than trusting
this extrapolation. Not started this session.

**Housekeeping note**: `git status` shows several scrnc files
("continued 17"-"continued 20"'s fold-fix, `Scrnc_get2ints.F` guard,
`Scrnc_jj_unpack_natural.F`, etc.) as still MODIFIED/untracked despite
earlier log entries saying "Committed" -- those commits apparently
never actually happened (or were reverted). Not addressed this
session (no commit requested); flagging so a future session doesn't
assume that work is safely in git history when it isn't yet.

**Committed**: nothing yet this session (not requested). Changed:
`scrnc_sigma_driver.F` (the bracket/bisection rootfind); this log
entry. Rerun outputs left in place, not archived to
`sigma_debug_archive/` yet: `xscrncsigma_rootfind.log`/
`qp_energies.dat` in each of `h2o_scrnc_sigma_debug{,_ccsd,_dz}` and
`h2o_scrnc_sigma_debug_symon_saveints`. Pilot scripts (scratch, not in
repo): `pade_mo5_probe.py`, `pade_dz_probe.py`.

**Next**: (1) get Ajith's read on whether/how to commit the backlog of
already-verified-but-uncommitted scrnc work (fold fix, save_ints
guard, this session's rootfind); (2) if the core-orbital extrapolation
gap matters for the manuscript's basis-set scope, widen the
`Omega_ext` grid per "continued 15" item 3; (3) the still-open
bracketed-root-find endorsement also separately resolves the "is the
continuation trustworthy vs does the QP equation have a stable root"
question the external review raised -- worth re-reading "continued
15" in full now that this piece is done.

---

## SESSION 2026-08-18 (continued 20) -- sym=on bug FULLY RESOLVED: second root cause was SAVE_INTS=ON missing, not a code bug at all; both fixes verified together, sym=on now matches sym=off to ~1e-4 or better across all 7 H2O/STO-3G MOs; standing FATAL guard added so this can't silently recur (RESUME HERE)

Picked up "continued 19"'s open residual (~0.11 max diff after the Iao_pair
permutation fix). Ajith's steer: "scrnc_prep_2ints calls scrnc_get2ints.
When sym=on all 4 files have integrals (abelian symmetry)... similar logic
can be found inside acesii" -- i.e. don't assume this is scrnc's own bug,
check the shared ACES2 machinery first.

**Traced the file lifecycle precisely**, stepping the pipeline one binary
at a time and checking `IIII`/`IIJJ`/`IJIJ`/`IJKL` existence after each:
present after `xvscf`, GONE after `xvtran`. Root cause, found via two
background research passes into `aces2/vmol/` and `aces2/vtran/`:
- `aces2/vmol/readin.F` (~line 519) and `aces2/vmol/twoel.F` (~line 45):
  `NSYMOP` (number of Abelian symmetry generators) gates which of the 4
  files even get OPENED by `xvmol` -- 0=only `IIII`, 1=+`IJIJ`/`IIJJ`,
  `>=2`=all four (C2v needs 2 generators for its 4 irreps, satisfied here;
  confirmed directly in the pipeline log: `xvmol` itself reports "45
  integrals of symmetry type I J I J", "39 ... I I J J" -- the data is
  genuinely computed).
- `aces2/vtran/load2.f` (lines 159-173), called from the standard
  `TRIIJJ`/`TRIJIJ`/`TRIJKL` AO->MO transform routines: once `xvtran` has
  consumed a raw integral file for the transform, it does
  `close(unit=iUnit,status='DELETE')` -- UNLESS one of five existing
  `IFLAGS`/`IFLAGS2` conditions holds, in which case `status='KEEP'`.
  Identified via a second background research pass (mapping indices to
  ZMAT keywords via `aces2/joda/asv_ctl.c`): `IFLAGS(26)` =
  `h_IFLAGS_save_ints`, ZMAT keyword **`SAVE_INTS=ON`** -- literally the
  existing, user-facing "keep the AO integral files" switch (the other
  four conditions are PROPS>=NMR-class, ABCDTYPE=AOBASIS, open-shell-
  reference-specific, and TDHF=ON -- none apply to scrnc's own runs).

**So the second bug isn't a scrnc code bug at all** -- `scrnc_get2ints.F`
reads the SAME raw integral files every other ACES2 "light"/pilot module
(`ccsd_light/get2ints.F`, `ccd_light/get2ints.F`, `cse/get2ints.F`,
`cse12/get2ints.F`) also reads, all sharing this identical 4-file/8-fold-
permutational-unfold pattern -- it's a standard, working ACES2 idiom.
scrnc's own pipeline (`xjoda;xvmol;...;xvtran;xintprc;xringcc;...;xscrnc`)
just needs `SAVE_INTS=ON` in the ZMAT, same as any other module that reads
raw AO integrals after `xvtran` has already run, which scrnc uniquely does
(most modules consume AO integrals BEFORE or DURING the transform, not
after).

**Verified together, decisively**: built `h2o_scrnc_sigma_debug_symon_saveints`
(same ZMAT as "continued 18"/"19"'s sym=on case, `+SAVE_INTS=ON`). Confirmed
all 4 files now survive past `xvtran`. Full driver run (`xscrncsigma`),
`qp_energies.dat` compared to the sym=off baseline by matched `eps_HF`:

```
 eps_HF     sym=off eps_QP   sym=on eps_QP    diff
-20.242       -20.25952       -20.25950       2D-5
 -1.269        -1.319428       -1.319382       5D-5
 -0.618        -0.616500       -0.616503       3D-6
 -0.453        -0.472807       -0.472815       8D-6
 -0.391        -0.408053(*)    -0.408248(*)    2D-4  (*not converged, BOTH --
                                                        same pre-existing
                                                        Pade issue as
                                                        "continued 17",
                                                        unrelated to sym)
  0.606         0.603966        0.603965       1D-6
  0.743         0.737865        0.737859       5D-6
```

All converged MOs now agree to 1e-6 to 1e-5 (residual fully consistent with
`cc_conv=9,scf_conv=10` truncation, not a structural bug); the one
unconverged MO fails to converge THE SAME WAY in both sym=off and sym=on,
confirming that non-convergence is the already-flagged, separate Pade-
robustness question, not a symmetry artifact. **The sym=on bug is closed.**

**Standing guard added** (`scrnc_get2ints.F`, right after `Nirrep` is read):
under `Nirrep>1`, `Inquire`s for all 4 files and calls `Errex` with a clear
message naming the missing file and pointing at `SAVE_INTS=ON` if any are
absent, instead of silently proceeding with zero-filled data -- this is
exactly how the bug went undetected for so long (see "continued 19"'s own
finding: `Scrnc_get2ints` just `Dzero`s and reads whatever files exist,
no completeness check). Verified both directions: without `SAVE_INTS=ON`
(files already deleted from a prior run), `xscrnc` now exits 1 with
`@-Scrnc_get2ints: FATAL -- symmetry AO integral file missing under
Nirrep>1: IIJJ` printed before `Scrnc_get2ints` even reaches its main
loop; with `SAVE_INTS=ON`, exits 0 with no FATAL line.

**Practical consequence for future scrnc test ZMATs under sym=on**:
`SAVE_INTS=ON` must be included, or the new guard will (correctly) abort
the run rather than let it proceed with silently-wrong physics. Sym=off
runs are unaffected either way (the guard only fires for `Nirrep>1`, and
`IIII` alone has been observed to survive `xvtran` regardless of this flag
in every sym=off case checked -- not investigated further since sym=off
was never broken).

**Committed**: `scrnc_get2ints.F` (the guard). This log entry.
**Archived**: `sigma_debug_archive/h2o_symon_qp_energies_saveints.dat`
(final, correct sym=on `qp_energies.dat`, for direct comparison against
`h2o_qp_energies_postfix.dat` from "continued 17"). Pipeline directory
left on disk, not committed: `/blue/ufhpc/perera/scrnc_test/
h2o_scrnc_sigma_debug_symon_saveints/`.

**Next**: per Ajith's own framing at the start of this whole thread --
"sym on off bug needs to be fixed first, then look at Pade, and in Pade
what we should follow is math not what happens to the results" -- the
sym=on/off work is now done. The Pade-robustness question (MO 5/4's
non-convergence, present identically in both symmetry settings, first
flagged in "continued 13" and sharpened in "continued 17") is next, with
the explicit instruction to reason from the mathematics of the
continuation/QP-equation rather than tuning against what "looks right."

---

## SESSION 2026-08-18 (continued 19) -- Iao_pair permutation fix IMPLEMENTED + VERIFIED (real, order-of-magnitude improvement); traced the REMAINING residual to a second, deeper, ARCHITECTURAL gap: Scrnc_get2ints.F's raw AO-integral files don't contain complete data under sym=on -- NOT fixed, needs Ajith's call on the right replacement

Direction 1 from "continued 18" implemented: `Scrnc_jj_unpack_natural.F` (new
file) rebuilds the packed<->natural permutation from `Iaopop`/`Dirprd` alone
(mirrors `Scrnc_reord_cols.F`'s own construction) and un-permutes
`JJ_lrespn` in place; wired into `scrnc_driver.F` right after either CC
path builds `JJ_lrespn`, before `Scrnc_wc_ao2mo_test`/`Scrnc_sigma_prep`
touch it. No-op under sym=off (early `If (Nirrep.Eq.1) Return`).

**Verified working, not just plausible**: rebuilt, reran H2O/STO-3G sym=on
and sym=off. `Scrnc_wc_ao2mo_bfcheck` now gives the IDENTICAL residual
between the two runs (`6.217249D-15` both). `Wc_reduced(p,q,m)` comparison
(matched by `eps_HF` permutation, same method as "continued 17"): 194/343
entries (57%) now match sym=off EXACTLY (<1e-6), and the worst-case
discrepancy dropped from ~4.1 (unbounded/wild) to ~0.11 -- an order of
magnitude. **This part of the fix is solid and should be kept.**

**Residual investigated -- root cause is NOT a leftover ordering bug**,
confirmed several ways before concluding this:
- Both matrices still internally p<->q symmetric to ~1e-17 (an ordering
  bug would show up here; it doesn't).
- The worst offenders have a distinctive signature: `sym=on` comes out
  EXACTLY `0.000000` where `sym=off` has a real nonzero value, including
  fully-diagonal entries like `Wc_reduced(4,4,4)` (same MO three times --
  no symmetry reason to vanish).

**Traced with the same "print and find where it breaks" method, three
layers deeper**:
1. `Scrnc_jj_rcc_assemble.F`'s own "Response-only JJ matrix" dump: packed
   row 17 is entirely zero, every other row in the same 21-entry irrep-1
   product block is populated. Packed position 17 = natural AO-pair (5,5)
   (verified via the now-exact `Iordr_ao` permutation) -- the self-pair of
   AO#5, the SOLE AO in irrep-2.
2. Compared `_DEBUG_LVL0` checksums already in the log for Npert_count=16
   (works) vs 17 (broken), in `Scrnc_prep_culomb_ints.F`'s existing prints:
   `<oo|j|xx>` nonzero in BOTH; `<vv|j|xx>`/`<vo|j|xx>` nonzero for 16,
   EXACTLY ZERO for 17. This is where the zero first appears in that
   routine's own output.
3. Added a new diagnostic (`Scrnc_dump_dble_array`, new small helper in
   `scrnc_xxxx_pqxx.F` -- needed because that routine's own
   `Implicit Integer (A-Z)` would silently bit-reinterpret a direct dump
   of `Icore` as integers, the exact same trap as the already-fixed
   Jj_lrespn bug from "continued 8"; routing through a genuinely
   Double-Precision-typed callee sidesteps it) dumping the raw `I(XX,XX)`
   input block itself, before ANY transform. Result: packed row/col 17
   (natural (5,5)) has exactly ONE nonzero entry (itself, the self-Coulomb
   diagonal) out of 21; rows/cols 18 and 21 (natural (6,6) and (7,7), the
   two self-pairs of the two irrep-3 AOs) show the SAME pattern -- zero
   cross-terms with columns 1-16 (all pairs native to irrep-1's own 4 AOs),
   nonzero only among themselves (columns 17-21, the "self-pairs of
   non-irrep-1 AOs" block). Rows 1-16 have rich cross-structure among
   themselves and are correspondingly zero for columns 17-21.

**This block-isolated pattern is NOT explained by Abelian point-group
selection rules** -- for an Abelian group, any AO self-pair's density is
ALREADY totally symmetric regardless of which irrep the underlying AO
belongs to, so a bare Coulomb integral like (1,1|5,5) between two
totally-symmetric densities is allowed to be nonzero by symmetry; there is
no group-theoretic reason for it to vanish. And it isn't a coincidental
physical zero either -- checked `Jj_bare` (built via an entirely SEPARATE
call chain, `Scrnc_sigma_bare_ao.F` -> its own fresh `Scrnc_get2ints` call,
no packing/reordering involved at all) and found the IDENTICAL zero
pattern there too (`bare(1,1;5,5)=0.0`, `bare(5,5;6,6)=bare(5,5;6,7)=
bare(5,5;7,6)=bare(5,5;7,7)=0.0`, etc.) -- two independent code paths
agreeing on the same wrong answer is exactly this project's own
established "both wrong the same way" failure signature (see "continued
8"), not evidence of correctness.

**Root cause traced to `Scrnc_get2ints.F` itself, one further level down**
-- both `Jj_bare` and `List_xxxx` (hence `JJ_lrespn`) ultimately source
their raw AO two-electron integrals from this ONE shared routine. It reads
four legacy, hardcoded filenames (`IIII`,`IIJJ`,`IJIJ`,`IJKL` -- the
classic ACES2/vmol integral-class sort, split by which index positions
coincide), using `Inquire(...,Exist=...)` to skip any that aren't present,
then reconstructs the full dense `Nbfns**2 x Nbfns**2` matrix via simple
8-fold PERMUTATIONAL symmetry only (lines 91-98: fill all 8
index-swap-equivalent slots from one stored value) -- no point-group
symmetry unfolding logic anywhere in this routine. Checked directly: in
BOTH the sym=on and sym=off H2O/STO-3G test directories, only the `IIII`
file exists (`IIJJ`/`IJIJ`/`IJKL` absent in both) -- same file, nearly
identical size (79964 vs 79864 bytes) -- but different CONTENT. Working
hypothesis (not yet independently confirmed against vmol's own source):
under real symmetry, `xvmol`'s integral sort applies genuine point-group
compression when writing this file (skips integrals derivable from others
via a symmetry OPERATION, not just an index permutation), and
`Scrnc_get2ints.F`'s reconstruction -- built and only ever exercised under
sym=off, where point-group compression is trivial/absent -- has no code to
undo that specific kind of compression, silently leaving the
un-reconstructed slots at their `Call Dzero` initial value.

**Not fixed. Deliberately stopped here** rather than attempting a
work-around alone: this is an architectural gap in how scrnc acquires raw
AO integrals under symmetry, not a local index bug -- the same class of
call this project's standing practice reserves for Ajith (see "continued
17"'s framing: physics-affecting research-pipeline changes, verify-before-
trusting). The right fix is very likely NOT extending
`Scrnc_get2ints.F`'s own unfolding logic (would mean re-deriving/
reimplementing point-group integral-symmetry expansion by hand, exactly
the kind of "why is this project's own machinery duplicated/reinvented"
risk flagged elsewhere in this log) but replacing this whole raw-file
read with whatever symmetry-aware AO two-electron integral retrieval the
REST of this codebase already trusts (e.g. however `List_xxxx`'s own
GETLST-based construction gets ITS starting data, or however vcceh/other
established ACES2 modules read full AO integrals under symmetry) -- but
which existing mechanism is right, and whether it is even called at a
point in the pipeline where scrnc could reuse it, needs Ajith's read of
ACES2's own integral-list conventions, not more code archaeology alone.

**Diagnostics left in the tree** (all `_DEBUG_LVL0`-gated, harmless,
intentionally kept rather than reverted -- see file-by-file list so none
of this needs re-deriving on resume):
- `scrnc_jj_rcc_assemble.F`: `Isympert`/`Iaopop` dump (confirms packed
  ordering); plain dumps of `JJ_lrespn` before bare-v add
  (`jjlrespn_resp_packed.dat`) and after (before Symmet2)
  (`jjlrespn_plusbare_packed.dat`).
- `scrnc_sigma_prep.F`: plain dump of `Jj_bare` in natural order
  (`jjbare_natural.dat`).
- `scrnc_xxxx_pqxx.F`: new `Scrnc_dump_dble_array` helper; call site
  dumping the raw `I(XX,XX)` block for `Iflag=4,Irrepxx=1,Ispin=1`
  (`ixxxx_iflag4_irrepxx1.dat`).
- `Scrnc_prep_culomb_ints.F`/`Scrnc_lrspn_rcc_main.F`: pre-existing
  `_DEBUG_LVL0` checksums (not new this session), just newly load-bearing
  for this investigation -- no changes needed, already sufficient.

**Fixed and installed**: `scrnc_jj_unpack_natural.F` (new),
`scrnc_driver.F` (wired the call in). **Not fixed**: the deeper
`Scrnc_get2ints.F` gap -- H2O/C2v `Wc_reduced`/QP energies are IMPROVED
but still NOT fully correct under sym=on; do not treat sym=on results as
trustworthy yet.

**Committed**: all of the above (fix + diagnostics), this log entry.
Archived: `sigma_debug_archive/h2o_symon_isympert_diag.txt`,
`h2o_symoff_isympert_diag.txt` (from "continued 18"). Pipeline directories
(with the raw `.dat` diagnostic dumps) left on disk, not committed, per
project convention: `/blue/ufhpc/perera/scrnc_test/
h2o_scrnc_sigma_debug_sym{on,off}_oneterm_isolation/`.

---

## SESSION 2026-08-18 (continued 18) -- sym=on bug ROOT-CAUSE FOUND via debug printing: Jj_lrespn's Iao_pair index is irrep-grouped/packed, NOT the natural (Mu-1)*Nao+Nu order Scrnc_wc_ao2mo.F's header comment claims -- exact break point pinned down, fix not yet applied

Per Ajith's instruction: use debug-level printing to find the FIRST place
things break, rather than continuing to reason forward from symptoms. This
is the SAME bug as the 2026-08-16 "real Abelian symmetry" session's "real
cause still OPEN" A1-block finding (see that entry, much further down this
file) -- today's work traces it to an exact, concrete mechanism that session
didn't reach.

**Ruled out first** (bisecting backward from the symptom): `SCFEVALA`/
`SCFEVECA` pairing is self-consistent even under `sym=on` -- traced
`vtran/dropvc.F`'s write of both records (lines ~316-384): both go through
the SAME `Do I=1,NBAS` masking loop with the SAME running `ITHRU` index, so
eigenvalue slot i and eigenvector column i are always the same physical
orbital, just in native-SCF (irrep-blocked) order rather than energy order.
Not the bug.

**Root cause, confirmed by direct printed evidence**: `Scrnc_wc_ao2mo.F`'s
own header comment claims `Jj_lrespn`'s `Iao_pair` index is
`(Mu-1)*Nao+Nu` (natural row-major AO-pair order), citing
"`scrnc_prep_2ints.F:68-77`" as authority. That citation is wrong -- those
lines only build `Irrep_ao_pairs`, a per-natural-position LABEL array;
the actual `Iao_pair` numbering `Jj_lrespn` is built with comes from
`Scrnc_lrspn_rcc_main.F` (lines ~127-135): `Do Irrepr=1,Nirrep; Do
Iao_pair=1,Irpdpdao(Irrepr); Npert_count=Npert_count+1;
Isympert(Npert_count)=Irrepx` -- i.e. Iao_pair position = a running
counter that visits ALL of one irrep's AO-pairs first, then the next
irrep's, etc. -- genuinely irrep-GROUPED, not natural raster order.
`scrnc_prep_2ints.F`'s own `List_xxxx` (the bare-integral piece added into
`JJ_lrespn` in `Scrnc_jj_rcc_assemble.F`, lines 125-141) uses this same
irrep-grouped convention via `Scrnc_reord_cols.F`/`Scrnc_sympack_aoints.F`.

Added a temporary `_DEBUG_LVL0`-gated print in `Scrnc_jj_rcc_assemble.F`
(kept in the tree, gated, not a permanent unconditional print) dumping
`Iaopop(1..Nirrep)` and `Isympert(1..Nao_pair)`. H2O/STO-3G/C2v,
`Iaopop=(4,1,2,0)`: printed `Isympert` sequence is `21` consecutive `1`s,
then `8` consecutive `2`s, then `16`ish consecutive `3`s (one early `1`
mixed at position 21, then clean), then `4` consecutive `4`s -- a perfectly
clean grouping into 4 contiguous same-irrep runs
(archived: `sigma_debug_archive/h2o_symon_isympert_diag.txt`). Cross-
checked the block sizes by hand: direct-product-symmetry pair counts
`Sum_irrep Iaopop(irrep)**2` (totally-symmetric/self-product pairs) `=
4**2+1**2+2**2+0**2 = 21` -- matches the first run exactly. **This
pattern is only possible if Iao_pair positions are sorted/grouped by
irrep; natural `(Mu,Nu)` raster order would interleave different
product-irreps within nearly every row** (confirmed by hand for `Mu=1`:
as `Nu` ranges `1..7`, the product irrep changes every few steps
depending on `Nu`'s own irrep, never staying constant for 21 straight
positions) -- so this is conclusive, not a coincidence of AO-block
ordering. Sym=off control run (`Iaopop=(7)`, one irrep) printed all-`1`s,
uninformative by construction (only one irrep exists, so the question
"packed vs natural" can't be distinguished there) -- exactly why this bug
has been invisible on every `sym=off` run in this project's history,
including the original 2026-08-16 investigation's own controls.

**Consequence**: `Scrnc_wc_ao2mo.F` (AO2MO transform) and
`Scrnc_sigma_prep.F`'s `Jj_corr = Jj_lrespn - Jj_bare` subtraction (`Jj_bare`
built via `Scrnc_get2ints`, genuinely natural `(Mu-1)*Nao+Nu` order, per
`Scrnc_sigma_bare_ao.F`'s own explicit `Do Mu...Do Nu...Iao_pair=
(Mu-1)*Nao+Nu` loop) BOTH treat `Jj_lrespn` as natural-order when it is
actually irrep-grouped-order. Under `sym=off` these two orderings coincide
(trivial, one irrep = no grouping needed), so every check this project has
ever run against `sym=off` passed cleanly regardless. Under `sym=on` they
diverge for the large majority of AO-pair slots, so `Jj_corr` -- and
everything built from it (`Wc_mo`, `Wc_reduced`, `Sigma_c`, `eps_QP`) --
silently combines mismatched AO pairs. This fully explains "continued 17"'s
finding (134/343 `Wc_reduced` entries wrong by orders of magnitude after
correcting for MO-index permutation, localized by orbital identity rather
than uniform) and matches 2026-08-16's finding (asymmetry confined to a
specific irrep block) as two different symptoms of the same root
mismatch.

**Not yet fixed.** Two candidate fix directions, not yet chosen or
attempted:
1. Reorder `Jj_lrespn` from irrep-grouped order back to natural
   `(Mu-1)*Nao+Nu` order right after `Scrnc_jj_rcc_assemble.F` builds it
   (or inside it, before return) -- everything downstream
   (`Scrnc_wc_ao2mo.F`, `Scrnc_sigma_prep.F`) then needs no changes at all,
   since their natural-order assumption becomes true again. Needs the
   actual inverse permutation, which `Scrnc_reord_cols.F`'s own `Iord`
   array (built in `scrnc_prep_2ints.F`, currently used one-way,
   natural->packed) may already encode, or can be reconstructed directly
   from `Isympert`+the per-irrep AO index lists (`Irrep_aos`/`Iaopop`/
   `Ioffao`).
2. Build `Jj_bare` in `Scrnc_sigma_prep.F` in the SAME irrep-grouped order
   instead (so the subtraction is apples-to-apples), then make
   `Scrnc_wc_ao2mo.F`'s AO2MO transform aware of the same permutation when
   applying the MO coefficient matrix `C` (whose row index needs to be
   confirmed as natural-AO-order or not -- not yet checked). More
   invasive, touches more files, higher regression risk to the
   already-verified `sym=off` path.

Direction 1 looks preferable (minimal blast radius, sym=off path
untouched) but the actual permutation array needs to be nailed down and
verified (e.g. against the `sym=off` case, where it should reduce to the
identity) before writing it -- next step on resume.

**Committed**: `scrnc_jj_rcc_assemble.F` (the `_DEBUG_LVL0`-gated
diagnostic print, kept for future reference); this log entry; new archive
files `sigma_debug_archive/h2o_symon_isympert_diag.txt` /
`h2o_symoff_isympert_diag.txt`. Rebuilt `xscrnc` via `gmake`/`gmake
install`, md5-verified installed.

---

## SESSION 2026-08-18 (continued 17) -- negative-frequency folding fix IMPLEMENTED + independently pilot-code-verified; QP results change substantially (new non-convergence + sign flips, not yet assessed); NEW independent bug found via sym=on/sym=off consistency check (pre-existing, unrelated to the fix)

Picked up "continued 16"'s confirmed gap. Per [[feedback_pilot_code_debugging]]
(standalone minimal code as ground truth) and
[[feedback_consistency_check_debugging]] (exhaust checks from every angle),
Ajith's instruction this session: complete the fix, verify with a pilot code,
and use the sym=on/off and RHF/UHF consistency checks that have worked before.

**Fix implemented**: `scrnc_sigma_driver.F`, fold loop (~line 227) -- added
the second propagator term `Dim2 = Omega_ext(J) - Xnode(K)` alongside the
existing `Dim = Omega_ext(J) + Xnode(K)`, both sharing the same
`Wc_reduced(.,.,Im)` and `Wnode(K)` weight, exactly as derived in "continued
16". `Tre`/`Tim` now sum both terms before multiplying into
`Sigma_re`/`Sigma_im`.

**Pilot-code verification (independent of the production code)**: wrote
`pilot_sigma_fold_check.f90` (scratchpad, archived to
`aces2/scrnc/sigma_debug_archive/pilot_sigma_fold_check.f90` +
its output `pilot_sigma_fold_check_output.txt`) -- a standalone program with
NO connection to ACES2/chemistry. Synthetic surrogate `g(w')=1/(1+(w'/b)**2)`
stands in for `Wc_reduced` at one fixed `(p,q,m)`; ground truth is
`Int_{-inf}^{inf} dw' [1/(i(omega+w')-eps_m)]*g(w')` evaluated via a
substitution `w'=b*tan(theta)` (maps the tail exactly, removes it
analytically) + fine composite Simpson quadrature -- a genuinely different
substitution and quadrature rule than production's Gauss-Legendre + rational
map, so this isn't the code checking itself. Reused the actual
`Scrnc_sigma_freq_grid`/`Scrnc_sigma_gaussl` routines verbatim (same
duplicated-per-module convention) to build the real production grid, then
compared the ONE-TERM (pre-fix) and TWO-TERM (post-fix) folded sums against
the Simpson ground truth, at 9 `(eps_m,omega)` combinations spanning
`eps_m={-0.5,0.3,-20.6}` (the last matching the DZ core-orbital scale from
"continued 13") and `omega={0,0.3,0.75}`.

Result: **two-term matches ground truth to `1e-9`-`1e-15` by N=32, at every
single test point**; **one-term is wrong by 5-90% at every test point and
does NOT converge to the true value as N grows** (converges cleanly to a
stable but incorrect number -- the signature of a missing term, not a
quadrature-precision issue). Hit one self-inflicted snag worth recording:
first run had the two `folded_sum` calls' `one_term` boolean swapped,
making the mislabeled output look like the OPPOSITE conclusion (one-term
matching, two-term wrong) -- caught by a second diagnostic program
(`diag_halfvsfull.f90`, not archived, pure scratch) that split the exact
integral into its `(0,inf)` and `(-inf,0)` halves explicitly and noticed the
mislabeled "two-term" column exactly equalled the `(0,inf)`-only half-integral,
not the full one -- fixed the call, reran, got the result above. Worth
remembering: a pilot code's OWN wiring needs the same "don't trust it blindly"
scrutiny as the thing it's checking.

**Rebuilt and reinstalled**: `module load intel/2025.1.0`; `ifx -i8 -o
scrnc_sigma_driver scrnc_sigma_driver.F`; installed to `aces2/bin/xscrncsigma`;
md5-verified installed binary matches the freshly-linked local one (per
[[feedback_aces_build_install_trap]]).

**Reran H2O/STO-3G (sym=off, REF=UHF) DRCCD and CCSD** in the existing
`h2o_scrnc_sigma_debug`/`h2o_scrnc_sigma_debug_ccsd` dirs (same
`sigma_input`/`ZMAT` as "continued 14"). Both exit 0. Self-checks unaffected
by the fix, as expected (fix is isolated to fold arithmetic downstream of
them): `Scrnc_sigma_bare_ao` still exactly `0.0`; `Scrnc_wc_ao2mo_bfcheck`
still `~1D-14`-`3D-14`.

**But the QP results themselves changed substantially** (post-fix
`qp_energies.dat` archived as `h2o_qp_energies_postfix.dat` /
`h2o_ccsd_qp_energies_postfix.dat`):
- DRCCD: MO 5 (`eps_HF=-0.391`) now FAILS to converge (51 iters, hits cap) --
  converged cleanly pre-fix. `Im[Sigma_c]` is now NEGATIVE for MOs 3,4,6,7
  (was positive for all 7 MOs pre-fix) -- per "continued 13"'s own flagged
  convention this was treated as a symptom of a bad Pade
  extrapolation/spurious pole, not dismissed as fine there, so it shouldn't
  be dismissed here either.
- CCSD: same pattern (MO 5 fails to converge, same sign flips), and MO 5's
  correction is enormous (`eps_QP=-0.629` vs `eps_HF=-0.391`, a ~0.24 Hartree
  shift) -- given non-convergence this number is not trustworthy as-is.
- Recomputing "continued 14"'s CC-vs-dRPA diff table with the fix (excluding
  the non-converged MO 5) gives differences up to ~1.8 eV (MO 2), several
  times larger than the pre-fix ~0.33 eV max -- a real, qualitative change to
  the internal CC-vs-RPA comparison this project is building toward as
  manuscript evidence.
- **Not yet assessed**: whether this is "the fix correctly revealing a
  bigger effect" or "the fix destabilized the existing Pade continuation
  in a new way" (or both) -- needs the same Pade-robustness discussion
  already flagged open in "continued 13"/"continued 15", now sharper since
  it affects more MOs than before. Deliberately not chased further or
  patched this session -- a numerics/physics call for Ajith, same standing
  practice as always.

**NEW finding, independent of the fold fix -- sym=on vs sym=off consistency
check (Ajith's suggested check) surfaces a real, separate bug.** Built
`h2o_scrnc_sigma_debug_symon` (identical ZMAT to the sym=off baseline,
only `sym=on`) and reran. `eps_HF` values matched the sym=off run exactly
(just reordered/re-blocked by irrep, as expected), but the QP corrections
for the SAME physical MOs (matched by `eps_HF` value) differ wildly --
not a small numerical difference: e.g. MO with `eps_HF=-20.24174` gets
`corr=-17.8 mHa` under sym=off vs `corr=-835 mHa` under sym=on; MO with
`eps_HF=-0.618` gets `+1.5 mHa` vs `+250 mHa`; sign flips on two other MOs
too. Magnitudes up to ~0.8 Hartree apart -- far outside anything explicable
as "different code path, same physics, small numerical difference."

**Isolated from the fold fix**: rebuilt a throwaway one-term (pre-fix
formula) control binary (`scrnc_sigma_driver_oneterm`, scratchpad only, not
committed/archived -- built by reverting the fold-loop patch on a scratch
copy of the file, real repo file untouched throughout, md5-verified after)
and reran BOTH sym=off and sym=on with it. **The sym=on/sym=off discrepancy
is present with the OLD one-term formula too**, same order of magnitude
(archived as `h2o_symoff_oneterm_qp_energies.dat` /
`h2o_symon_oneterm_qp_energies.dat`, plus the two-term sym=on result as
`h2o_symon_twoterm_qp_energies.dat`). **Conclusion: this is a pre-existing,
independent bug, not caused by or related to today's fold-formula fix.**
Leading suspicion (not confirmed, not investigated further this session):
a MO-index-ordering mismatch between `Wc_reduced(p,q,m)` and `Eps(m)` as
written to `wcreduc.dat` by `Scrnc_sigma_prep.F` -- irrep-block MO ordering
under `sym=on` vs plain energy ordering under `sym=off` -- but this needs
its own dedicated investigation, not guessed at further here.

**REF=RHF consistency check attempted, inconclusive**: `h2o_scrnc_sigma_debug_rhf`
(same ZMAT, `REF=RHF` instead of `UHF`, `sym=off`) -- the DRCCD/CC step itself
completes fine (`Rcl_check_uhf` T-amplitude checksums print normally), but
the pipeline dies with `@ACES_EXIT` right after, before `xscrnc` finishes,
at node 1. Never previously exercised -- every prior scrnc test in this whole
project uses `REF=UHF` even for closed-shell H2O (established convention,
not investigated as to why). Not chased further -- appears to be a
pre-existing limitation of this code path with `REF=RHF`, unrelated to the
fold fix, orthogonal to both bugs above. Left open.

**Committed**: `scrnc_sigma_driver.F` (the fold fix only); this log entry;
new archive files under `aces2/scrnc/sigma_debug_archive/` (post-fix
`qp_energies.dat` for H2O DRCCD/CCSD, the sym=on/sym=off one-term/two-term
comparison outputs, the pilot code + its output). Full pipeline directories
left on disk, not committed, per project convention:
`/blue/ufhpc/perera/scrnc_test/h2o_scrnc_sigma_debug{,_ccsd,_symon,_rhf,
_symoff_oneterm_isolation,_symon_oneterm_isolation}/`.

**Next steps, both needing Ajith's call, neither started**: (1) Pade
robustness given the fix now destabilizes MO 5 and flips several
`Im[Sigma_c]` signs -- same open numerics question as "continued 13"/
"continued 15", now more pressing; (2) investigate the sym=on MO-index-
ordering bug in `Scrnc_sigma_prep.F` (or wherever the actual mismatch turns
out to be) -- brand new, not previously known, found only because the
sym=on/off consistency check was run as a check on the fold fix and
happened to catch something unrelated.

---

## SESSION 2026-08-17 (continued 16) -- verified item 1 from "continued 15": the omega'-integral folding IS missing a term (CONFIRMED bug); grid-conflation worry (item 3/4) is a false alarm; Ajith is reconsidering whether to frame this as GW at all

**Housekeeping note first, corrected**: Ajith clarified same day: "Similarity
to GW is undisputed. One improvement or difference is linear response. We
want something new. That is all I am saying." He is NOT rejecting the GW
framing (an earlier version of this note overread his first remark as
that -- wrong, corrected here). The GW resemblance stays; CC-linear-response
W is already one acknowledged difference from RPA-GW, but he doesn't yet
consider that sufficient novelty on its own -- still studying GW more,
still deciding what the additional angle is. So the Pade/contour-
deformation/GWGamma-hierarchy material in "continued 15" is NOT moot, still
live. The frequency-reduction check below is unaffected either way -- it's
just whether the existing quadrature correctly represents the omega'
convolution it claims to.

**Further refinement same day**: "Nothing we did is moot. They
improve/debug our code. If we get big improvement with CC, then CCLR-GW is
worth a manuscript." All debugging/correctness work stands on its own
value regardless of the novelty question. And the novelty bar itself is
sharper than "need one more new ingredient beyond CC-response screening":
if CC screening gives a BIG improvement (not just a nonzero effect) over
RPA-GW, that's sufficient on its own -- CCLR-GW would be the manuscript.
The concrete next test this points to is benchmarking CC-vs-dRPA screened
QP energies against a trusted reference (CCSD(T) IPs or GW100-style set),
not searching for a separate new angle -- see
[[project_scrnc_gw_manuscript_goal]] for full framing.

**Item 3/4 from "continued 15" (quadrature-vs-continuation grid conflation)
-- FALSE ALARM, already correctly separated in code.** Checked
`scrnc_sigma_driver.F` directly: `Xnode(K)`/`Wnode(K)` (the omega' quadrature,
built by `Scrnc_sigma_freq_grid`, N nodes mapped from Gauss-Legendre(-1,1) to
(0,infinity) via the standard `Omega_c*(1+x)/(1-x)` substitution) and
`Omega_ext(J)` (the external omega_j points, read directly from
`sigma_input`, completely independent count `Nomega`) are already two
distinct arrays, and the main fold loop (`Do J.../Do Im.../Do K` structure,
lines ~227-244) already reuses one pipeline run's `Wc_reduced` across every
J with no extra pipeline cost -- i.e. the "many cheap external omega_j
evaluations from one expensive omega' quadrature" idea the external review
recommended is already exactly how the code works. Nothing to fix here.

**Item 1 (negative-frequency folding) -- CONFIRMED, real gap, not yet
fixed.** Re-derived the reduction: starting from the full-axis convolution
`Sigma_c(iw) = -(1/2pi) Int_{-inf}^{inf} dw' G0(iw+iw') Wc(iw')` and using
`Wc(-iw')=Wc(iw')` (Wc is an even function of the response frequency, true
for any real-Hamiltonian response function), substituting w'->-w' on the
negative half and folding gives

  Sigma_c(iw) = -(1/2pi) Int_0^inf dw' Sum_m
     [ 1/(iw+iw'-eps_m) + 1/(iw-iw'-eps_m) ] Wc(pm,qm;iw')

i.e. TWO propagator terms per quadrature node, not one. Read the actual
fold loop in `scrnc_sigma_driver.F` (lines ~227-244): for each node K and
external point J it computes `Dim = Omega_ext(J) + Xnode(K)`,
`Dre = -Eps(Im)`, then `T = 1/(Dre+i*Dim)` -- exactly
`1/(i*omega_j+i*omega'_K-eps_m)`, accumulated as
`Sigma += (-Wnode(K)/2pi) * T * Wc_reduced`. **There is no companion term
using `Omega_ext(J) - Xnode(K)`** -- the second propagator from the
derivation above is simply absent from the sum. No occ/virt iη branch or
other mechanism was found elsewhere in the file that could be silently
supplying it. This is a genuine, previously-undetected gap in the Sigma_c
construction itself (independent of the Pade/continuation questions), not
a documentation/notation issue as first suspected.

**Not yet fixed** -- deliberately left for Ajith's return given (a) it's a
physics-affecting change to the research pipeline (standing "verify before
trusting/committing" practice) and (b) the open GW-vs-new-direction question
above may change whether this code path matters at all. If/when resuming:
the fix is mechanical -- inside the `Do J / Do Im` block, add a second
`Dim2 = Omega_ext(J) - Xnode(K)` term with the same `Eps(Im)`/`Wnode(K)`
weight and accumulate it into `Sigma_re`/`Sigma_im` alongside the existing
term. Worth re-running the H2O/STO-3G DRCCD-vs-CCSD comparison
("continued 14") after the fix to see whether it changes the magnitude/sign
of the already-observed differences, since that comparison was run before
this gap was found.

---

## SESSION 2026-08-17 (continued 15) -- external theoretical review of the Sigma_c/Pade implementation + modern-GW positioning discussion, archived not yet acted on

Ajith pasted a discussion document (`~/claude_docs/CC_screened_GW_discussion.md`,
an external review of `SCRNC_SIGMA_C_EQUATIONS.md` and
`scrnc_sigma_prep.F`/`scrnc_sigma_driver.F`) to be archived. **Not yet
implemented or verified against code -- this is a set of recommendations to
evaluate, not a confirmed diagnosis.** Full text in that file; key points:

**Part 1 -- concrete implementation review, ordered "do not change AO->MO/Wc
machinery yet, proceed in this order":**
1. *Check first*: whether the working Sigma_c integral (0->infinity in
   omega') already implicitly contains the negative-frequency contribution
   that a symmetric W_c(-iw')=W_c(iw') reduction from [-inf,inf] should
   produce (would normally give two propagator denominators, not one) --
   flagged as possibly more fundamental than the Pade convergence problem.
2. The current pipeline is imaginary-axis evaluation + Pade continuation, not
   true contour deformation (which would evaluate Sigma_c(omega_real)
   directly via Green's-function residues, no continuation step).
3. **Likely actionable fix, low cost**: the omega' quadrature (expensive, CC-
   response) and the external omega_j evaluation grid (cheap, algebraic) are
   two independent grids. Suspects the `0.075i, 0.493i, 2.030i, 13.40i` nodes
   in the notes are the N=4 omega' quadrature, possibly conflated with the
   omega_j continuation points in code/notation -- worth auditing
   `scrnc_sigma_driver.F` for this. If they're entangled, decoupling them
   means the core-orbital Pade domain-of-validity gap (DZ, 2/14 MOs, "continued
   13") could be fixed by cheaply sampling many more external omega_j points
   over the relevant energy scale (e.g. log/rational grid out to ~20 Ha for
   the eps_HF~-20.6 core case) WITHOUT rerunning the expensive omega'
   quadrature -- reassess Pade only after this decoupling is confirmed/fixed.
4. Separately, suggests replacing the bare fixed-point QP iteration
   omega_{n+1}=eps_HF+Re(Sigma_c(omega_n)) with a bracketed root find on
   F_p(omega)=omega-eps_HF-Re(Sigma_c(omega)), to separate "is the
   continuation trustworthy" from "does the QP equation have a stable root" --
   endorses NOT adding damping just to force convergence (matches Ajith's
   existing "not silently patched" stance from "continued 13").
5. If the decoupled-grid approach still leaves Pade unreliable, implementing
   true contour deformation w/ residues becomes attractive precisely because
   the expensive ingredient (W(iw')) is already computed and reusable.
6. Endorses the existing HF-exchange argument (no separate Sigma_x module
   needed since F=h+J-Sigma_x already for canonical HF orbitals) as worth
   keeping in the writeup, not accidental.

**Part 2 -- modern-GW positioning for the novelty question**
([[project_scrnc_gw_manuscript_goal]]): frames scrnc as CC-improved
*screening* only (chi_CC -> W_CC -> Sigma=iG0*W_CC), explicitly NOT yet a
vertex-corrected theory (Sigma=iGWGamma) -- proposes the hierarchy
`G0W_RPA -> G0W_CC -> G0W_CC*Gamma_CC` as the natural framing, with the
2026-08-17 DRCCD-vs-CCSD comparison ("continued 14") being the first data
point on step 1->2 of that hierarchy. Names recent vertex-corrected/modern-GW
directions as context (G3W2, SOSEX, self-consistent GWGamma, BSE vertex
corrections, contour-deformation/imaginary-time alternatives to Pade,
low-scaling/THC methods) -- flags that vertex effects in the self-energy and
polarizability can partially cancel, so improving W alone need not improve QP
energies monotonically, i.e. the "real differences" found in "continued 14"
demonstrate an effect but not yet that it's an improvement.

**Not done**: no code changes, no verification of the [-inf,inf] reduction
claim, no audit of whether omega'/omega_j are actually conflated in
`scrnc_sigma_driver.F`. Natural next step per this doc's own recommended
order: verify the frequency-reduction identity first, then check the
quadrature-vs-continuation grid separation in code, before touching Pade
itself. Needs Ajith's read/call before acting, per the standing "verify
before trusting" and "positioning is his call" practice in this project.

---

## SESSION 2026-08-17 (continued 14) -- CCSD-path bug found+fixed (directory not cleared between nodes); first internal CC-vs-RPA screening comparison run -- REAL, non-trivial differences found

Ajith's "check internally" instruction, in response to "do we have
something novel to add": since ring-CC (RCCD/DRCCD) IS RPA, reformulated
through the same CC-response machinery ([[project_rpa_module_tests]]:
"RPA=ringcc/dRPA=dringcc confirmed in code"), the cleanest possible test
of whether going to genuine CCSD-level screening changes anything is
entirely internal to scrnc -- run the IDENTICAL sigma-driver/AO2MO/
Padé/QP pipeline on the SAME geometry/basis, changing only
`CALC=DRCCD` vs `CALC=CCSD` in the ZMAT. No new code needed for the
comparison itself; this was the first time the CCSD path was ever run
through the sigma-driver chain, though.

**Bug found and fixed first**: CCSD path crashed on the driver's
SECOND node (`@GETLST: Assertion failed. List (1,127) does not
exist.`), never on the first. Root cause: `scrnc_sigma_driver.F`
re-runs `xjoda...xscrnc` fresh each node but never actually DELETED
the previous node's leftover generated files first -- `xjoda`
apparently doesn't fully reset every record-container file (e.g.
`MOABCD`) that the CCSD linear-response path writes dynamic list
records into between separate `xscrnc` invocations in the same
directory. Same root cause as the already-documented "rerunning a
pipeline step on top of stale JOBARC is NOT idempotent" lesson, just
manifesting on a path never previously driven repeatedly in one
directory. Never surfaced on the ring-CC/DRCCD path (different list-
record usage pattern, apparently not colliding the same way).

**Fix**: added an explicit `find . -maxdepth 1 -type f ! -name ZMAT
! -name GENBAS ! -name sigma_input ! -name pipeline_cmds ! -name
'*.log' -delete` at the top of each node's loop iteration in
`scrnc_sigma_driver.F`, before writing `frequency`/running
`pipeline_cmds`. Still the SAME single shared directory (Ajith's
design choice, "continued 10") -- only the generated-file CONTENTS
get cleared each iteration, never the directory itself. **Verified
no regression**: reran the H2O/STO-3G DRCCD case (which never had
this bug) with the fix in place -- bit-identical `qp_energies.dat` to
before, digit for digit.

**CCSD path itself, once fixed**: clean on H2O/STO-3G. Exit 0,
bare-v symmetry exactly `0.0`, AO2MO bfcheck `3.6D-15` to `3.0D-14` at
all 4 nodes, all 7 MOs converged (unlike the DZ case -- STO-3G's
similar-magnitude core orbital, eps_HF~-20.2, converges fine on both
paths; only DZ's core, eps_HF~-20.6, hit the Padé domain-of-validity
issue from "continued 13" -- consistent, not a new puzzle).

**The comparison** (H2O/STO-3G, same geometry/ZMAT settings, same
sigma_input `N=4,Omega_c=1.0,Nomega=6`, only `CALC=DRCCD` vs
`CALC=CCSD` differs):

```
 p     eps_HF   dRPA corr   CCSD corr  diff(mHa)  diff(eV)
 1  -20.24174     -26.35m     -32.66m      -6.31    -0.172
 2   -1.26854       9.58m      19.91m      10.33     0.281
 3   -0.61801       3.65m       7.56m       3.91     0.106
 4   -0.45307       1.15m      -0.50m      -1.65    -0.045
 5   -0.39127      -2.52m      -7.68m      -5.16    -0.140
 6    0.60591      -3.63m      -7.78m      -4.14    -0.113
 7    0.74258      -9.88m     -21.96m     -12.08    -0.329
```

("corr" = eps_QP-eps_HF, the QP correction itself; diff = CCSD corr
minus dRPA corr.) **Real, non-trivial differences** -- up to 0.33 eV
(MO 7), several times the typical ~0.1-0.2 eV "good agreement" bar
used in GW100-style benchmark comparisons, and MO 4 even flips SIGN
between the two screening levels (dRPA: opens the gap slightly; CCSD:
closes it slightly). Going from RPA-equivalent to genuine CCSD-level
dynamical screening is NOT a small perturbative correction here, at
least for this toy system.

**What this does and doesn't establish, stated plainly**: this is
strong evidence the CC-response screening construction captures
physics that plain RPA misses and that it propagates through to
chemically/spectroscopically meaningful QP shifts -- necessary
groundwork for a novelty claim. It does NOT yet establish which
screening level is more ACCURATE (no comparison against CCSD(T) IPs,
experiment, or GW100-style references attempted) -- that comparison,
plus resolving the Padé robustness gap so the comparison can be
trusted beyond STO-3G-scale toy systems, is the natural next step
before treating this as manuscript evidence.

**Committed**: `scrnc_sigma_driver.F` (directory-cleanup fix); this
log entry; CCSD debug artifacts (small text outputs) under
`aces2/scrnc/sigma_debug_archive/h2o_ccsd_*`. Full pipeline directory
left on disk, not committed:
`/blue/ufhpc/perera/scrnc_test/h2o_scrnc_sigma_debug_ccsd/`.

---

## SESSION 2026-08-17 (continued 13) -- larger-basis (H2O/DZ) check: core-machinery clean, but surfaces a REAL Pade-robustness gap (2 of 14 MOs fail to converge) -- NOT silently patched, flagged for discussion

Per Ajith's "check everything with large basis" instruction. Built a
NEW H2O/DZ ZMAT (`h2o_scrnc_sigma_debug_dz`) rather than reusing the
pre-existing `h2o_scrnc_drccd_uhf_dz` test dir -- that one uses
`sym=on` and no `CC_EXTRAP=OFF`, settings never exercised through this
session's sigma-driver code path, so reusing it would have confounded
"does the basis-size scale-up work" with "does a different symmetry/
extrapolation setting work". New ZMAT matches the already-verified
H2O/STO-3G run's flags exactly (`sym=off, CC_EXTRAP=OFF,
DENSITY=RESPONSE, REF=UHF`), only `BASIS=DZ` changed -- confirmed
`H:DZ`/`O:DZ` present in GENBAS. Nao=14 (`Nao_pairs=196`, confirmed via
pipeline's own "number of ao pairs" printout), Nmo=14 -- 2x the AO/MO
count of every prior test.

**Core machinery: clean, no regressions.** Same `N=4, Omega_c=1.0,
Nomega=6` settings as the H2O/STO-3G run. Exit 0.
`Scrnc_sigma_bare_ao` symmetry check exactly `0.0` at all 4 nodes;
`Scrnc_wc_ao2mo_bfcheck` regression check `3.6D-14` to `6.5D-14` at all
4 nodes (larger than the 7-MO case's `~1D-14` but still comfortably
machine-precision-clean, expected growth with matrix size) -- the AO2MO
transform, bare-v subtraction, and frequency-orchestration machinery
all scale to 2x the basis with no code changes and no correctness
regression.

**QP solve: 12 of 14 MOs converged cleanly** (4-12 iterations,
well-behaved corrections of a few hundredths of a Hartree, same
qualitative pattern as every prior run). **2 of 14 (p=1, p=7) hit the
50-iteration cap without converging** -- investigated rather than
ignored or silently patched:
- `p=1` (O 1s core, `eps_HF=-20.559`): `Sigma_c(1,1;i*omega)` on the
  fitted `omega in [0,1.5]` grid is small and nearly flat
  (`Re` ranges only `-9.396D-3` to `-9.281D-3` across all 6 points) --
  evaluating the resulting Pade rational approximant at the on-shell
  point `z=eps_HF=-20.559`, roughly 14x past the edge of the fitted
  domain, is a genuine Pade EXTRAPOLATION-DOMAIN VIOLATION, not a bug
  in the continuation algorithm itself (already validated exact to
  `~1D-15` against a known analytic function in "continued 11", but
  only tested there for evaluation points near/inside the fitted
  domain). Result came out wildly different in scale from every other
  MO (`eps_QP=-19.58`, a ~1 Hartree shift vs `~0.01-0.06` Hartree for
  every valence/virtual MO) with a spurious NEGATIVE `Im[Sigma_c]`
  (unphysical for occupied-state removal energies in this convention --
  every converged MO has positive `Im[Sigma_c]`) -- both symptoms
  consistent with a spurious near-pole of the fitted rational function
  lying between the sampled domain and `z=-20.559`, not a real physical
  result.
- `p=7` (`eps_HF=0.311`, WELL INSIDE the fitted `[0,1.5]` domain --
  not an extrapolation issue): `Sigma_c(7,7;i*omega)`'s real part
  drops by ~8x between the first two sample points
  (`1.443D-3` at `omega=0` -> `1.832D-4` at `omega=0.3`) then stays
  flat (`~1-2D-4`) for the rest -- this sharp early kink in only 6
  sampled points is exactly the kind of feature a low-order Thiele
  fit can represent with a spurious pole landing close to the real
  axis near `z~0.3`, destabilizing the fixed-point iteration even
  in-domain. Also produced a spurious negative `Im[Sigma_c]`.

**Deliberately NOT patched with iteration damping**: damping (mixing
the fixed-point update) would very likely make both cases LOOK
converged, but would not fix an underlying bad Pade extrapolation --
it would just converge smoothly to a still-untrustworthy number. The
real fix needs a numerical-methods decision (wider/better-scaled
`Omega_ext` sampling per orbital-energy scale, higher Pade order more
points, an explicit domain-of-validity check that flags/excludes
QP corrections requested far outside the fitted imaginary-axis range,
or some combination) -- exactly the kind of call to raise with Ajith
rather than deciding alone. Left open, written up in the new equations
doc (`SCRNC_SIGMA_C_EQUATIONS.md`, this session) for discussion.

**Committed**: this log entry; DZ debug artifacts (small text outputs
only) under `aces2/scrnc/sigma_debug_archive/h2o_dz_*`. Full pipeline
directory left on disk, not committed:
`/blue/ufhpc/perera/scrnc_test/h2o_scrnc_sigma_debug_dz/`.

---

## SESSION 2026-08-17 (continued 12) -- full sigma-driver chain (orchestration + W_c reduction + Pade + QP solve) verified on H2O/STO-3G/DRCCD, not just H2

Ran the exact same `xscrncsigma` chain ("continued 10"/"continued 11")
unmodified on a second, larger test system -- same incremental-
verification-first pattern as every other piece of this project (H2
first, then H2O), per "continued 11"'s own flagged open item.

Setup: fresh `/blue/ufhpc/perera/scrnc_test/h2o_scrnc_sigma_debug/`,
`ZMAT`/`GENBAS` seeded from the already-validated
`h2o_scrnc_drccd_uhf_sto3g` (7 AOs/7 MOs, STO-3G), same
`sigma_input`/`pipeline_cmds` as the H2 `Nomega=6` run (`N=4,
Omega_c=1.0`).

**Result**: exit 0. A `file` check on the captured log showed it as
"data" (binary) rather than "text" -- turned out to be harmless
(some subprocess mixes non-text bytes into stdout somewhere in the
pipeline; `grep -a` recovers plain-text matching fine) but is worth
remembering next time a `grep` on one of these logs comes back
suspiciously empty.
- `Scrnc_sigma_bare_ao` symmetry self-check: exactly `0.0` at all 4
  nodes (same as H2).
- `Scrnc_wc_ao2mo_bfcheck` (pre-existing AO2MO regression check):
  `4.4D-15` to `1.1D-14` at all 4 nodes -- clean, no regression, and
  properly exercised on a non-trivial 7x7 MO case (H2's own 2x2 case
  is nearly too small to stress-test the block-Xgemm transform
  machinery).
- `qp_energies.dat`, all 7 MOs converged (4-12 iterations, well under
  the 50-iteration cap):
  ```
    p      eps_HF          eps_QP       Im[Sigma_c(eps_QP)]   iters  converged
    1   -2.024174D+01   -2.026809D+01    1.425449D-02     4       1
    2   -1.268536D+00   -1.258960D+00    2.173615D-02     5       1
    3   -6.180077D-01   -6.143544D-01    1.387520D-02     5       1
    4   -4.530670D-01   -4.519204D-01    4.198156D-02    10       1
    5   -3.912734D-01   -3.937945D-01    5.288300D-02    12       1
    6    6.059097D-01    6.022754D-01    1.504180D-02     5       1
    7    7.425789D-01    7.326967D-01    1.772998D-02     6       1
  ```
  Physically sensible: MO 1 (O 1s core) barely perturbed by valence
  correlation (small, expected); valence MOs 2-5 shift by a few
  hundredths of a Hartree; virtuals 6-7 shift down slightly (narrowing
  gap, same qualitative direction as H2); all `Im[Sigma_c]` positive
  and reasonable-magnitude on-shell broadenings, no wild outliers or
  non-convergence.

**Verdict**: the full chain (external per-point orchestration -> bare-v
subtraction -> MO reduction -> Sigma_c accumulation -> Pade continuation
-> QP solve) is no longer H2-only -- confirmed working end-to-end on a
second, structurally different (lower symmetry, more MOs, core+valence
orbital spread) system with no code changes needed.

**Not yet done**: same open items as "continued 11" (degenerate-orbital
treatment, off-diagonal Sigma_c unused, Pade robustness under sparser/
noisier data not stress-tested) -- none newly opened by this H2O run.

**Committed**: this log entry; H2O debug artifacts (small text outputs
only, same convention as the H2 archive) under
`aces2/scrnc/sigma_debug_archive/h2o_*`. Full pipeline directory left
on disk, not committed:
`/blue/ufhpc/perera/scrnc_test/h2o_scrnc_sigma_debug/`.

---

## SESSION 2026-08-17 (continued 11) -- Pade continuation (imaginary->real axis) + on-shell QP equation solve built and verified end-to-end

Continuation of "continued 10" in the same session. Built the deferred
piece flagged there: real-axis Pade continuation of `Sigma_c(p,p;i*omega)`
and the quasiparticle-equation solve.

**Physics simplification applied (flagged, not silently assumed --
standard textbook result, not new physics)**: since this project's
reference is HF (not DFT), the usual "Sigma_x - v_xc" double-counting
correction in the QP equation is EXACTLY ZERO for canonical HF orbitals
-- eps_HF(p) already equals `<p|h+J-Sigma_x|p>` exactly (that's the
definition of "canonical": the Fock operator h+J-Sigma_x is diagonal in
this basis), and the GW self-energy's static-exchange piece is the SAME
Sigma_x operator. So this driver builds NO separate Sigma_x machinery at
all -- the on-shell QP equation used is simply
`eps_QP(p) = eps_HF(p) + Re[Sigma_c(p,p;eps_QP(p))]`, solved by plain
fixed-point iteration from `eps_QP^(0)=eps_HF(p)`.

**New code** (appended to `scrnc_sigma_driver.F`, same standalone
non-ACES2-linked program as "continued 10"):
- `Scrnc_sigma_pade_coeff`/`Scrnc_sigma_pade_eval`: Thiele reciprocal-
  difference continued-fraction Pade continuation (Vidberg-Serene 1977
  algorithm, standard numerical-GW technique) -- coefficient recursion
  and the standard forward two-term (`A(n)`,`B(n)`) evaluation
  recurrence.
- `Scrnc_sigma_qp_solve`: for each diagonal MO index p, builds the
  Pade fit from the already-computed `Sigma_c(p,p;i*Omega_ext(j))`
  points (pure post-processing on data the "continued 10" architecture
  already produces -- no new pipeline runs), then fixed-point-iterates
  the QP equation (50-iteration cap, `1D-8` convergence tolerance).
  Writes `qp_energies.dat` (`eps_HF`, `eps_QP`, `Im[Sigma_c(eps_QP)]`
  as the QP lifetime/broadening, iteration count, converged flag).
  Only runs if `Nomega>=2` (skipped, same as before, for a single-point
  debug run).

**Verified**:
- Pade routines validated FIRST, in isolation, against a known analytic
  function (`f(z)=1/(z-2)`, sampled at 8 imaginary points, evaluated
  back at 4 real points) before trusting them on real physics data --
  matched the exact closed form to `~1D-15` at every test point
  (`/tmp/.../scratchpad/pade_test.F`, not committed, throwaway
  verification harness).
- Full pipeline rerun on the same `h2_scrnc_sigma_debug` directory,
  `sigma_input` extended to `Nomega=6` (`Omega_ext=0.0,0.3,...,1.5`,
  same `N=4,Omega_c=1.0` inner quadrature as before -- confirms the
  "only N pipeline runs regardless of Nomega" design claim from
  "continued 10": still exactly 4 pipeline subprocess calls, now
  producing 6x the external-frequency data for free). Exit 0.
  `qp_energies.dat`:
  ```
    p      eps_HF          eps_QP       Im[Sigma_c(eps_QP)]   iters  converged
    1   -5.777270D-01   -5.720521D-01    2.758885D-03     4       1
    2    6.690811D-01    6.634062D-01    2.758829D-03     4       1
  ```
  Both MOs converged in 4 iterations. Physically sensible: small
  correlation shifts (`~6D-3` Hartree) in the expected direction
  (occupied MO less bound, virtual MO lower -- HOMO-LUMO gap narrows
  slightly under correlation, a commonly-seen GW correction direction),
  reasonable-magnitude on-shell lifetime broadening, and near-identical
  broadening for both MOs -- consistent with H2's symmetric bonding/
  antibonding MO structure.

**Not yet done**: multi-reference/degenerate-orbital QP treatment (not
needed for H2's non-degenerate case, may matter for larger/higher-
symmetry test systems later); off-diagonal `Sigma_c` is computed but
not used in the QP solve (standard for a non-degenerate on-shell
treatment); no test yet on H2O (only H2/STO-3G so far, same
incremental-verification-first pattern as every other piece of this
project); Pade continuation is known to be numerically delicate with
noisy/sparse input data in general -- this session's `Nomega=6`,
evenly-spaced, noise-free (analytic pipeline output, not experimental)
grid worked cleanly, but robustness with fewer/more points or a
noisier `Wc_reduced` has not been separately stress-tested.

**Committed**: `scrnc_sigma_driver.F`, plus this log.

---

## SESSION 2026-08-17 (continued 10) -- per-frequency-point orchestration + Sigma_c(p,q) accumulation built and verified end-to-end on H2/STO-3G/DRCCD debug case

Picked up right at "continued 9"'s resume point after a context reset.
Confirmed its open question first: no saved driving script exists
anywhere (`aces2/scrnc/`, `/blue/ufhpc/perera/scrnc_test/*/`, searched
for `run_scrnc*`/`*scrnc*driv*`) -- prior dynamic-omega validation runs
really were driven by hand, one point at a time (edit `frequency`,
rerun the full `xjoda...xscrnc` chain, repeat), confirming that
session's working assumption.

**Design decisions (Ajith, this session)**:
- Orchestration: external (not an internal Fortran loop inside
  Scrnc_driver.F), reusing a SINGLE shared working directory across all
  N frequency-point pipeline runs -- not a fresh subdirectory per point.
- Language: Fortran, not shell/python -- a standalone external driver
  program that spawns the existing, unmodified pipeline binaries as
  subprocesses via the Fortran 2008 `EXECUTE_COMMAND_LINE` intrinsic.
- Bare-Coulomb subtraction (`W_c = W - v`): do it now, not deferred --
  needed for `Sigma_c` to be physically meaningful, not just a wiring
  placeholder.
- External evaluation frequency `omega` in `Sigma_c(p,q;i*omega)`:
  `omega=0` for this session's debug/wiring pass; the real target is
  a proper real-axis Pade-continuation frequency scan, deferred as its
  own next step. The driver's architecture (see below) is built so that
  scan can reuse the SAME N pipeline runs already cached in this
  session's design -- no rearchitecture needed later, just populating
  `sigma_input`'s `Omega_ext` list with more points.

**New files**:
- `scrnc_sigma_prep.F` (new, part of the normal `xscrnc` build, called
  unconditionally from `Scrnc_driver.F` right after the existing
  `Scrnc_wc_ao2mo_test` diagnostic call): builds the bare AO 2-electron
  Coulomb matrix `Jj_bare` (via `Scrnc_get2ints`, the same routine
  `scrnc_prep_2ints.F` uses to prep `List_xxxx`, re-indexed into
  scrnc's own `Iao_pair=(Mu-1)*Nao+Nu` convention -- safe despite
  `Scrnc_get2ints`'s own different internal `indx2` convention because
  bare 2-electron AO integrals are already fully permutationally
  symmetrized at that routine's own construction, confirmed by
  inspection, not assumed), computes `Jj_corr = Jj_lrespn - Jj_bare`
  in place (`W_c = W - v`, valid since the AO->MO transform is linear),
  reuses the already-verified `Scrnc_wc_ao2mo` unchanged on `Jj_corr`,
  reduces the result to `Wc_reduced(p,q,m) = Wc_mo(p,m,q,m)` (the
  diagonal-in-second/fourth-index slice the self-energy formula
  actually needs -- `Nmo**3` instead of the full `Nmo**4` tensor), and
  writes `Wc_reduced` + the alpha-spin `SCFEVALA` orbital energies to a
  plain formatted file `wcreduc.dat`. Alpha-spin only -- inherited
  scope from `Scrnc_wc_ao2mo_test.F`, not a new decision.
- `scrnc_sigma_driver.F` (new, standalone `PROGRAM`, deliberately NOT
  linked into the normal ACES2/JOBARC-aware build -- compiled and
  installed separately, see below): reads `sigma_input` (N, Omega_c,
  Nomega, then Nomega external `Omega_ext` values) and `pipeline_cmds`
  (one line, the shell command for one pipeline pass -- kept external/
  path-agnostic so the SAME driver works for either the ring-CC or
  CCSD pipeline chain without editing Fortran source). Builds the
  quadrature grid via a local copy of `Scrnc_freq_grid`/`Scrnc_gaussl`
  (same multiply-duplicated-per-module convention as the original).
  Loops `K=1,N`: writes the `frequency` file for node `Xnode(K)`
  (unsquared -- matches the existing `h2_scrnc_drccd_imagfreq_0.3`
  convention, `Scrnc_driver.F` squares it internally), runs
  `pipeline_cmds` via `EXECUTE_COMMAND_LINE`, checks the exit status,
  reads back `wcreduc.dat`, and immediately folds node K's contribution
  into the running `Sigma_c(p,q;i*Omega_ext(j))` accumulation for EVERY
  requested `j` in one pass (only `Nmo**3` transient storage needed per
  node, never all N nodes' `Wc_reduced` simultaneously) -- since the
  expensive pipeline only depends on the inner quadrature frequency
  `omega'`, while the outer evaluation frequency `omega` only enters
  algebraically in the G0 denominator, exactly N pipeline runs suffice
  regardless of how many external `omega` points are requested, ready
  to support the deferred Pade-continuation scan without rearchitecture.
  Writes `sigma_c_output.dat` (real+imaginary `Sigma_c(p,q)` per
  requested `Omega_ext`).

**Build**: `scrnc_sigma_prep.F` picked up automatically by the normal
`gmake`/`gmake install` in `aces2/scrnc/` (wildcard-discovered like
every other `.F` file there) -- confirmed it does NOT get pulled into
the `xscrnc` binary's own duplicate-main hazard: `scrnc_sigma_driver.F`
also compiles in that same pass (it's a `.F` file in the same
directory) but its `PROGRAM` is never referenced by `scrnc_main`'s call
graph, so the linker never pulls its object's `MAIN__` into `xscrnc` --
confirmed via `md5sum`/size delta (only ~8KB growth, consistent with
just `Scrnc_sigma_prep`'s own new code, not an entire duplicate
program). `scrnc_sigma_driver` itself built separately, standalone, no
ACES2 library linkage: `ifx -i8 -o scrnc_sigma_driver
scrnc_sigma_driver.F`, installed as `aces2/bin/xscrncsigma` (new
executable name, doesn't collide with `xscrnc`). Re-hit the
`make`-vs-`make install` PATH trap ([[feedback_aces_build_install_trap]])
immediately -- first plain `gmake` left the OLD `aces2/bin/xscrnc`
untouched even though the local module-dir binary and `.o` files were
fresh; `gmake install` was required to actually update
`aces2/bin/xscrnc`, confirmed via `md5sum` matching between the local
and installed copies afterward.

**Verified end-to-end**, fresh directory
(`/blue/ufhpc/perera/scrnc_test/h2_scrnc_sigma_debug/`, seeded from
`h2_scrnc_drccd_uhf_sto3g`'s `ZMAT`/`GENBAS`, debug settings `N=4,
Omega_c=1.0, Nomega=1, Omega_ext(1)=0.0`,
`pipeline_cmds="xjoda;xvmol;xvmol2ja;xvprops;xvscf;xvtran;xintprc;
xringcc;xringcl;xscrnc"`):
- `xscrncsigma` ran all 4 subprocess pipeline passes and exited 0, no
  manual intervention.
- New `Scrnc_sigma_bare_ao` self-check: max bare-`v` AO-integral
  asymmetry = exactly `0.0` at every node.
- Pre-existing `Scrnc_wc_ao2mo_bfcheck` regression check unaffected:
  still `~2.2D-16`, confirming `Scrnc_sigma_prep.F` didn't disturb the
  already-verified AO2MO transform it reuses.
- Free independent sanity check noticed in the pipeline's own existing
  debug output: `Scrnc_get2ints`'s raw bare-AO-integral printout
  (e.g. `(11|11)=0.774605825`) matches the already-printed "RPA JJ
  matrix" (full `W`, e.g. `JJ(1,1)=0.77461131149`) to 5-6 significant
  figures, with the residual difference (`~5.5D-6`) matching the
  separately-printed "Response-only JJ matrix (before bare v added)"
  matrix's own magnitude (`~5.5D-6`) almost exactly -- independent
  confirmation, from data the pipeline was already printing for
  unrelated reasons, that `W = v + (small correlation piece)` for this
  system and that `Jj_bare`'s magnitude/convention is right, without
  needing to add any new instrumentation to check it.
- `sigma_c_output.dat`: `Sigma_c(p,q;i*0)` came out purely diagonal
  (off-diagonal elements exactly `0.0`) -- correct qualitative
  signature for H2/STO-3G's symmetric bonding/antibonding MO pair not
  mixing under this system's symmetry. Diagonal magnitudes
  (`Sigma_c(1,1)=7.95D-3 + 5.48D-3i`, `Sigma_c(2,2)=-8.44D-3 +
  6.34D-3i`) are physically reasonable in scale for an STO-3G
  correlation self-energy; nonzero imaginary part at `omega=0` is
  expected (not a bug) -- the contour integral here runs one-sided over
  `omega' in [0,Infinity)`, not symmetrized over `+/-omega'`, so there
  is no structural reason for the imaginary part to cancel at this
  stage (that would only happen after the real-axis analytic
  continuation, not yet built).

**Not yet done / next step**: the deferred real-axis Pade-continuation
frequency scan (populate `sigma_input`'s `Omega_ext` list with a proper
scan of points, fit `Sigma_c(p,q;i*omega_j)` to a Pade approximant, and
continue to the real axis for the quasiparticle equation
`eps_QP = eps_HF + Sigma_c(eps_QP) + Sigma_x - v_xc`-style solve) --
this session's driver architecture is ready for that (see "only N
pipeline runs regardless of Nomega" note above) but the actual Pade
fit + QP solve code does not exist yet. The bare-exchange `Sigma_x`
term also still not scoped at all (flagged, unchanged, since
"continued 5").

**Committed**: `scrnc_driver.F`, `scrnc_sigma_prep.F`,
`scrnc_sigma_driver.F`, plus this log. Debug test artifacts left on
disk (not committed, per project convention):
`/blue/ufhpc/perera/scrnc_test/h2_scrnc_sigma_debug/`.

---

## SESSION 2026-08-17 (continued 9) -- frequency quadrature grid implemented and verified (Gauss-Legendre + rational map to [0,Infinity)); NOT yet wired to drive per-frequency pipeline runs / Sigma_c accumulation

Next roadmap item after the AO2MO transform fix (continued 8, below):
the frequency quadrature grid for the contour-deformation Sigma_c(i*
omega) integral, previously flagged "not yet researched/decided"
(continued 5 entry). Ajith's call, asked directly: Gauss-Legendre with
a rational map to [0,Infinity), over a minimax grid (minimax needs
external tabulated points/weights this codebase doesn't have; GL+map is
simple, no new dependency, fits this project's whole correctness-first
approach so far).

Found `aces2/vdint/gaussl.f` already has a working Gauss-Legendre root/
weight finder (`GAUSSL(X1,X2,N,X,W)`, Numerical-Recipes-derived, JG
4/93) -- reused it rather than writing one from scratch. Per Ajith's
Fortran-typing correction earlier this session
([[feedback_aces_fortran_implicit_typing]]): no explicit Integer*8
anywhere in the new code, plain implicit typing throughout, matching
this file's own convention.

New files: `scrnc_freq_grid.F` (`Scrnc_freq_grid(N,Omega_c,Xnode,
Wnode)` -- builds the N-point [0,Infinity) grid via `omega'=Omega_c*
(1+x)/(1-x)` applied to GL's native [-1,1] nodes, with the matching
Jacobian on the weights; `Scrnc_gaussl` is a local copy of
`aces2/vdint/gaussl.f`'s routine, not a cross-module link -- same
multiply-duplicated-per-module convention as SYMEXP/DODIIS elsewhere in
this codebase, see GNUmakefile.src) and `scrnc_freq_grid_test.F`
(`Scrnc_freq_grid_test`, `_DEBUG_LVL0`-gated like
`Scrnc_bcktran_test.F` -- genuinely optional pure-numerics self-check,
unlike the AO2MO transform which Ajith said must always run
unconditionally). Wired into `Scrnc_driver.F` right after the AO2MO
transform call.

**Verified**: built clean via `gmake`/`gmake install` (matching the
lesson from continued 7/8 -- always the real build system, never a
hand-rolled compile), ran on the H2O/STO-3G/UHF/DRCCD case. Checked
against two closed-form integrals at N=8,16,24 and Omega_c=0.5,2.0:
`Integral 1/(1+omega**2) = Pi/2` converges to `1.7D-13` by N=24 (both
Omega_c scales); `Integral Exp(-omega) = 1` converges more slowly
(expected -- exponential decay is the harder case for a rational-map GL
grid) but still monotonically, down to `4.1D-9` at N=24,Omega_c=2.0.
Both integrals converge cleanly with increasing N at both Omega_c
scales -- no bugs, grid is correct. Also reran the AO2MO
`Scrnc_wc_ao2mo_bfcheck` in the same run to confirm no regression:
still `max|brute force - fast| = 1.065814D-14`, unchanged.

**Not yet done**: this grid is NOT yet used to actually drive anything
-- no code yet loops over `Xnode(1..N)`, runs the pipeline (or the
relevant sub-piece) at each imaginary frequency, collects `Wc_mo` at
each point, and accumulates the `Sigma_c(p,q)` sum with the GL weights
and the `1/(i*omega+i*omega'-eps_m)` denominator. Investigated how
multi-frequency runs work in this codebase first: `scrnc_main.F` calls
`Scrnc_driver` exactly ONCE per `xscrnc` invocation (no outer loop);
frequencies are read from a `"frequency"` file
(`Scrnc_set_vars.F:59-76`: `No_Freqs`, a real/imaginary flag, then
`No_Freqs` frequency values into the `Xfreq(200)`/`Freqinfo` common) --
but only the LAST-read value populates the scalar `Freq` actually
consumed downstream; `Xfreq`/`No_freqs` don't appear to be looped over
anywhere else in this module. So the working assumption (not yet
confirmed against how the earlier 2026-08-16 imaginary-axis W(omega)
validation was actually driven) is: multi-frequency point-by-point runs
happen via EXTERNAL orchestration (write a `"frequency"` file with a
single imaginary `Freq` value, run the full `xjoda...xscrnc` chain,
repeat per point) -- not an internal Fortran loop. Next step: confirm
this against how the 2026-08-16 dynamic-omega H2/H2O validation runs
were actually driven (check old pipeline.log/test-dir provenance for a
driving script), then build the actual per-point orchestration +
Sigma_c(p,q) accumulation on top of `Scrnc_freq_grid`'s nodes/weights
and the now-fixed `Scrnc_wc_ao2mo`.

---

## SESSION 2026-08-17 (continued 8) -- ROOT CAUSE FOUND AND FIXED: Jj_lrespn was implicitly Integer*8, not Double Precision (classic Fortran I-N implicit-typing trap on a name starting with 'J'); AO2MO transform now verified to machine precision and wired in unconditionally (CLOSED)

Used gdb (hardware watchpoints + direct DWARF address-of, both cross-
checked against raw RDI-register reads) to find the actual mechanism,
per Ajith's "complete the project" instruction. First attempt was a
DEAD END worth recording so it isn't repeated: setting a breakpoint at
`scrnc_wc_ao2mo_do_`'s function name (gdb auto-skips the prologue) and
reading `$rdi` immediately gave a STABLE, plausible-looking address --
but this was a false trail. At -O0, by the time gdb's default
prologue-skip lands, argument registers have already been spilled to
their home stack/Work-array slots and are being REUSED as scratch for
the first statement's own address computation -- `$rdi` at that point
had been silently repurposed to hold `Blk_mo`'s address (the routine's
own last internal scratch array), not `Jj_lrespn`'s. Chasing this gave
an entirely spurious "MKL dgemm writes into Jj_lrespn's memory"
watchpoint hit that looked like conclusive proof of aliasing but wasn't
-- a lesson for next time: never trust a raw register read at a
breakpoint as "the argument value" without cross-checking against a
DWARF-based `&variable` (or explicit pre-prologue address) lookup.

Redid it properly: `set language c` (needed -- Fortran-mode expression
parsing rejected C-style pointer casts) then `(long)&jj_lrespn` at
THREE levels of the call chain (`Scrnc_wc_ao2mo_test`, `Scrnc_wc_ao2mo`,
`Scrnc_wc_ao2mo_do`) all agreed exactly: `0x155554d2c628`, Fortran index
127498918 -- EXACTLY `Ijj_lrspn_bgn` (`=Memleft+1`, independently
confirmed via `maxcor=127498917` read directly off the live argument).
No aliasing, no overlap, address correct and stable throughout -- fully
retracts the "continued 7" watchpoint finding above.

**So the address was always right; the VALUE was wrong.** Broke exactly
at the read statement (`scrnc_wc_ao2mo.f:104`,
`Blk_ao(Mu,Nu) = Jj_lrespn(Iao_pair,Jao_pair)`) for the very first
iteration (`Mu=Nu=Iao_pair=Jao_pair=1`), confirmed the address seen
there matches the verified-true address exactly (delta=0), then did
`set language fortran; print jj_lrespn(1,1)` -- got back
`4605170097222381441`, a huge INTEGER, not a small real value. This is
the exact int64 bit-reinterpretation of a normal double: `struct.unpack
('<q', struct.pack('<d', 0.777650825))` in Python gives
`4605179675476665134` -- same leading digits, matching magnitude and
sign (~4.6D+18), i.e. the compiled code was reading the 8 raw bytes at
that address as an INTEGER, not a REAL*8. This is EXACTLY the same
"~4.6D+18" magnitude seen in the very first "DBG Xtmp" print from the
2026-08-17 "continued 6" entry, and the same order of magnitude as
every garbage `Wc_mo`/bfcheck value seen throughout this whole
investigation -- one single root cause explains all of it.

**Root cause**: `Jj_lrespn`'s name starts with `J`, which falls in
Fortran's implicit I-N range (default INTEGER). This file's own
`Implicit Double Precision (A-H,O-Z)` statement covers ONLY A-H and
O-Z -- it does NOT touch I-N. A bare `Dimension Jj_lrespn(Nao_pairs,
Nao_pairs)` (present in all 4 routines that take it as a dummy arg:
`Scrnc_wc_ao2mo`, `Scrnc_wc_ao2mo_do` in `scrnc_wc_ao2mo.F`, and
`Scrnc_wc_ao2mo_test`, `Scrnc_wc_ao2mo_bfcheck` in
`scrnc_wc_ao2mo_test.F`) only sets the ARRAY SHAPE, not the type -- so
`Jj_lrespn` was silently `Integer*8` in every one of these new routines.
Every element READ converted the raw IEEE-754 bit pattern into its
integer value, then (since the destination, e.g. `Blk_ao`, IS correctly
Double Precision under A-H,O-Z) an INTEGER-TO-REAL numeric conversion
ran on that huge integer, producing the ~1e18-1e19-scale garbage seen
everywhere. Because -i8 makes `Integer` and `Double Precision` both
8 bytes, array-indexing/stride arithmetic was completely unaffected --
which is exactly why every address-based check across this whole
investigation (hand-verified offsets, gdb address-of, watchpoints on
the correct address) kept coming back clean: the bug was never about
WHERE the data lived, only about how its bits were being interpreted
once read. This also fully explains why "brute force" and "fast"
`Wc_mo` always agreed with each other despite both being wrong --
`Scrnc_wc_ao2mo_bfcheck` has its OWN independent `Jj_lrespn` dummy arg
with the identical bare-`Dimension` bug, so both paths independently
misread the same underlying bits the same way.

The codebase's own established convention for this exact array already
protects against this everywhere else -- `scrnc_jj_rcc_assemble.F:6`:
`Double Precision JJ_lrespn(Nao_pair,Nao_pair)`, an explicit type
declaration, not a bare `Dimension`. The four new routines in this
session's two new files simply missed that convention.

**Fix**: changed `Dimension Jj_lrespn(Nao_pairs,Nao_pairs)` to
`Double Precision Jj_lrespn(Nao_pairs,Nao_pairs)` in all four routines
(`scrnc_wc_ao2mo.F` x2, `scrnc_wc_ao2mo_test.F` x2). Rebuilt clean via
`gmake clean && gmake && gmake install` in `aces2/scrnc/` (production
flags, no debug/-O0 overrides), md5-verified the installed binary.
Reran the H2O/STO-3G/UHF/DRCCD ring-CC case fresh
(`h2o_wc_ao2mo_fixed/`): `@-Scrnc_wc_ao2mo_bfcheck: max|brute force -
Scrnc_wc_ao2mo| = 1.065814D-14` (machine precision), `brute value =
fast value = 5.061380D+00` (a physically sane magnitude, not garbage).
**Fixed and verified.**

**Also this session (Ajith's instruction)**: the `#ifdef _DEBUG_LVL0`
guard around the `Scrnc_wc_ao2mo_test` call in `Scrnc_driver.F` was
removed -- this transform is a required step of the self-energy
pipeline (every downstream step needs `Wc_mo`), not an optional
diagnostic, so it now always runs unconditionally. `Scrnc_wc_ao2mo_test`
still bundles the brute-force cross-check every call, which is fine at
STO-3G toy-system scale (this project's whole scope so far) but would
not scale to a real system -- worth splitting the always-needed
transform from the diagnostic brute-force check before this goes beyond
toy systems; not done this session (not asked, and out of scope of the
immediate fix).

**Status: this specific bug is CLOSED.** The AO->MO transform is
production-ready and wired into the real pipeline. The broader
self-energy roadmap from the "continued 5" entry above is NOT complete
-- still not started: the frequency quadrature grid choice, the Pade
continuation implementation, G0 pole placement/eta convergence, the
bare-exchange Sigma_x term, and the final quasiparticle equation solve.
"Complete the project" (Ajith, this session) refers to clearing this
specific blocker, not the full GW self-energy implementation -- confirm
with Ajith which of those pieces to tackle next.

---

## SESSION 2026-08-17 (continued 7) -- ruled out the build-flags theory for the Scrnc_wc_ao2mo corruption; a same-binary reproducibility check confirms it's fully deterministic (100% reproducible), not a heisenbug -- clear to use a debugger next

Ajith's instruction after the terminal freeze: stop hand-linking a standalone
test harness and instead compile/link `Scrnc_wc_ao2mo`/`Scrnc_wc_ao2mo_test`
through this tree's own `GNUmakefile`/`GNUmakefile.src` (`64BIT=1` ->
`-i8` -Zp8, `-Bstatic -lmkl_intel_ilp64` ILP64 MKL), same as every other
scrnc file, so nothing about the standalone "verified correct in isolation"
tests' build could be silently mismatched vs. the real pipeline. Rationale
was sound (this project's whole bug history is ifx-flag-sensitivity), so
worth checking properly before more manual debugging.

**Result: theory ruled out.** Re-enabled the `Scrnc_wc_ao2mo_test` call in
`Scrnc_driver.F` (`_DEBUG_LVL0`-guarded, was `CSSS`-disabled), rebuilt via
plain `gmake` + `gmake install` in `aces2/scrnc/` (needed `module load
intel/2025.1.0` first -- `ifx` wasn't on `PATH` in this shell; confirmed via
`gmake -n` that the real compile line is
`ifx -c -Wall -traceback -O2 -Zp8 -zero -i8 -fPIC scrnc_driver.f`, i.e.
`-i8` really is always on, exactly as Ajith said, no manual flag-wrangling
needed), md5-verified `aces2/bin/xscrnc` matches the freshly-linked local
binary (per [[feedback_aces_build_install_trap]]). Reran the ORIGINAL
ring-CC (DRCCD) H2O/STO-3G/UHF test case (`ZMAT` from
`h2o_scrnc_wc_ao2mo_check/`, fresh dir `h2o_scrnc_wc_ao2mo_check2/`) via the
exact same manual `xjoda,xvmol,xvmol2ja,xvprops,xvscf,xvtran,xintprc,
xringcc,xringcl,xscrnc` chain the earlier verification runs used (NOT
`xaces2` -- tried that first, it routes `CALC=DRCCD` through `xvcc` instead
of `xringcc` and dies in `xvcc` with an unrelated `GETLST` list-size error;
this whole DRCCD path apparently isn't reachable via the standard top-level
driver, only the manual per-binary chain -- worth a separate look someday
but not this bug). Corruption reproduced identically in kind: `Wc_mo`
elements ~1e19 (should be O(1-100)), brute-force and fast-transform paths
agree with each other (both give the same garbage value), confirming the
corruption is upstream of both algorithms, in what `Jj_lrespn` reads as.

**Apparent non-determinism claim RETRACTED (was a bad comparison, not a
real finding)**: initially compared this run's worst-mismatch location
`(p,q,r,s)=(7,1,7,6)` against the earlier session's `(6,6,6,6)` and
concluded run-to-run variance -> uninitialized memory. Ajith immediately
pushed back, correctly: this is single-threaded serial code, identical
input MUST give identical output, full stop -- if it doesn't, something
external is wrong (different binary, corrupted files), not "expected
randomness." Ajith's instruction: `gmake clean`, `gmake install` fresh
from `aces2/scrnc/`, then run the INSTALLED binary twice and compare
directly, rather than comparing across sessions/builds.

Did exactly that: `gmake clean` (full), `gmake` (fresh compile of all of
scrnc including `scrnc_wc_ao2mo.F`/`scrnc_wc_ao2mo_test.F`), `gmake
install`, md5-verified `aces2/bin/xscrnc` matches the local build
(`3efb006494128b2431ee5becdcd396fe`). Ran the SAME installed binary,
same `xjoda...xscrnc` chain, in two independent fresh directories
(`h2o_wc_ao2mo_repro_A`, `h2o_wc_ao2mo_repro_B`) from identical
`ZMAT`/`GENBAS`. `diff -a` of the two full pipeline logs: the ONLY
differences anywhere in ~20k lines are `@ACES_INIT-DEBUG`/`@INTS_-DEBUG`
raw heap-address printouts (`base addr=`, `i0=`, `ICORE(...)  addr=` --
ordinary ASLR noise, a fresh process gets a different heap base every
launch, expected and harmless) and sub-hundredth-second wallclock timing
noise. Every physics number, including the `Scrnc_wc_ao2mo_bfcheck`
block, is BIT-IDENTICAL between the two runs: `worst (p,q,r,s) = 7 1 7
6`, `brute value = fast value = -2.459861D+19`, both runs, exactly.

**Corrected conclusion**: the corruption is fully deterministic and
100% reproducible given the same binary + same input -- there is no
run-to-run randomness. The earlier "(6,6,6,6) vs (7,1,7,6)" mismatch
was comparing against a stale run from BEFORE this session's rebuild
(different code -- the old run still had the temporary `DBG` print
instrumentation compiled in, and was very likely built via a different,
not-fully-`gmake`-clean path pre-freeze) -- not a real property of the
current code. The "uninitialized memory" theory and the associated
`Nao`/`Nao_pairs` out-of-bounds-read speculation are RETRACTED --
spot-checked `Nao_pairs` threading through `Scrnc_driver.F` anyway
(`Nao_pairs=Nao*Nao` set once at line 115, passed unchanged through
`Ijj_lrspn_bgn`'s sizing and the `Scrnc_wc_ao2mo_test` call site) and
found no mismatch, consistent with this being a dead end. Back to
believing this is a genuine deterministic bug (wrong index, aliasing,
or similar) -- good news, since deterministic + 100% reproducible means
a debugger watchpoint (the (a) option from the entry below) will show
the same thing every time, no need to worry about it being a heisenbug.

---

## SESSION 2026-08-17 (continued 5) -- self-energy architecture DECIDED (Ajith: contour-deformation route), initial data-flow analysis done, implementation not yet started (RESUME HERE)

Ring-CC is closed (see entry below). Moving to the "next level" per
Ajith's direction: build `Sigma_c(omega)` on top of the now-verified
`W(omega)`. The 2026-08-16 entry (below) had already identified this as
a genuine fork needing Ajith's call and left it unresolved -- asked
explicitly this session.

**Decision (Ajith, 2026-08-17): route 1, contour-deformation** --
imaginary-axis `W(i*omega')` (already validated on both CCSD and
ring-CC paths, static+dynamic, real+imaginary axis, H2 and H2O) +
frequency quadrature + Pade analytic continuation to the real axis.
Rejected route 2 (sum-over-poles via Hbar diagonalization) specifically
because it would diagonalize/upfold Hbar, moving closer to the
diagrammatic/upfolding GW-CC formalism style the manuscript's novelty
claim needs to stay distinct from (see
[[project_scrnc_gw_manuscript_goal]]).

**Standard equation** (imaginary axis, orbital/MO basis, G0 diagonal
in the MO basis since G0_mn(i*nu) = delta_mn/(i*nu-eps_m)):
```
Sigma_c(p,q; i*omega) = -(1/2*pi) * Integral d(omega') Sum_m
     [1/(i*omega+i*omega'-eps_m)] * W_c(p,m,q,m; i*omega')
```
where `W_c = W - v` (correlation part only, bare Coulomb subtracted)
and `W_c(p,m,q,m'; i*omega')` is a 4-MO-index screened-interaction
tensor -- structurally analogous to a 2-electron integral, but with the
bare Coulomb kernel replaced by the frequency-dependent `W_c`.

**Key finding this session -- the AO->MO transform is 4-index, not
2-index, but collapses to two big matmuls.** Checked how `JJ_lrespn`
(scrnc's `W(omega)` in AO-pair basis) is actually dimensioned:
`Nao_pairs = Nao*Nao` (`scrnc_driver.F:115`) -- the FULL, unpacked
square of AO indices (not a symmetry-triangle-packed set), and
`JJ_lrespn` is genuinely `Nao_pairs x Nao_pairs` = `Nao**2 x Nao**2`
(confirmed via `scrnc_jj_rcc_assemble.F:6`,
`I_pairs=J_pairs=Nao_pair`). So `JJ_lrespn(mu*nu, lam*sig)` already
carries all four AO indices via two compound Kronecker-style pair
indices -- this means the needed AO->MO transform is the SAME
computational class as a standard 2-electron-integral AO2MO transform
(4 AO indices -> 4 MO indices), NOT the simpler 2-index density-like
transform `scrnc_bcktran_left.F`/`scrnc_i1tran.F`'s ICASE=3 handle
(those were built for a different, now-apparently-wrong hypothesized
use case -- the 2026-08-16 "bcktran target"/"AO-vs-MO self-energy
direction" open question is resolved by this: neither existing
transform kernel is directly the right shape for this job).

Because `JJ_lrespn` is stored UNPACKED (`Nao_pairs=Nao**2`, no AO
symmetry packing), the 4-index transform does NOT need the full
standard ERI quarter-transform algorithm's index-permutation
machinery -- it collapses to two ordinary matrix multiplies using a
Kronecker-product transformation matrix `T = C tensor C`
(`T((mu,nu),(p,q)) = C(mu,p)*C(nu,q)`, size `Nao**2 x Nao**2`):
```
Step 1 (row/left pair):  W'(pq, lam*sig) = T^T * JJ_lrespn
Step 2 (col/right pair): W_c_mo(pq, rs) = W'(pq, :) * T
```
i.e. `W_c_mo = T^T * JJ_lrespn * T`, two `Nao**2 x Nao**2` Xgemm calls
(reusing the same BLAS scrnc already calls elsewhere), NOT a bespoke
index-permutation transform. `T` itself is cheap to build (Kronecker
product of the existing SCF eigenvector matrix `C` with itself).

**Scaling caveat, flagged not yet resolved**: `Nao_pairs=Nao**2` means
`JJ_lrespn` is `O(Nao**4)` in memory -- fine for STO-3G H2/H2O (this
whole project's test systems so far) but will not scale to any
realistically-sized molecule the eventual manuscript needs. Decided
(consistent with how every other piece of this project has been built)
to prove correctness on the existing STO-3G toy systems first and treat
storage-scaling as a separate, later optimization -- NOT to redesign
the storage convention now. Revisit before attempting any real-system
GW numbers.

**AO-pair index convention, pinned down (was the real blocker for
writing the transform correctly)**: `scrnc_prep_2ints.F:68-77` builds
the authoritative `Irrep_ao_pairs` map (shared via
`Common/Sym_pqxx_ints/`, reused by `scrnc_lrspn_rcc_main.F` and others)
via `Do k=1,Nbfns; Do l=1,Nbfns; Ithru=Ithru+1; ...`, i.e.
**`Iao_pair = (Mu-1)*Nao + Nu`**, with `Mu` (=`k`) the OUTER/slow loop
and `Nu` (=`l`) the INNER/fast loop. This is the reshape needed to
treat `JJ_lrespn(Iao_pair,Jao_pair)` as a genuine 4-AO-index tensor
`W(Mu,Nu,Lam,Sig)` for the transform above -- confirmed by direct
inspection, not assumed (guessing this backwards would silently
transpose the transform, same failure class as several already-found
bugs in this project).

**Not yet started**: actually writing the `T=C(x)C` construction and
the two-Xgemm transform routine; the frequency quadrature grid choice
(minimax vs Gauss-Legendre-after-transform, standard choices in
numerical-GW contour-deformation codes -- not yet researched/decided);
the Pade continuation implementation; the G0 pole placement/eta
convergence; and the final quasiparticle equation solve
(`eps_QP = eps_HF + Sigma_c(eps_QP) + Sigma_x - v_xc`-style, needs the
bare-exchange `Sigma_x` term too, which is a separate, simpler,
frequency-independent piece not yet scoped at all). This entry is the
resume point for all of that.

---

## SESSION 2026-08-17 (continued 4) -- RE-VERIFIED CLEAN: ring-CC H2O JJ symmetry re-checked fresh against current code, no open bug (correcting a stale mischaracterization made earlier this same session)

Earlier in this session I described the H2O ring-CC JJ asymmetry
(2026-08-14) as "still open" -- that was wrong, based on only having
read the stale `## RESUME POINT (as of 2026-08-14...)` section near the
end of this file without having yet read the 2026-08-15 entry above it
that actually fixed it (missing `Call Symmet2` in `scrnc_jj_rcc_lrspn.F`/
`scrnc_jj_cc_lrspn.F`, commit `35cdb853`, "This closes out the entire
scrnc JJ-asymmetry investigation, this time for real").

On Ajith's instruction to "clean out ringcc" before moving on, did NOT
just trust that old log entry -- re-verified against the CURRENT
codebase (the fix has since moved into `scrnc_jj_rcc_assemble.F`, split
out during the 2026-08-15-evening dynamic-omega refactor):

- Confirmed `Call Symmet2(JJ_lrespn,Nao_pair)` present in both
  `scrnc_jj_rcc_assemble.F:147` and `scrnc_jj_cc_lrspn.F:385`.
- Ran 3 fresh, full pipelines (`xjoda...xringcc,xringcl,xscrnc`) on
  H2O/STO-3G/UHF with the CURRENT code, parsed the real 49x49 "RPA JJ
  matrix" block programmatically (4-cols-per-header ACES2 `output`
  format), checked every symmetric pair:
  - DRCCD static (`h2o_scrnc_drccd_postrefactor_check`): max
    |JJ(r,c)-JJ(c,r)| = 0.0
  - RCCD static (`h2o_scrnc_rccd_postrefactor_check`): 0.0
  - DRCCD dynamic, omega=0.3 a.u. real
    (`h2o_scrnc_drccd_omega03_check`): 0.0 -- **this closes a real gap**:
    dynamic-omega correctness had only ever been verified on H2 before
    (2026-08-14 entry, "Verified (H2/STO-3G/UHF/DRCCD...)"), never on
    H2O. Now confirmed on H2O too.
- No errors/non-convergence in any of the three runs (checked directly,
  not just assumed).
- Confirmed the two minor indexing-typo fixes (`scrnc_driver.F:149-150`,
  `scrnc_form_lrspn_rcc_lists.F:61`) are still present, no regression.
- `git status` on `aces2/scrnc/`, `aces2/ringcc/`, `aces2/ringcl/`:
  clean, nothing uncommitted (only an unrelated untracked file from the
  2026-08-16 vcceh independent-pilot thread,
  `aces2/scrnc/independent_pilot/ci_like_vs_lr_polarizability.py` --
  not ring-CC, not touched).

Test artifacts left on disk (not committed, per project convention):
`/blue/ufhpc/perera/scrnc_test/h2o_scrnc_drccd_postrefactor_check/`,
`h2o_scrnc_rccd_postrefactor_check/`, `h2o_scrnc_drccd_omega03_check/`.

**Verdict: ring-CC has no open correctness issues.** Static and dynamic
paths, both DRCCD and RCCD variants, all exactly symmetric on H2O with
the current code. Clear to move on.

---

## SESSION 2026-08-17 (continued) -- RESOLVED by Ajith's read: only the CI-like piece is required to be self-symmetric; CI-like's raw asymmetry in vcceh is a confirmed real bug (RESUME HERE)

Read `spnspn_jcp.pdf` directly (via Dropbox) and pulled the actual text
of Eqs. (26)-(31), not just the earlier paraphrase. Precise statement:

- **Eq. (30)** splits the quadratic model into three terms: term A =
  `<0|(1+L)(H_N'-E_N')T_N|0>`, term B = `<0|(1+L)(H_N-E_N)T_N'|0>`
  (A+B = the CI-like piece), plus a third "quadratic correction" term
  `<0|(1+L)[H,T_N T_N']|0>`. Swapping N<->N' swaps A and B into each
  other, so **A+B is symmetric by algebraic construction alone** --
  no physics/derivative argument needed, holds explicitly "if the
  excitation space is restricted to singles and doubles" (line
  ~620-624 of the paper), i.e. at CCSD truncation.
- The quadratic correction term (third term) has **no such structural
  guarantee on its own**. Only the FULL Eq. (30) = Eq. (26) (all three
  terms combined) is proven symmetric, and that proof is the separate,
  weaker "elementary considerations" argument right after Eq. (26)
  (line 547): it's an exact mixed second energy derivative, and mixed
  partials of energy commute. That argument does not apply to any one
  piece in isolation.

**Conclusion**: Ajith's recollection was correct -- only the CI-like
piece is required to satisfy `f(N,N')=f(N',N)` (i.e. JJ(N,N')=JJ(N',N))
on its own. This settles the flip-flop below: the "corrected"/PySCF-
docstring framing (SESSION 2026-08-16 night, below) was the
overcorrection; the original [[feedback_eomcc_lr_eigenvectors_never_equal]]
reading stands, but scoped precisely to the CI-like piece only.

Cross-checked against the actual `vcceh` code the same session (not
just the paper): confirmed in `aces2/vcceh/secprp.F:368-421` that terms
A and B (`XPROP(IPERT,JPERT)` and `XPROP(JPERT,IPERT)`) are BOTH
independently computed via a full `DO 40 IPERT/DO 50 JPERT` double loop
over the entire NOPERT x NOPERT range (no `J>=I` restriction) -- this is
NOT a "compute one triangle, mirror it" shortcut like NWChem's (see
CALCQUAD1-vs-NWChem entry below). `SYMMET2`
(`aces2/libr/symmet2.f:20-24`) then does
`A(I,J)=0.5*(A(I,J)+A(J,I))` -- an average of two independently-solved
numbers, not a mirror-copy. Since terms A and B are proven to be exactly
equal under exchange by the paper's own algebra, any nonzero raw
pre-`SYMMET2` difference between `XPROP(IPERT,JPERT)` and
`XPROP(JPERT,IPERT)` for the CI-like piece is a **genuine implementation
bug** in `secprp.F` (or upstream in `OPRBAR`/the perturbed-T-amplitude
lists it reads), not expected/tolerable behavior.

Same logic does NOT extend to `calcquad1.f`'s raw asymmetry (the
~34.5% CALCQUAD1 number below) -- the paper gives no term-by-term
symmetry guarantee for the quadratic correction term alone, so its raw
asymmetry is not on its own diagnostic of a bug the way CI-like's is.
The already-established CALCQUAD1-vs-NWChem post-`SYMMET2` agreement
(4-5 sig figs, entry below) remains the best evidence that piece is
correctly implemented overall.

**Next step (not yet started)**: chase the CI-like term A/B mismatch
itself. Both terms are computed via the same `OPRBAR` call
(`secprp.F:350-353`, `DOLEFT`/`DORIGHT` branches) reading perturbed
T-amplitude lists `372+ISIDE` -- worth checking whether this is the
same list/scratch-record collision family as the already-fixed
ring-CC (`scrnc_jj_rcc_lrspn.F`) and CCSD-response
(`scrnc_lrspn_cc_main.F`) bugs, or something new in `secprp.F`/`OPRBAR`
itself. `secprp.F` is the general vcceh property-tensor path (not
scrnc-specific), so a fix here has broader blast radius than the
scrnc-local fixes -- confirm scope before touching it.

---

## SESSION 2026-08-17 (continued 2) -- CHECKED, NOT RESOLVED: `OPRBAR`/`LINEQY` wiring itself has no structural bug, but the "just numerical residual" explanation for the 7.6% gap is now RULED OUT (RESUME HERE)

Two things verified this pass, on Ajith's direct ask to check for a real
coding error before trusting `Symmet2` as the explanation:

**1. No structural bug in the `OPRBAR`/`LINEQY`/`Symmet2` wiring.**
`secprp.F` calls `OPRBAR` twice, not once (missed the first call two
entries ago): once inside the amplitude-solving loop (~line 197,
`DORIGHT=ISIDE.EQ.1`/`DOLEFT=ISIDE.EQ.2`) which writes the raw RHS to
list 373 and is then solved in place by `LINEQY`/`PLINEQY` into the
converged perturbed-T amplitude; once more in the property-construction
section (~line 350, flags deliberately **reversed**:
`DOLEFT=ISIDE.EQ.1`/`DORIGHT=ISIDE.EQ.2`), which writes the L-bar
operator (a direct, non-iterative construction -- confirmed `LINEQY`/
`PLINEQY` never touch list 374 anywhere) straight to list 374. With
`ISIDE` hardcoded to 1, this reversal is what lets a single fixed
`ISIDE` populate both 373 (solved) and 374 (direct) correctly. Confusing
to read, but deliberate and correct -- not a record-collision bug like
the two already-fixed scrnc-local ones.

**2. Checked whether the amplitude solve behind the ~7.6%
`h2o_ccsd_polar_check_asym` gap (two entries above) actually converged.**
Raw pipeline log
(`/blue/ufhpc/perera/scrnc_test/h2o_ccsd_polar_check_asym/pipeline.log`):
all three perturbations (DIPOLE_X/Y/Z) hit
`@LINEQY-I, Derivative amplitudes converged after 8 iterations.`
(lines 12597/12709/12821), with the intervening `@CHECKSUM` trace
showing the residual shrinking cleanly from ~0.5 to ~1e-13 beforehand --
consistent with genuinely reaching `LINEQY`'s own `CUTOFF=1.D-14`
(`lineqy.F:63`). No non-convergence, no max-iteration cutoff hit.

**This is decisive against, not for, "Symmet2 is just cleaning up
ordinary numerical residual."** A solve converged to 1e-14 should leave
term-by-term noise at roughly that scale. A 7.6% raw asymmetry is six-
plus orders of magnitude too large to be explained by finite iterative-
solve precision. So: confirmed no coding error in the wiring checked so
far, but ALSO no longer have a benign explanation for the gap's size --
"just numerical noise, nothing to fix" is not yet an earned conclusion
either. Status is genuinely open, not resolved either direction.

## SESSION 2026-08-17 (continued 6) -- Scrnc_wc_ao2mo (AO->MO transform) written and verified correct in isolation; wiring it into the real pipeline hit an unexplained, reproducible corruption -- disabled, not fixed (RESUME HERE)

Implemented `aces2/scrnc/scrnc_wc_ao2mo.F` (new): `Scrnc_wc_ao2mo` /
`Scrnc_wc_ao2mo_do`, the two-block-Xgemm AO->MO transform of `Wc`
(`JJ_lrespn`) designed in the entry above, using the confirmed
`Iao_pair=(Mu-1)*Nao+Nu` convention. All loop/index variables kept in
the I-N implicit-integer range on purpose (see the file's own header
comment) -- an earlier attempt using `P,Q,R,S,Sig` (which fall in the
O-Z implicit-DOUBLE-PRECISION range under this file's
`Implicit Double Precision (A-H,O-Z)`) was a REAL bug, caught via a
brute-force cross-check (`Scrnc_wc_ao2mo_bfcheck` in the new
`scrnc_wc_ao2mo_test.F`, same methodology as `Scrnc_bcktran_left.F`'s
own brute-force check) and fixed.

**Verified correct in isolation**: two standalone test programs
(dedicated, non-shared arrays, linked directly against
`Scrnc_wc_ao2mo.F` + `aces2/crust/xgemm.f` + MKL, bypassing all ACES2
Work-array/GETREC machinery) at Nao=Nmo=2 and Nao=Nmo=7 (matching the
real H2O/STO-3G production size) both agree with an independent
quadruple-nested-loop brute-force reference to ~1e-11 (machine
precision), for every element. The core transform algorithm is
believed correct.

**NOT resolved**: wiring it into the real pipeline (via
`Scrnc_wc_ao2mo_test`, called from `Scrnc_driver.F` guarded by
`_DEBUG_LVL0`, using the REAL `Jj_lrespn`/SCF eigenvectors) reproducibly
corrupts the result -- `Wc_mo` elements come out ~1e18-1e19 in
magnitude (should be O(1-100)). Chased extensively this session and
ruled out, in order: (1) a genuine implicit-typing bug (`P,Q,R,S,Sig`
as O-Z-range doubles -- real bug, fixed, but corruption persisted
unchanged); (2) `-i8`-vs-explicit-`Integer`-declaration kind mismatch
(renamed all loop vars to the I-N range instead of declaring
`Integer` explicitly -- ALSO had zero effect on the corruption, same
exact wrong value at the same location before and after, ruling this
out as the cause despite being a plausible mechanism in general); (3)
`ifx -O2` miscompilation, this project's own established failure class
per `feedback_aces_debug_try_low_opt_first` -- rebuilt
`scrnc_wc_ao2mo.o` at `-O0` by hand (matching the real `cpp`+`ifx`
invocation from `make -n`), corruption persisted identically; (4)
memory offset/aliasing between `Jj_lrespn`, `C`, `Wc_mo`, and
`Scrnc_wc_ao2mo`'s internal scratch, in the real `Scrnc_driver.F`
Work-array layout -- hand-verified the absolute offsets of every array
involved, found no overlap; (5) a naive "large local test array
overflows the 8MB stack" false lead in one intermediate reproduction
attempt -- confirmed via `ulimit -s` (8192 KB) vs a 16MB test array,
irrelevant to the real pipeline since ACES2's Work array is a COMMON
block, not stack-allocated.

**Most specific finding, still unexplained**: instrumented directly
inside `Scrnc_wc_ao2mo_do`'s Step-1 extraction loop. `Jj_lrespn(1,1)`
printed via a `Write` statement shows the correct value (0.8131535...,
matching the independently-printed "RPA JJ matrix" block earlier in
the same log) -- but assigning that SAME expression,
`Xtmp = Jj_lrespn(Iao_pair,Jao_pair)`, into a plain scalar `Xtmp` ONE
LINE EARLIER in strictly sequential `-O0` code gives ~4.6D+18 for the
identical `(Iao_pair,Jao_pair)=(1,1)`. I.e. the identical array-element
read expression evaluates to two different values a few nanoseconds
apart in the same loop iteration, with nothing in between that should
change `Iao_pair`/`Jao_pair`/`Jj_lrespn`. Confirmed via a PRE/POST
check that the destination (`Blk_ao`/`Xtmp`) genuinely starts at 0.0
(correctly zero-initialized by `-zero`) and becomes ~4.6D+18 only after
the assignment executes -- so it is the READ side that's wrong at that
call site specifically, not stale/uninitialized destination memory.
Never got further than this before deciding to stop rather than keep
guessing.

**Decision, this session**: rather than continue blind trial-and-error
(already 5 ruled-out hypotheses), disabled the driver-level call
(`Scrnc_driver.F`'s `Scrnc_wc_ao2mo_test` hook, now `CSSS`-commented
with an explicit "DO NOT re-enable" note) and cleaned all temporary
debug instrumentation out of `Scrnc_wc_ao2mo.F` itself, leaving it in
its believed-correct, standalone-verified state. Confirmed NO
regression from any of this session's changes: reran the H2O DRCCD
static/dynamic and CCSD control cases fresh
(`h2o_scrnc_drccd_postrefactor_check`, `h2o_scrnc_ccsd_postclean_check`)
-- both still exactly symmetric (max diff 0.0), matching the "ring-CC
is clean" verdict from earlier this session.

**Not yet done / real next step**: root-cause the read-inconsistency
above. Candidate approaches not yet tried: (a) use an actual debugger
(gdb/idb) with a watchpoint on `Jj_lrespn`'s memory address rather than
more printf-style instrumentation -- would show definitively whether
something writes to that address between the two reads, or whether
this is a codegen artifact; (b) write `Jj_lrespn` out to a small
dedicated JOBARC record via `PUTREC`/`PUTLST` right where it's known
good, then have `Scrnc_wc_ao2mo_test` read it back via `GETREC` fresh
(bypassing whatever call-chain/argument-passing path is implicated)
and compare; (c) check whether this is specific to `_DEBUG_LVL0` builds
specifically (i.e. try reproducing with a non-`_DEBUG_LVL0` build, in
case some OTHER `_DEBUG_LVL0`-guarded code earlier in `Scrnc_driver.F`
is the actual culprit writing where it shouldn't, and this transform
code is just an innocent bystander that happens to read the clobbered
region). Do not re-attempt the already-ruled-out hypotheses above
without new evidence.

---

**Next step (not yet started, this is now the concrete blocker)**:
directly inspect what `XPROP(IPERT,JPERT)` (`FACT * SDOT(T(IPERT)[373],
Y(JPERT)[374])`) is actually supposed to represent term-by-term. Two
live hypotheses, not yet distinguished:
1. A genuine implementation bug specifically in how `T(IPERT)` (373) or
   `Y(JPERT)` (374) is built/read for CROSS pairs (IPERT != JPERT) --
   the wiring check above only confirmed 373/374 get populated at all,
   not that the SDOT contraction combines them correctly for every
   ordered pair.
2. The paper's symmetry proof (Eq. 30's term structure, valid
   at CCSD truncation) may be a claim about the SUM
   `term_A(N,N')+term_A(N',N)` being invariant under simultaneous
   exchange, not a claim that `term_A(N,N')` alone equals `term_A(N',N)`
   alone -- in which case a raw single-term `XPROP(I,J)` before
   `Symmet2` is not expected to already be close to `XPROP(J,I)`, and
   the "7.6%" framing itself (comparing two individually-meaningless
   half-quantities) may be measuring the wrong thing. Needs re-deriving
   exactly what `FACT*SDOT(T(I),Y(J))` corresponds to in Eq. (30)'s
   term structure before concluding either way.
Best concrete test: a finite-field CCSD energy check (vary DIPOLE_X and
DIPOLE_Y field strengths independently, finite-difference the energy)
for just the H2O-asym case's alpha_xy vs alpha_yx pair, as an
independent ground truth outside the OPRBAR/Symmet2 machinery entirely
-- would settle whether 7.6% raw asymmetry is expected-but-mislabeled
or a real bug, without more Fortran-side guessing.

---

## SESSION 2026-08-17 (continued 3) -- RESOLVED: the raw pre-Symmet2 asymmetry is expected by construction, not a bug -- Symmet2 supplies the missing second term of the formula, not numerical cleanup

Re-derived from the paper's OWN sum-over-states form (Eqs. 36-40, not
just the operator-form Eq. 30) and it explains every observation
cleanly.

Eq. (40)'s CI-like expression, written as a dynamic (frequency-
dependent) sum-over-states property, has the standard two-resonance-
term structure:
```
E_NN'(w) = Sum_l [ <0|(Hbar_N-E_N)|l><l|H_N'|0> / (E-E_l-w) ]   <- "+w" term = term_A(N,N')
         + Sum_l [ <0|(Hbar_N'-E_N')|l><l|H_N|0> / (E-E_l+w) ]   <- "-w" term = term_A(N',N)
```
At w=0 (the static case, where the asymmetry is seen) both denominators
become identical. `secprp.F` computes only ONE resonance term per
(IPERT,JPERT) pair -- `XPROP(I,J) = FACT*SDOT(T(I)[373, converged
perturbed-T], Y(J)[374, direct L-bar operator])` -- and sets `FACT=2`
(`secprp.F:123`, the w=0/`ITOP=1` branch), implicitly assuming the two
resonance terms are equal at w=0. That assumption has NO theoretical
basis: only the SUM `term_A(N,N')+term_A(N',N)` is proven symmetric
(Eq. 40 is built as exactly that sum, matching Eq. 30's term-A/term-B
split), never the individual terms.

Because the outer double loop in `secprp.F` genuinely computes
`XPROP(I,J)=2*term_A(I,J)` AND `XPROP(J,I)=2*term_A(J,I)`
independently (confirmed two entries above -- full NOPERT x NOPERT
loop, no shortcut), the raw pre-`Symmet2` mismatch is exactly
`2*(term_A(I,J)-term_A(J,I))` -- a real, physically meaningful gap
between the two resonance terms for that specific system, not noise,
not a convergence artifact, not a bug. `Symmet2`'s averaging,
`0.5*(XPROP(I,J)+XPROP(J,I)) = term_A(I,J)+term_A(J,I)`, is therefore
not defensive cleanup -- **it supplies the missing second resonance
term**, obtained "for free" from the mirror loop position instead of
being computed twice per pair. Deliberate efficiency trick, correctly
implemented.

Consistent with every prior data point without exception: H2/STO-3G
(small, near-symmetric system) has the two resonance terms nearly equal
already -> near-zero raw pre-`Symmet2` gap in earlier tests; distorted
H2O has them genuinely differ more -> the observed 7.6% gap; the
quadratic correction term (established two entries above to have NO
swap-pair structural protection at all) shows an even larger 34.5% raw
gap, consistent with having no such "assume equal, then correct"
mechanism backing it; CALCQUAD1's post-`Symmet2` agreement with NWChem
to 4-5 sig figs makes sense because post-`Symmet2` IS the complete,
mathematically correct formula.

**Verdict**: [[feedback_eomcc_lr_eigenvectors_never_equal]] does NOT
apply here after all -- reinstating it two entries above was premature;
this is the "expected/accepted behavior" case, not a missing-symmetrize
bug, but for a different and more precise reason than blind trust in
the raw framing would have given ("it's just how the approximation
works" is right, but specifically because of the ITOP=1/FACT=2
resonance-term-doubling shortcut, not because of vague "sum-over-states
structural" hand-waving). No code changes needed in `secprp.F`,
`OPRBAR`, `LINEQY`, or `Symmet2`. CLOSED -- clear to move on to
scrnc/GW-work proper.

---

## SESSION 2026-08-17 -- RE-OPENED (needs Ajith's read): the landmark paper's own derivation says quadratic AND CI-like are both exactly symmetric at CCSD truncation, not just full CI

Read the actual text of the landmark paper (now in hand via Dropbox,
`spnspn_jcp.pdf` = Perera, Nooijen, Bartlett, JCP 104, 3290 (1996)) to
verify claims made in this log about symmetry. Eqs. (26) and (30)
derive the quadratic EOM-CCSD model as an exact energy second
derivative `E_NN' = d^2E/dk_N dk_N'`, and state it "is symmetric in N
and N', as required by elementary considerations" -- and explicitly
that this holds "if the excitation space is restricted to singles and
doubles substitution operators," i.e. at CCSD truncation, NOT only
for full CI. Splitting Eq. (30), the CI-like piece (first two terms)
is separately, algebraically symmetric under N<->N' exchange on its
own (it has the form f(N,N')+f(N',N)), independent of the derivative
argument.

This is exactly the claim the PySCF pilot's docstring made ("LR
construction guaranteed exactly symmetric by Schwarz's theorem
regardless of CC truncation") that the entry below (SESSION 2026-08-16
night) retracted as "the wrong standard." Reading the primary source,
that retraction looks like an overcorrection: exact symmetry requires
only that the perturbed-amplitude (response) equations are solved
consistently/converged, not a full-CI/variational wavefunction. Under
this reading, the raw pre-Symmet2 asymmetry seen numerically in ACES2
(CI-like ~7.6%, CALCQUAD1 ~34.5%, see two entries below) is a
numerical/implementation artifact of the ACES2 code path, not expected
structural behavior of either approximation -- which would instead
re-confirm [[feedback_eomcc_lr_eigenvectors_never_equal]] (the earlier
"missing Symmet2 step" finding) over the "corrected" framing below.

**Not yet done**: this needs Ajith's own read/judgment before treating
it as settled -- it reopens a question that flip-flopped multiple
times already in this log (RESOLVED -> CORRECTED -> RESOLVED again
below). The CALCQUAD1-vs-NWChem post-Symmet2 agreement (4-5 sig figs,
see two entries below) is unaffected either way and stands as the
strongest evidence CALCQUAD1 itself is correctly implemented; what's
back open is only the *interpretation* of the raw pre-Symmet2 split.

---

## SESSION 2026-08-16 (night, continued) -- RESOLVED: CALCQUAD1 is NOT buggy, it's validated against NWChem to 4-5 sig figs; CI-like's larger raw symmetry is a "sum-over-states" structural artifact, not evidence of correctness (RESUME HERE)

Direct test of Ajith's question ("does NWChem's lower triangle agree
with ACES2 quadratic's lower triangle?"): since NWChem and ACES2 use
different molecular orientations, compared full-tensor EIGENVALUES
(rotation-invariant) instead of raw axis-labeled elements. Using
NWChem's one genuinely-computed off-diagonal element (`alpha(Z,X) =
1.6268696`, not the mirrored copy -- see the correction two entries
above) plus its three genuine diagonal elements, vs ACES2 CALCQUAD1's
post-Symmet2 tensor from the `vcceh_run.out` run:

| | lambda1 | lambda2 | lambda3 |
|---|---|---|---|
| NWChem TCE-LR (genuine elements only) | 5.593206 | 2.338374 | 0.052156 |
| ACES2 CALCQUAD1 (post-Symmet2) | 5.593252 | 2.338436 | 0.052150 |

Agreement to 4-5 significant figures on ALL THREE eigenvalues -- not
just the isotropic trace (already known, ~0.001% agreement), but the
full anisotropic shape of the tensor. CI-like's isotropic value alone
was off by ~9.4% from NWChem by comparison.

**This reverses the standing conclusion from the entries below.** The
finding that CALCQUAD1's raw pre-Symmet2 left/right split (-0.073442
vs -0.103871, ~34.5%) is LARGER than CI-like's raw split (~7.6%) had
been read as "CALCQUAD1 makes things worse -> implementation bug."
That reasoning implicitly assumed raw pre-symmetrization asymmetry is
itself a correctness metric. It isn't, and Ajith confirmed why: **CI-
like is a sum-over-states-type construction, which is structurally
MORE symmetric by its own form, independent of whether the underlying
physics/numbers are more accurate.** CALCQUAD1's larger raw split is
just what restoring the missing quadratic/size-extensivity term (the
Perera-Nooijen correction) looks like before the two sides are
averaged -- it does not indicate a defect. What actually matters is
the POST-Symmet2 answer, and CALCQUAD1's post-Symmet2 answer matches
an independent, trusted code (NWChem TCE-LR) almost exactly, while
CI-like's does not (9.4% isotropic error). **Conclusion: CALCQUAD1 is
correctly implemented and validated.** The real story of this whole
multi-session investigation was never "find the bug in CALCQUAD1" --
it was confirming that CI-like is a known, expected, structurally-
more-symmetric-but-less-accurate approximation, exactly as the
Perera-Nooijen paper already lays out.

**Reference for the manuscript / for anyone picking this up cold**:
Ajith says the Perera-Nooijen paper has extensive discussion of
exactly how CI-like / LINEAR / QUADRATIC fit together (including the
sum-over-states symmetry point above) -- that paper, not this debug
log's earlier speculative reasoning, is the authoritative source for
the theory-background section. Full citation not yet pulled into this
log; do that before drafting the manuscript's theory section (see
also the Stanton & Bartlett / Koch & Jorgensen lineage noted in the
entry below).

**Open item (2) resolved, same session**: Ajith's call is that
CI-like is the RIGHT construction for scrnc/GW specifically, not a
compromise -- GW itself is built on variational/RPA-type (sum-over-
states) machinery, so CI-like's sum-over-states structural symmetry
(the same property that makes it "more symmetric but less exact" for
a generic CC property, per the entry above) is actually a MATCH for
GW's own theoretical structure, not a mismatch. scrnc should keep
using CI-like (IEOMPROP=0) -- no change needed there. The
quadratic/CALCQUAD1 validation above still matters for the manuscript
(establishes CI-like's known, quantified deviation from the "exact"
CC-LR answer) but is not a reason to switch scrnc's own construction.

Ajith also proposed a specific alternative construction worth noting
for the manuscript's formalism section: instead of the standard CC
similarity-transformed Hamiltonian `H-bar = e^{-T} H e^{T}` (the basis
of the EOM L/R-vector "CI-like" approach), consider `e^{Lambda} H
e^{T}` -- i.e., replacing the EOM left-eigenvector L with the
Lambda-deexcitation operator directly. He noted this may already be
what the PySCF pilot's `mu_relaxed()` function computes (solves
Lambda self-consistently via `mycc.solve_lambda()`, then contracts the
Lambda-relaxed one-particle density via `mycc.make_rdm1()`) -- that
reading looks right: `make_rdm1(t1,t2,l1,l2)` is exactly the
`<0|(1+Lambda) e^{-T} p^dagger q e^{T}|0>`-type density, which is the
Lambda-based (not L-eigenvector-based) construction. Worth revisiting
this correspondence explicitly when drafting the formalism section --
not yet fully worked out here, just flagged.

**Not yet done / open**: ~~(1) pull exact citations (Perera & Nooijen;
Stanton & Bartlett; Koch & Jorgensen) with year/journal for the
manuscript.~~ RESOLVED 2026-08-17 via Dropbox search, see below. (2)
revert the still-uncommitted `secprp.F` debug-print flip
(`if(.false.)` -> `.true.`) and rebuild `xvcceh` clean, per the
not-yet-done note two entries below (still applies).

**Citations pulled 2026-08-17** (Dropbox, `/Private/Documents/Downloads/`):
- **THE landmark paper (confirmed by header scan, supersedes the
  Rozyczko entry below as the primary "Perera & Nooijen" citation):
  Perera, S. A.; Nooijen, M.; Bartlett, R. J.** "Electron correlation
  effects on the theoretical calculation of nuclear magnetic resonance
  spin-spin coupling constants." *J. Chem. Phys.* **104**, 3290
  (1996). File: `spnspn_jcp.pdf`. This paper derives the EOM-CCSD
  CI-like, quadratic, and linked-diagram forms for general
  second-order properties -- exactly scrnc's CI-like/LINEAR/QUADRATIC
  framework -- applied to NMR spin-spin coupling constants. This is
  Ajith's own 1996 paper and the direct theoretical ancestor of scrnc
  (see [[project_scrnc_gw_manuscript_goal]]: scrnc is "a direct
  linear-response construction, generalized from Ajith's own 1994 NMR
  spin-spin-coupling response theory"). Its 1994 predecessor is also
  in Dropbox: **Perera, S. A.; Sekino, H.; Bartlett, R. J.**
  "Coupled-cluster calculations of indirect nuclear coupling
  constants: The importance of non-Fermi contact contributions." *J.
  Chem. Phys.* **101**, 2186 (1994). File: `spnspn_first_jcp.pdf`.
- **Rozyczko, P. B.; Perera, S. A.; Nooijen, M.; Bartlett, R. J.**
  "Correlated calculations of molecular dynamic polarizabilities."
  *J. Chem. Phys.* **107** (17), 6736 (1 November 1997). Secondary
  citation -- same CI-like/LINEAR/QUADRATIC framework applied to
  dynamic polarizabilities instead of NMR couplings. File:
  `Rozyczko_pol.pdf`.
- **Stanton, J. F.; Bartlett, R. J.** "A coupled-cluster based
  effective Hamiltonian method for dynamic electric polarizabilities."
  *J. Chem. Phys.* **99** (7), 5178 (1 October 1993). The original
  CI-like/CC-EH construction paper. File: `Stanton_Bartlett_EOMPOL.pdf`.
- **Koch, H.; Jørgensen, P.** "Coupled Cluster Response Functions."
  *J. Chem. Phys.* **93**, 3333 (1990). The foundational CC-LR paper.
  Not present as a standalone PDF in the Dropbox `Downloads` folder,
  but the exact citation is now confirmed independently from the
  reference lists of two related papers that *are* in Dropbox (both
  under `/Private/Documents/Downloads/Excited states/`):
  `Chris_koch_jorgensen_cc2.pdf` (Christiansen, Koch, Jørgensen,
  "The second-order approximate coupled cluster singles and doubles
  model CC2," *Chem. Phys. Lett.* **243**, 409 (1995)) and
  `Jorgensen_EE2.pdf` (Christiansen, Koch, Jørgensen, "Perturbative
  triple excitation corrections to CCSD excitation energies," *J.
  Chem. Phys.* **105**, 1451 (1996)). Either of these two is usable
  as a secondary citation for the CC-response lineage if the primary
  1990 paper's PDF is never tracked down.

---

## SESSION 2026-08-16 (night) -- CORRECTION (Ajith): the whole "true tensor is exactly symmetric" premise was wrong; only full CI/variational methods guarantee exact symmetry -- CALCQUAD1's regression (7.6% -> 34.5%) is still the real bug (SUPERSEDED, see entry above)

Ajith corrected the framing of every entry below that treated "exact
symmetry" as the ground-truth target: **CI-LIKE and LINEAR
(IEOMPROP=0/1) are both approximations to QUADRATIC (IEOMPROP=2);
"LR" is later-generation terminology, and it was subsequently proven
that ACES2's QUADRATIC construction *is* LR.** CI-like can even look
*more* symmetric than quadratic by construction, coincidentally, and
that means nothing about correctness. **Only full CI or a genuinely
variational method gives an exactly symmetric response tensor** --
truncated CC-LR (quadratic), even implemented perfectly, is not
expected to be exactly symmetric either, just much closer than the
cruder CI-like/LINEAR approximations it improves on.

This retracts the framing (not just the NWChem numbers, see the entry
directly below) of several earlier entries this session: the PySCF
`ci_like_vs_lr_polarizability.py` pilot's docstring claim that the
Lambda-relaxed "LR" construction is *guaranteed exactly symmetric* by
Schwarz's theorem regardless of CC truncation level is not the right
standard to hold ACES2's CALCQUAD1 to. Do not cite that pilot, or the
now-retracted NWChem numbers, as proof that "the true answer is
exactly symmetric" in the manuscript or in future debugging framing.

**What survives this correction, unchanged**: the concrete,
tractable bug. CALCQUAD1 was found to make the alpha_xy/alpha_yx
asymmetry go from ~7.6% (CI-like, IEOMPROP=0) to ~34.5%
(IEOMPROP=2/QUADRATIC) on the same distorted-H2O dipole case -- a
~5x regression, reproduced independently on two separate test cases
(the ACES2-native dipole case and scrnc's own AO-pair case, see
entries below). Since QUADRATIC is supposed to be the *better*
approximation (closer to LR/closer to full-CI-level symmetry, per
Ajith's correction above), it making things dramatically worse is
still strong evidence of a genuine implementation bug in
`CALCQUAD1` -- that conclusion does not depend on "the true tensor is
exactly symmetric" being true. A truncation-level residual asymmetry
would be expected to shrink or stay flat going CI-like -> QUADRATIC,
not grow 5x.

**Resume point, concretely**: proceed to the already-identified next
step -- read `vcceh/calcquad1.F` in full and trace why it produces a
LARGER off-diagonal gap than the cruder CI-like construction it's
supposed to improve on, using the saved `vcceh_run.out` matrices
(both CI-like and CALCQUAD1, full off-diagonal + diagonal) as
reference data. The isotropic/diagonal cross-check against NWChem
(CALCQUAD1 trace matches NWChem to ~0.001%, CI-like trace off by
~9.4%, see two entries below) is unaffected by this correction and
still stands as independently-computed, valid evidence that
CALCQUAD1's diagonal/trace physics is basically right -- it's
specifically the off-diagonal construction that's broken.

**Deeper why, from Ajith, this session -- the real defect in CI-like
is size-extensivity, not tensor symmetry per se**: symmetry was
always a downstream SYMPTOM, not the actual defining flaw. The
EOM-CC property formulation built from L and R vectors (Stanton &
Bartlett's EOM-CC framework, "CI-like" in ACES2's naming) is not
size-extensive -- a strict requirement of proper CC theory. Perera &
Nooijen identified that the L/R-vector EOM construction is MISSING A
QUADRATIC TERM, and derived the correction for it (this is the origin
of `CALCQUAD1`/IEOMPROP=2, "the Perera-Nooijen quadratic
correction"). Koch & Jorgensen independently reached an equivalent
construction from a completely different starting point -- the CC
response-theory/Lagrangian viewpoint, which is where "LR" terminology
comes from -- confirming from the outside that ACES2's quadratic
correction and later-generation "linear response" theory are the same
underlying physics. **Citations to use for the manuscript's
theory-background section**: Stanton & Bartlett (EOM-CC formulation),
Perera & Nooijen (quadratic correction / size-extensivity fix), Koch &
Jorgensen (independent CC-response-theory derivation of the same
result) -- exact citation details (year/journal) not yet pulled, do
that before drafting that section.

---

## SESSION 2026-08-16 (evening, continued yet further) -- CORRECTION (Ajith): NWChem's "exact off-diagonal symmetry" was a code artifact, NOT independent evidence -- retracted; diagonal/isotropic finding still stands

Ajith pushed back on the "NWChem independently confirms the true
CCSD-LR tensor is exactly symmetric" claim two entries below: "I do
not think off diagonal wrong. NWchem might have symmetrized the
tensor to match with physics." Checked the actual NWChem TCE source
(`/dist/apps/nwchem/7.2.3/src/tce/ccsd_lr_alpha.F`) instead of taking
the printed output at face value -- he was right, I was wrong.

The main compute loop (`ccsd_lr_alpha.F:387-388`):
```fortran
do axisA = 1, 3
  do axisB = 1, axisA        ! only axisB <= axisA ever computed!
```
NWChem's TCE-LR module only ever computes ONE triangle of the tensor
(lower triangle including diagonal). It NEVER independently evaluates
the `axisB > axisA` elements. The "full" matrix used for
diagonalization is then built by literally mirroring
(`ccsd_lr_alpha.F:519-524`):
```fortran
do i = 1, 3
  do j = 1, i
    alphacopy(i,j)=alpha(i,j)
    alphacopy(j,i)=alpha(i,j)   ! same value copied into the transpose slot
  enddo
enddo
```
and the print statement (lines 541-546) reuses the single stored
`alpha(2,1)` value for BOTH printed off-diagonal positions. So what
the entry below reported as "alpha_XZ = alpha_ZX to 7 digits,
independent confirmation of symmetry" was NWChem printing the exact
same stored number twice -- not an independent numerical check.
This is a legitimate performance optimization on NWChem's part
(exploiting that a properly Lagrangian-relaxed CC linear-response
tensor IS provably symmetric via Schwarz's theorem on mixed partial
derivatives of a stationary functional -- see the PySCF pilot
docstring, `scrnc/independent_pilot/ci_like_vs_lr_polarizability.py`),
but it means NWChem's off-diagonal numbers carry ZERO independent
evidentiary weight for the ACES2 asymmetry question. **Retracted**:
"THIRD independent line of evidence... exactly symmetric" (entry
below). The only genuine independent evidence that a fully-relaxed
CC-LR construction gives symmetric alpha_xy/alpha_yx remains the
earlier PySCF pilot script, which actually re-solves Lambda
independently at each perturbed field (no symmetry assumed anywhere
in that code).

**What still stands, unaffected by this correction**: the DIAGONAL
elements (`axisA.eq.axisB` branch, line 434-435: `alpha2=alpha1`) ARE
genuinely, independently computed in NWChem -- that branch is a
different, legitimate simplification (a diagonal element is trivially
self-referential: computing "response of mu_i to F_i" only ever needs
one direction, there's no second independent quantity to compare it
against, unlike i!=j where two genuinely different quantities
alpha_ij and alpha_ji COULD in principle disagree and NWChem's code
just never checks). So the isotropic/trace comparison from the entry
below is still valid, independently-computed cross-code evidence:
NWChem isotropic 2.6612454 au vs ACES2 CALCQUAD1 isotropic 2.661279 au
(~0.001% agreement) vs ACES2 CI-like isotropic 2.912252 au (~9.4% off).

**Where this leaves the investigation**: the "CALCQUAD1's diagonal is
right, bug is isolated to the off-diagonal construction" conclusion
(entry below) still holds -- it was derived from the isotropic/trace
comparison, which is unaffected by this correction. But the framing
of NWChem as a third independent code proving the TRUE tensor is
symmetric is gone; that claim now rests on the PySCF pilot and the
Schwarz-theorem argument alone, not on this NWChem run. Do not cite
the NWChem off-diagonal numbers as evidence in the manuscript.

---

## SESSION 2026-08-16 (evening, continued even further) -- CALCQUAD1's ISOTROPIC (trace) polarizability matches NWChem to ~0.001%; CI-like's is off by ~9.4% -- bug narrows to the off-diagonal construction specifically

Ajith asked whether the diagonal elements are identical (pre- vs
post-`Symmet2`, within the CALCQUAD1/quadratic result). Ran ACES2's
own `xvcceh` directly (NOT `xscrnc`) to get real numbers instead of
guessing: flipped the `if(.false.)` debug print at
`aces2/vcceh/secprp.F:474` to `.true.` (temporary, uncommitted, same
as before), rebuilt+installed `xvcceh` (md5-verified fresh), and ran
a plain `xaces2` job in `/home/perera/nwchem_scrnc_crosscheck/
aces2_vcceh_run/` with `SYMMETRY=OFF, BASIS=STO-3G, REF=UHF,
CALC=CCSD, PROPS=EOM_NLO, EOMPROP=QUADRATIC`, same distorted
geometry (R1=0.957, R2=1.05, A=104.5) as the NWChem cross-check
above. `_DEBUG_LVL0` is already on by default in this build
(`GNUmakefile.src`), so the "CI-LIke Pol. matrix" print fired too --
got BOTH full matrices in one run, saved to `vcceh_run.out`:
```
 CI-LIke Pol. matrix
 ROW 1   6.03267129332  -0.01893062086   0.00000000000
 ROW 2  -0.03111777202   2.65062243908   0.00000000000
 ROW 3   0.00000000000   0.00000000000   0.05346127221

 Quadratic contribution, ifreq 1
   1   5.590835  -0.073442   0.000000
   2  -0.103871   2.340853   0.000000
   3   0.000000   0.000000   0.052150

 Molecular Polarizabilities (post-Symmet2)
   1   5.590835  -0.088657   0.000000
   2  -0.088657   2.340853   0.000000
   3   0.000000   0.000000   0.052150
```
Confirms the trivial fact directly: `Symmet2` leaves the diagonal
completely untouched (pre- and post- CALCQUAD1 diagonal identical,
5.590835/2.340853/0.052150 both times) -- it only averages the (1,2)
and (2,1) off-diagonal pair. Diagonal DOES differ between CI-like and
CALCQUAD1 themselves (6.03267/2.65062/0.05346 vs
5.590835/2.340853/0.052150) -- CALCQUAD1 is not a no-op on the
diagonal, it's a real correction there too.

**New, more decisive finding**: the polarizability tensor's trace
(3x isotropic value) is orientation-invariant, so it's a valid
apples-to-apples cross-code check regardless of ACES2 vs NWChem axis
convention differences. Comparing:
- NWChem TCE-CCSD-LR isotropic (its own printed value) = **2.6612454** au
- ACES2 CALCQUAD1 isotropic = (5.590835+2.340853+0.052150)/3 = **2.661279** au -- agrees with NWChem to ~0.001% (diff 0.000034)
- ACES2 CI-like isotropic = (6.03267129332+2.65062243908+0.05346127221)/3 = **2.912252** au -- off from NWChem by ~9.4%

**Implication**: CALCQUAD1's diagonal/trace-level physics is
essentially CORRECT (matches the independent NWChem reference almost
exactly), while CI-like's is not. This narrows the bug: CALCQUAD1 is
NOT globally broken -- its isotropic/diagonal contribution is right,
and the "makes asymmetry worse" failure mode from the entry below is
specifically an OFF-DIAGONAL (cross-term, i != j) construction bug.
This is a much more tractable target for the planned
`vcceh/calcquad1.F` read-through than "the whole routine is wrong."

**Not yet done**: same as before -- read `vcceh/calcquad1.F` in full,
now specifically focused on how it builds/combines the i!=j terms
differently from i==j, using this run's saved `vcceh_run.out` (full
matrices, both CI-like and CALCQUAD1, diagonal AND off-diagonal) as
the reference data. Also revert the `secprp.F` debug-print edit and
rebuild `xvcceh` back to clean before moving on (temporary/uncommitted,
same convention as prior sessions -- do this before ending the
session if not already done).

---

## SESSION 2026-08-16 (evening, continued further) -- NWChem independent cross-check now WORKING; off-diagonal "exact symmetry" claim below is CORRECTED/RETRACTED by the entry above (code only computes one triangle) -- diagonal/isotropic numbers still valid

Set up an independent third-code cross-check in
`/home/perera/nwchem_scrnc_crosscheck/` (NWChem 7.2.3 on HPG, module
`intel/2025.1.0 openmpi/5.0.7 nwchem/7.2.3`): `h2o_ccsd_lr.nw` runs
UHF/STO-3G CCSD static dipole polarizability via
`tce ccsd` + `set tce:lineresp T`, on the exact same distorted H2O
geometry (R1=0.957, R2=1.05 Angstrom, A=104.5) as ACES2's
`h2o_ccsd_polar_check_asym_quad2` case referenced in the entry above.

Hit and fixed FOUR unrelated infrastructure bugs to get it running
(all in `/home/perera/nwchem_scrnc_crosscheck/`, not physics bugs):
1. `NWCHEM_BASIS_LIBRARY` (module-set) pointed at a path that doesn't
   exist on HPG compute nodes (only resolves on login/interactive
   nodes) -- fixed by `cp`ing the real `sto-3g` file from
   `/dist/apps/nwchem/7.2.3/src/basis/libraries/sto-3g` into a local
   `basis_lib/` and exporting `NWCHEM_BASIS_LIBRARY=.../basis_lib/`
   (trailing slash) in the sbatch scripts before `mpirun`. (Passing an
   absolute path directly to `library` in the `.nw` basis block does
   NOT work -- NWChem mangles `/` to `_` and treats the whole string
   as a lookup tag against the compiled-in default dir.)
2. `noorient` is not a real NWChem geometry keyword (valid ones:
   `nocenter`/`noautoz`/`noautosym`) -- NWChem silently misparsed it
   as the geometry's *name* instead of a flag, leaving the
   default-named geometry object empty -> `scf_init: no geometry`.
   Fixed by removing it (`symmetry c1` already disables autosym).
3. `run.sbatch` had `--ntasks=4` with no `--nodes`, so Slurm spread
   the 4 MPI ranks across 3 separate nodes and the job died instantly
   (0-byte output, exit 250) -- fixed with `--nodes=1
   --ntasks-per-node=4` to match the working minimal-test layout.
4. The `tce` block's `2eorb` keyword requires an RHF/ROHF reference
   (`tce_init: 2eorb needs RHF/ROHF ref.`) and is incompatible with
   the UHF reference this case needs -- removed it.

**Result, decisive**: with all four fixed, the job runs clean end to
end (UHF SCF converges to -74.959514194573 Ha, CCSD-LR converges in
20 iterations). Because the `.nw` geometry was given as literal
Cartesian coordinates with `nocenter noautoz` (no reorientation), the
molecule sits in NWChem's XZ-plane, making Y the out-of-plane axis --
so the XY/YX couplings are trivially zero by molecular symmetry and
NOT the informative pair (this is why the printed tensor shows
Y-row/Y-column all zero). The physically meaningful off-diagonal pair
for this orientation is XZ/ZX, and NWChem gives:
```
 X            3.9236260     -0.0000000      1.6268696
 Y           -0.0000000      0.0521562     -0.0000000
 Z            1.6268696     -0.0000000      4.0079540
```
i.e. **alpha_XZ = alpha_ZX = 1.6268696 bohr^3, identical to all 7
displayed digits** -- exactly symmetric, no fix-up/Symmet2-equivalent
step needed anywhere in NWChem's TCE-LR pipeline.

Compare to ACES2 on the equivalent geometry (previous entries):
- IEOMPROP=0 (CI-like): alpha_xy=-0.01893, alpha_yx=-0.03112 (~49% relative asymmetry)
- IEOMPROP=2/CALCQUAD1: alpha_xy=-0.073442, alpha_yx=-0.103871 (~34.5% asymmetry, WORSE)
- NWChem TCE-CCSD-LR (independent code, this entry): exactly symmetric

This is now a THIRD independent line of evidence (after the earlier
PySCF pilot, which confirmed the exact CCSD response function is
symmetric) that the true CCSD linear-response polarizability tensor
IS symmetric, and that BOTH of ACES2's schemes (baseline CI-like
*and* CALCQUAD1) have a genuine implementation bug in ACES2's own
oprbar/Lineqy/CALCQUAD1 machinery -- not an accepted/expected
characteristic of the CI-like formalism. The open question from the
entry above ("is the pre-Symmet2 asymmetry already-known/accepted
physics?") is now answered: no, it is not accepted physics, it is a
bug, and a real one worth chasing down for the manuscript.

**Not yet done**: run the h2o_ccsd_lr.nw case a second time with the
geometry reoriented (or added noautoz/nocenter removed) so the
in-plane axes land on X/Y instead of X/Z, to get a literal
alpha_xy/alpha_yx pair directly comparable to ACES2's own printed
labels (cosmetic only, not expected to change the conclusion). Then
return to the original plan: read `vcceh/calcquad1.F` in full and
start the backward-trace/instrumented-checksum methodology against
this now-working NWChem reference.

---

## SESSION 2026-08-16 (evening, continued) -- CORRECTION: EOMPROP=QUADRATIC IS reachable (Ajith); Calcquad1 itself confirmed buggy, makes asymmetry WORSE not better (RESUME HERE)

Ajith corrected the previous entry's "IEOMPROP=2 is unreachable"
conclusion: "EOM_PROPS=QUADRATIC should do CCLR." He was right, and
the earlier conclusion was wrong -- the error was searching only
Fortran source for an assignment to `IFLAGS2(106)`/`h_IFLAGS2_eomprop`
and finding none. The real keyword-to-flag path goes through
`joda/asv_ctl.c` (the C-based ASV/keyword parser), which is invisible
to a Fortran-only search. Found it there: `h_IOPPAR_eomprop=206` maps
handle strings `CILIKE=0/LINEAR=1/QUADRATIC=2/COMBO=3`
(`asv_ctl.c:988`), and confirmed via `flags.h`
(`h_IOPPAR_print=1 /* START IFLAGS */`) that "IOPPAR" indices above
100 are the SAME underlying array as `IFLAGS2`, just offset by 100 --
i.e. `IOPPAR` index 206 IS `IFLAGS2(106)` in Fortran-side terms, no
separate array, no missing translation. The real ZMAT keyword is
`EOMPROP=QUADRATIC` (test case `tests/fast/cases/046/ZMAT` already
uses it, `EOMPROP=QUADRATIC` alongside `PROP=JSC_ALL`). Confirmed
working: rerunning the ACES2-native distorted-H2O CCSD polarizability
case (previous entry) with `EOMPROP=QUADRATIC` in the ZMAT shows
`Updating (EOMPROP,QUADRATIC): 'QUADRATIC' converts to 'QUADRATIC' ->
2` and "Fully quadratic evaluation of properties." in the output --
the keyword genuinely reaches `IEOMPROP=2` and `CALCQUAD1` genuinely
runs.

**Decisive result, now that the correction is tested with real data**:
`Calcquad1`'s own debug print (`secprp.F` line ~474, guarded by an
`if(.false.)` that was temporarily flipped to `.true.` to see it, then
reverted -- never committed) shows the matrix immediately AFTER
`CALCQUAD1` runs, still BEFORE `Symmet2`:
```
 Quadratic contribution, ifreq 1
   alpha_xy = -0.073442
   alpha_yx = -0.103871
```
Absolute difference 0.0304, ~34.5% relative to the average magnitude
-- WORSE than both the plain CI-like (IEOMPROP=0) raw asymmetry
(7.6% relative, previous entry) AND the pre-CALCQUAD1 "CI-LIke Pol.
matrix" printed just before it in this SAME run (alpha_xy=-0.01893,
alpha_yx=-0.03112, ~49% relative but smaller absolute gap of 0.0122).
**`CALCQUAD1` makes the asymmetry substantially WORSE, not better**,
reproducing EXACTLY the same qualitative result found earlier today
when `Ieomprop=2` was force-hacked into `scrnc_set_vars.F` for scrnc's
own AO-pair case (asymmetry went from 5.17e-4 to 3.58e-3, ~7x worse)
-- now confirmed on a SECOND, INDEPENDENT, ACES2-native test case
using the real keyword, not a hack. Two independent reproductions of
the same qualitative failure mode is strong evidence this is a genuine
implementation bug in `CALCQUAD1` (or shared infrastructure it calls),
not a physics-concept problem or a scrnc-specific issue.

**Where this leaves the investigation, now on much firmer ground**:
- The BASELINE (`IEOMPROP=0`, "CILIKE") pre-`Symmet2` asymmetry
  (0.14% for scrnc's H2O case, 7.6% for this dipole case) is confirmed
  present in ACES2's own core property code, not scrnc-specific --
  still an open question whether this baseline-level gap itself is a
  known/accepted characteristic of the CI-like scheme or a separate,
  deeper bug (unresolved either way).
- The DOCUMENTED FIX for that gap (`EOMPROP=QUADRATIC`/`CALCQUAD1`) is
  reachable and does run, but is ITSELF demonstrably broken --
  confirmed on two independent test cases -- making things worse
  rather than closing the gap it's meant to fix. This is now the
  clearer, more tractable bug to actually chase: `CALCQUAD1.F` is a
  much smaller, more self-contained routine than the full oprbar/
  Lineqy machinery, and now has TWO independent, reproducible failing
  test cases to debug against (the ACES2-native dipole case is
  simpler/cleaner than scrnc's AO-pair case and doesn't need scrnc at
  all to reproduce).

**Not yet done**: read `vcceh/calcquad1.F` (or wherever `CALCQUAD1` is
implemented) in full and start the same backward-trace/instrumented-
checksum methodology used throughout this whole project, using the
ACES2-native distorted-H2O dipole polarizability case
(`h2o_ccsd_polar_check_asym_quad`/`_quad2` ZMATs, already built and
working) as the primary reproducer -- it's simpler and doesn't require
scrnc's own AO-pair machinery to be correct first.

---

## SESSION 2026-08-16 (evening, continued) -- MAJOR REFRAME: ACES2's own core property machinery has the SAME pre-Symmet2 asymmetry (SUPERSEDED, see entry above -- "IEOMPROP=2 unreachable" was WRONG)

Ajith's direction: "since you have a pyscf reference, check every step
and match, find out what this bug is all about." Started by trying to
match the `Mexpect` convention exactly (see the entry below --
inconclusive, a real spin-block gap remained). Redirected to a more
powerful test: rather than fighting to reproduce scrnc's exact AO-pair
integral/density conventions in Python, checked whether the SHARED
"oprbar" machinery ITSELF (used identically by scrnc for AO-pair
operators and by ACES2's own standard `vcceh/secprp.F` for dipole
properties) has this asymmetry AT ALL, using ACES2 itself as the
reference this time -- no PySCF, no convention-matching risk.

**Method**: ran ACES2's own standard CCSD static dipole polarizability
(`PROPS=EOM_NLO`, triggers `vcceh/secprp.F` -> the exact same
`oprbar1r.F`/`oprbar2u.F`/`Lineqy.F` routines scrnc's
`Scrnc_form_mubar_s.F`/`_d_uhf.F` were adapted from) for H2O/STO-3G/
UHF/CCSD, sym=off. First attempt used the SAME symmetric H2O geometry
as all along -- gave an exactly-diagonal polarizability tensor
(alpha_xy=alpha_yx=0.0), which turned out to be a bad test: H2O's own
C2v symmetry (and/or ACES2's automatic reorientation to a principal-
axis frame) forces the off-diagonal elements to vanish trivially,
telling us nothing about whether they'd AGREE if nonzero. Switched to
a deliberately ASYMMETRIC H2O geometry (R1=0.957, R2=1.05 Angstrom,
breaking C2v) to get genuinely nonzero, independently-computed
alpha_xy and alpha_yx.

**Result, decisive**: ACES2's own raw output, literally labeled
"CI-LIke Pol. matrix" (same "CI-like" naming convention scrnc uses for
its own pre-`Symmet2` matrix!) in the `PROPS=EOM_NLO` pipeline log:
```
 CI-LIke Pol. matrix
 ROW 1  ...   -0.09222300617   ...
 ROW 2  -0.08519804618   ...
```
`alpha_xy=-0.09222300617` vs `alpha_yx=-0.08519804618` -- these
DIFFER by ~7.6% relative, LARGER than scrnc's own H2O AO-pair
discrepancy (~0.14% relative). The final displayed "Molecular
Polarizabilities" table shows `-0.088711` for BOTH (1,2) and (2,1) --
exactly the simple average of the two raw values, i.e. Symmet2 forcibly
symmetrizing them, identical in spirit to what scrnc's own
`Call Symmet2(JJ_lrespn,Nao_pair)` does.

**Conclusion, tentative but significant**: this pre-symmetrization
asymmetry is NOT specific to scrnc's AO-pair adaptation -- it is
apparently a property of the shared "CI-like"/EOM-effective-Hamiltonian
(IEOMPROP=0, "John's original implementation" per the code's own
comments, as distinct from "Hideo's modification" IEOMPROP=1 and the
"quadratic Perera-Nooijen correction" IEOMPROP=2) construction ITSELF,
which ACES2's own core, decades-established property module has
apparently ALWAYS relied on `Symmet2` to paper over, for properties in
general (not just scrnc's screened-Coulomb AO-pair case). This is
CONSISTENT with, not contradicted by, the earlier independent-PySCF-
pilot finding: that pilot confirmed the EXACT, TRUE CCSD linear
response function (via finite-difference of the fully relaxed T/Lambda
density matrix, a different and exact construction) IS symmetric to
high precision -- meaning the CI-like/EOM-effective-Hamiltonian
SCHEME is an APPROXIMATION to the true response function that happens
to not exactly satisfy the reciprocity theorem, with the size of the
violation apparently varying by system/property (0.14% here, 7.6% for
this particular distorted-water polarizability case) -- not a
mechanical coding defect with a single fixed root cause to find and
patch.

**This does NOT mean "nothing to do here"** -- it reframes what the
open question actually is. Possibilities, not yet distinguished:
1. This IS simply how the IEOMPROP=0 CI-like scheme behaves, always
   has, and Symmet2 is the correct, intentional, long-standing
   remedy -- in which case scrnc doing the same thing is not a defect
   at all, just consistent with established ACES2 practice.
2. IEOMPROP=2 ("the quadratic correction... documented fix for a gap
   like this," per the earlier-cited note) is SUPPOSED to reduce or
   eliminate this asymmetry for genuine CCSD calculations, and the
   earlier test showing it made scrnc's OWN asymmetry worse (not
   better) means something in HOW IEOMPROP=2 was invoked was
   incomplete/wrong (recall it's currently unreachable via any ZMAT
   keyword in this codebase -- untested, unvalidated machinery) --
   worth testing IEOMPROP=2 on THIS SAME ACES2-native dipole-
   polarizability case (where it presumably CAN be reached via a
   documented keyword) to see whether it actually shrinks the
   alpha_xy/alpha_yx gap there, which would validate or refute
   hypothesis 2 independent of any scrnc-specific code question.
3. There's a genuine, separate implementation bug that happens to
   affect BOTH the dipole and AO-pair cases identically, since they
   share the same "oprbar" machinery -- possible, but nine-plus prior
   clean-vs-reference checks make this the least likely of the three.

**IEOMPROP=2 reachability checked -- it's dead code throughout the
ENTIRE codebase, not just scrnc**: confirmed `vcceh.F` itself reads
`IEOMPROP = IFLAGS2(106)` (line 95), the exact same flag scrnc reads.
The earlier `grep -rl h_IFLAGS2_eomprop` search (only the read sites in
`scrnc_set_vars.F`/`rpa/scrnc_set_vars.F`) already established nothing
anywhere WRITES this flag -- this now confirms it's not a scrnc-
specific gap, it's that IEOMPROP=2 ("the quadratic Perera-Nooijen
correction," the code's own documented fix for a CI-like-vs-exact gap)
is genuinely UNREACHABLE from any ZMAT keyword for ANY property
calculation in this whole ACES2 build, CCSD or otherwise. This means
`Calcquad1` (the routine that implements it) has essentially never
been exercised in production, in this codebase's current state --
consistent with the earlier finding that force-hacking `Ieomprop=2` in
`scrnc_set_vars.F` made scrnc's asymmetry WORSE rather than better:
that's much more likely to mean `Calcquad1`'s implementation itself is
incomplete/unvalidated (never having had a real caller) than that the
underlying quadratic-correction PHYSICS idea is wrong. Testing
IEOMPROP=2 "the real way" on the ACES2-native dipole case would
require the identical kind of source hack (there is no legitimate
path), so it would not actually be a cleaner test than the one already
done -- not pursued further for now.

**Where this leaves the investigation**: the asymmetry is confirmed
present in ACES2's own core, decades-established CI-like/EOM-effective-
Hamiltonian property machinery (IEOMPROP=0), not just scrnc's
adaptation of it -- and the one available attempt at the documented
"fix" (IEOMPROP=2) is unreachable/unvalidated dead code throughout the
whole build. This is squarely a question for Ajith's own domain
knowledge: is a small pre-Symmet2 asymmetry in CI-like/EOM-effective-
Hamiltonian properties (of the magnitude seen here, ~0.1-8% relative
depending on system/property) ALREADY KNOWN, documented, and accepted
ACES2 behavior for this scheme (in which case scrnc doing the same
thing is not a defect, just consistent with established practice) --
or is this itself news, meaning there may be a genuine, deeper issue
in the shared oprbar/Lineqy machinery that predates scrnc entirely and
has simply never been noticed because Symmet2 always hides it in the
final published numbers.

---

## SESSION 2026-08-16 (evening, continued) -- Mexpect convention reconciliation attempted, partial progress, root cause still OPEN

Attempted step 1 of the two next-steps named in the entry below:
extend the independent pilot to compute `Mexpect` directly (via
PySCF's `mycc.make_rdm1()`, the standard Lambda-CCSD relaxed 1-RDM)
and diff against scrnc's own printed "mu_0/mu_bar" values for AO-pair
9/41 (sym=off run, `h2o_scrnc_ccsd_uhf_sto3g_crosscheck2`: AO-pair 9
Mexpect=-0.021827/Mbar00=-0.008283, AO-pair 41
Mexpect=-0.016931/Mbar00=-0.005255).

**First attempt (raw `Tr[D_relaxed*M]`) way off**: 2.01 and 0.64 vs
scrnc's -0.022 and -0.017 -- wrong order of magnitude entirely (these
look like genuine electron-population-scale density-matrix diagonal
values, not scrnc's much smaller numbers).

**Correlation-only density gets much closer**: subtracting the UHF
reference density (`D_corr = D_CCSD - D_HF`) before contracting with
`M` gives 0.0427 and 0.0411 -- same order of magnitude as scrnc's
values now (consistent with `Mexpect` representing a CORRELATION
correction to a density-matrix element, not the raw expectation value,
which makes physical sense for a "screened Coulomb via CC response"
construction). But the ratio to scrnc's values isn't a clean, constant
factor across the two AO pairs (-1.96 for AO-pair 9, -2.43 for AO-pair
41, a real ~20% spread) -- meaning there's still a genuine convention
detail not yet matched (most likely how the AO-pair "M" operator
should be applied per spin-block: `dm_ao = dm1a+dm1b` sums both spins
uniformly, but scrnc's `Scrnc_i1tran.F`-built `MIJ`/`MAB`/`MAI` blocks
are built and summed per-`ISPIN` INSIDE `Scrnc_form_mubar_00.F`'s own
loop, which may not be equivalent to my "sum densities first, then
contract once" approach if the AO-pair operator itself needs to be
spin-resolved differently, e.g. only counted once total rather than
implicitly doubled by being applied to both alpha and beta density
blocks identically).

**Decision, given the scope already covered this session**: this
reconciliation is turning into its own open-ended sub-investigation
(spin-block convention matching) with uncertain payoff, layered on top
of an already extremely long session (2 confirmed+committed fixes,
independent-pilot confirmation the bug is real, nine-plus clean
code-level checks). Stopping the convention-chase here rather than
continuing to sink effort into an increasingly uncertain path --
flagging this as the concrete resume point rather than pushing further
unprompted.

**On resume**: the most promising concrete next step is fixing the
spin-block handling in the pilot's `Mexpect` calculation (build and
contract `MIJ`/`MAB`/`MAI`-equivalent quantities PER SPIN, matching
`Scrnc_form_mubar_00.F`'s own `DO 10 ISPIN=1,1+IUHF` loop structure
exactly, rather than summing `dm1a+dm1b` before contracting) --
if that closes the remaining ratio gap to a clean constant, the
resulting reconciled pilot could then directly reproduce scrnc's own
`Mbar00`/`Mexpect`/full response numbers element-by-element, which
would be a far sharper localization tool than anything tried so far
this session.

---

## SESSION 2026-08-16 (evening, continued) -- Lineqy invocation + oprbar.f orchestration audited, both clean; root cause STILL open

Continued from the independent-pilot session below (which confirmed
the bug is real, not missed physics). Went back to code-level auditing
with that confirmation in hand, targeting the two leads named at the
end of that entry.

**Lineqy invocation checked against `vcceh/secprp.F`'s own call**:
arguments (`IOFFSET`, `IPERT`, `FREQ`, `ISIDE`, `REALFREQ`) match in
order and role. Found `IPERT` is passed with swapped semantics
relative to secprp.F (`Ioffset=Npert_count`/`Ipert=Iao_pair` in scrnc
vs `IOFFSET=IOFFPRT(IPERT)`/`IPERT`=global counter in secprp) --
harmless: confirmed by direct inspection that neither `lineqy.F` nor
`plineqy.f` ever actually USES the `IPERT` argument anywhere in their
bodies (dead parameter, present in the signature only). For `Nirrep=1`
(our sym=off test), `Npert_count` and `Iao_pair` are numerically
identical anyway, so `IOFFSET` itself is unaffected regardless. (Note
for whenever the sym=on bug is revisited: `Ioffset` should really be a
within-irrep-block position like secprp's `IOFFPRT(IPERT)`, not the
raw global `Npert_count` -- a separate, real latent bug for the
multi-irrep case, though moot for sym=off.)

Checked `Lineqy.F`'s own convergence: hardcoded `CUTOFF=1.D-14`
(machine precision, not a loose tolerance), and non-convergence within
`KMAX` iterations hits `CALL ERREX` (hard error, run terminates) rather
than silently returning a partially-converged result -- so "loose/
incomplete Davidson convergence for the richer H2O subspace" is not a
viable explanation either (a real convergence failure would have
crashed the run, which it didn't -- confirmed via
`grep "Derivative amplitudes converged"` showing clean 7-8 iteration
convergence for every perturbation, no anomalies).

**`vcceh/oprbar.f` (the actual trusted orchestration routine
`secprp.F` calls, which scrnc split into `scrnc_lrspn_cc_main.F` +
`scrnc_jj_cc_lrspn.F`) read in full and traced line-by-line against
scrnc's own Iside=1-then-Iside=2 sequencing** -- this is exactly the
kind of scrnc-specific adaptation logic (not a near-verbatim routine
copy) that could hide a bug invisible to the six previous per-routine
diffs. Specifically verified the "MBAR00-MEXPECT for Iside=2" argument
convention (matches exactly: `oprbar.f` lines 170/184 use
`MBAR00-MEXPECT` for `ISIDE=2`, plain `MBAR00` for `ISIDE=1`, scrnc's
`scrnc_jj_cc_lrspn.F` does the identical thing) and the L-side
array-seeding pattern -- `oprbar.f` seeds its L-side MIJBAR/MABBAR/
MAIBAR/MIABAR OUTPUT arrays fresh from the BARE M(ij)/M(ab)/M(ai)
(`IMLIJBAR` etc., lines 154-158), while passing the ALREADY Iside=1-
DRESSED R-bar arrays (`IMRIJBAR`/`IMRABBAR`/`IMRAIBAR`) as the
formula's own "M" INPUT slot for the Iside=2 call (line 166-170) --
traced scrnc_jj_cc_lrspn.F's own Dcopy/argument-passing pattern
(`Work(I000/I010/I020)` = still-bare Ioff_oo/vv/vo copied fresh into
the Iside=2 output slots `Work(I070/I080/I090/I100)`, while
`Work(I030/I040/I050)` = the Iside=1-dressed values from the earlier
call, passed as the Iside=2 formula's own "MIJ,MAB,MAI" input) and
confirmed it matches this exact pattern precisely.

**Where this leaves things**: essentially the entire response-
construction pipeline has now been audited against trusted references
at BOTH the individual-routine level (six clean diffs: mubar_00,
mubar_s, mubar_d_uhf, load_vec, plus direct empirical checks of raw
MIJ symmetry and Nsizec sizing) AND the orchestration level (Lineqy
invocation, oprbar.f sequencing) -- nine-plus independent checks, all
clean. The independent PySCF pilot (below) already confirms the bug is
real. At this point continued line-by-line code auditing has a very
low remaining hit rate. The two genuinely different angles still worth
trying, in priority order:
1. Reconcile the independent pilot's sign/scale convention with
   scrnc's `Ddot(373,374)*Fact` construction (noted as not-yet-done in
   the pilot entry below) to get an ELEMENT-BY-ELEMENT numerical
   comparison against the raw list-373/374 dumps already captured for
   `Npert_count=9`/`41` -- this would localize the bug far more
   precisely than either the reciprocity check or code-reading alone.
2. Extend the independent pilot itself to compute `Mbar_00`/`Mexpect`
   (or the full Xi-vector construction, now that the exact `oprbar.f`
   formulas have been read and traced carefully enough to reduce the
   re-derivation risk that motivated avoiding this route originally)
   and diff element-by-element against scrnc's own raw checksums/prints
   for the SAME perturbations.

---

## SESSION 2026-08-16 (evening) -- real Abelian symmetry turned on for the first time ever; bug found, one false lead found+reverted, real cause still OPEN

Ajith's direction: "let's get abelian sym as rest of ACES II." Every
single scrnc test case in this entire project (all sessions above) had
used `sym=off`. Built the first-ever real-symmetry scrnc test:
H2O/STO-3G/UHF/CCSD/DENSITY=RESPONSE, C2v (same ZMAT as the trusted
`h2o_scrnc_ccsd_uhf_sto3g` control case, just with `sym=off` removed --
`/blue/ufhpc/perera/scrnc_test/h2o_scrnc_ccsd_uhf_sto3g_symon`). ACES2
auto-detects C2v correctly (`Iaopop(1..4)=(4,1,2,0)`,
`Irpdpdao(1..4)=(21,8,16,4)`, summing to Nao=7/Nao_pairs=49). Pipeline
runs to completion, no crash.

**Real bug found, immediately**: parsed the pre-`Symmet2` "Response-only
JJ matrix" with a small Python script (`file` reports this pipeline.log
as "data" not "ASCII text" -- some stray non-UTF8 byte trips plain
`grep`'s binary-file heuristic silently; `grep -a` or Python's
`encoding='latin-1'` are needed, NOT a race condition -- this explains
several earlier "grep returned nothing" false alarms this session,
worth remembering permanently). Found 203/210 nonzero off-diagonal
pairs asymmetric (~1e-4 to 4e-4), all confined to AO-pair range 1-21
(the totally-symmetric A1 irrep block) -- an internal ordering/
consistency issue *within* that block, invisible on every prior
`sym=off` (Nirrep=1) run since there was only ever one such block.

**Ruled out, empirically**: `scrnc_prep_2ints.F` /
`Scrnc_reord_cols.F` / `Scrnc_sympack_aoints.F` (the AO-pair
irrep-grouping machinery, and the `Irrep_ao_pairs` array both CCSD and
ring-CC drivers rely on for `Npert_count`/`Irrepr` bookkeeping) --
isolated the BARE (uncorrelated) AO integral part by subtracting
"Response-only" from "CI-like" (=bare+response); the bare part is
symmetric to ~1e-11 (parser rounding noise), 100% of the real asymmetry
is in the response part alone. A plausible-looking hand-derived theory
about a row/column reordering mismatch in `Scrnc_reord_cols.F` (kept
in git history, not repeated here) turned out to be wrong once checked
this way -- first instance this session of "verify empirically before
trusting a plausible reading," see the second instance right below.

**FALSE LEAD, found then reverted** -- worth keeping the story, not
just the correction, per this log's own established practice:
1. Found `scrnc_form_mubar_d_uhf.F`'s `IRREPX.EQ.1 .AND. ISIDE.EQ.2
   .AND. IEOMPROP.EQ.0` correction (inside `DO 300 IRREPZR=1,NIRREP`)
   lacked an `IRREPZR.EQ.IRREPX` gate that a structurally similar
   correction just above it (lines 398-405) has. Patched it to add
   that gate; committed as `03dceae1`.
2. **This was WRONG.** Found the actual origin of this whole routine --
   `vcceh/oprbar2u.F`, the general, decades-old ACES2 EOM-CC property
   routine (`grep -rl OPRBAR2U`) this was clearly adapted from
   (`INSMEM('OPRBAR2U',...)` literal string still present in scrnc's
   own error messages). `oprbar2u.F` line 440 has the EXACT SAME
   `IF(IRREPX.EQ.1 ...)` with NO `IRREPZR.EQ.IRREPX` gate, in the
   identical loop position -- the trusted, presumably symmetry-tested-
   for-decades reference does it the same "ungated" way scrnc did.
   Re-examined why: the `SAXPY` input (`ICORE(I010)`) is freshly read
   via `GETLST(...,IRREPTR,LISTT)` where `IRREPTR=DIRPRD(IRREPZR,
   IRREPX)` -- a DIFFERENT, irrep-appropriate T1 sub-block on every
   `IRREPZR` iteration, not the same value repeated. This term
   legitimately distributes a correction across every irrep block of
   the 2-particle quantity; the ADJACENT `IRREPZR.EQ.IRREPX`-gated term
   (a genuinely different, LOCALIZED correction using a fixed `MIABAR`
   offset) was not a valid analogy for it. Treating "structurally
   similar-looking code two blocks apart" as "must follow the same
   gating rule" was the error.
3. Reverted cleanly (`git revert --no-edit 03dceae1`, commit
   `3635808b`) rather than a destructive reset, since `03dceae1` was
   already on `master`.

**Systematic reference comparison done instead** (the right way to
have started, and what to do FIRST next time a scrnc routine looks
suspect): `scrnc_form_mubar_d_uhf.F`, `scrnc_form_mubar_s.F`, and
`scrnc_form_mubar_00.F` all have clean, findable vcceh ancestors --
`oprbar2u.F`, `oprbar1r.F`, `oprbar0.f` respectively (1:1
correspondence confirmed via `wc -l` + `diff -u`). All three diffs are
CLEAN: only cosmetic changes (renamed subroutine/variables, removed
antisymmetric-perturbation handling scrnc never needs since its
perturbations are always symmetric Coulomb-type, added debug
checksums, and scrnc's own genuinely-new `PUTLST`-to-lists-480/482/456
persistence). No formula or indexing bug found in any of the three.

**Also checked, directly, empirically**: added a targeted diagnostic in
`scrnc_lrspn_cc_main.F` (guarded `_DEBUG_LVL0`, fires only for
`Npert_count.EQ.9`, one of the worst-offending AO pairs in the
sym=off H2O regression -- see the separate entry below) dumping every
`MIJ(i,j)` vs `MIJ(j,i)` element of the raw alpha OO block (input to
`Scrnc_form_mubar_s`/`_00`/`_d_uhf`, built by
`Scrnc_prep_culomb_ints.F` from `Scrnc_i1tran.F`'s AO2MO transform).
Result: exactly symmetric to machine precision (~1e-16) for all 10
off-diagonal pairs of the 5x5 alpha OO block. Rules out the raw
integral/AO2MO prep chain too.

**Where this leaves things**: the three "mubar" formula routines are
clean (verified against trusted references) and their raw input is
clean (verified directly). `JJ(p,q)=JJ(q,p)` is fundamentally a
CROSS-perturbation identity (`Ddot(373(p),374(q))` must equal
`Ddot(373(q),374(p))`), not a property of any single perturbation's own
R-side/L-side construction in isolation -- every check done so far
(Mexpect/Mbar00 R==L, mubar_d_uhf R==L, MIJ symmetry) is inherently
single-perturbation-scoped and structurally CANNOT catch a bug that
only shows up comparing two DIFFERENT perturbations. This is the right
next angle on resume: instrument list 373/374 raw values for a specific
pair (e.g. Npert_count=9 and 41) and directly verify or refute the
cross identity by hand.

**Committed**: `scrnc_form_mubar_d_uhf.F` (revert `3635808b`), this log.
**Left in the tree, uncommitted**: the `Npert_count.EQ.9` MIJ-symmetry
diagnostic in `scrnc_lrspn_cc_main.F` -- harmless
(`_DEBUG_LVL0`-guarded, single-perturbation-scoped), kept as a template
for the next targeted check.

---

## SESSION 2026-08-16 (evening) -- SEPARATE bug: H2O/STO-3G/CCSD asymmetric even with sym=off, NOT a convergence artifact (OPEN, bigger than the symmetry bug above)

While regression-testing the (now-reverted) IRREPZR change above,
parsed the `h2o_scrnc_ccsd_uhf_sto3g` control case's own "Response-only
JJ matrix" (the SAME system the 2026-08-15 evening session's RESUME
POINT claimed was "exactly symmetric, max diff = 0.0") with the same
rigorous programmatic parser used for the sym=on case. **It is not
symmetric**: 600 of 1176 nonzero off-diagonal pairs differ, up to
5.17e-4 (JJ(9,41)=0.36434921829 vs JJ(41,9)=0.36486641546, both
independently verified by hand against the raw printed matrix, not a
parser artifact). The worst-offending pairs disproportionately involve
AO-pair index 9 (=natural AO pair (2,2), O's 2s orbital self-pair) and
41 (=natural AO pair (6,6), one H's 1s self-pair) -- both DIAGONAL
self-pairs. This is almost certainly why every prior "exactly
symmetric" verification in this whole project's history never caught
it: EVERY SINGLE ONE used H2/STO-3G (1 occupied orbital per spin -- the
OO block is a trivial 1x1 scalar, so an OO-block off-diagonal bug is
structurally invisible there). H2O has 5 occupied orbitals per spin --
a genuinely non-trivial 5x5 OO block -- and is apparently the first
time this project has checked JJ symmetry for a system where that
block isn't degenerate.

**Ruled out, decisively**: numerical convergence noise. Reran with
`cc_conv=12,scf_conv=12` (vs the control's `cc_conv=9,scf_conv=10`) --
the asymmetry is bit-identical to the last printed digit (5.172e-04
both times). A real convergence-noise effect would shrink when
tightening convergence by 3+ orders of magnitude; this didn't move at
all, meaning it's a deterministic, reproducible bug, not numerical
residue.

**Ruled out**: the same R-side/L-side "computed twice" double-build
check used (initially, wrongly, see above) to chase the IRREPZR false
lead -- `Scrnc_form_mubar_00` Mexpect/Mbar00 and
`Scrnc_form_mubar_d_uhf` Iside=1 checksums both agree exactly between
R-side and L-side's internal redundant recompute for this sym=off H2O
case too. Whatever is wrong is NOT the "same quantity computed twice
inconsistently" pattern that explains most of this project's prior
bugs -- it's something that's WRONG THE SAME WAY both times it's
computed (invisible to an R==L cross-check on a single perturbation).

**Confirmed cross-path**: reran the SAME H2O/STO-3G system through the
ring-CC/DRCCD path instead of CCSD (fresh ZMAT,
`/blue/ufhpc/perera/scrnc_test/h2o_scrnc_drccd_uhf_sto3g_check`) --
ALSO asymmetric, even worse (1009/1176 pairs), with the SAME
disproportionate AO-pair-index-9 signature. CCSD and ring-CC share
`Scrnc_form_mubar_s.F`/`Scrnc_form_mubar_00.F`/
`Scrnc_prep_culomb_ints.F` but NOT `Scrnc_form_mubar_d_uhf.F`
(CCSD-only) or ring-CC's own Hbar-based doubles machinery -- since
BOTH paths show the same signature, the shared root cause is most
likely in the shared code, NOT in CCSD-specific doubles handling. (This
prediction held up: the systematic reference-diff work in the entry
above found `Scrnc_form_mubar_s.F`/`_00.F` clean too, and the direct
MIJ-symmetry check found the shared `Scrnc_prep_culomb_ints.F`/AO2MO
input clean as well -- so the shared-code hypothesis is confirmed
correct in scope, even though the specific culprit inside that shared
scope is still unidentified.)

**Why this matters more than the symmetry bug above**: this affects the
CCSD AND ring-CC paths' correctness for ANY system with a non-trivial
OO or VV block, independent of point-group symmetry -- i.e. it's not
really a "symmetry" bug at all, despite being found while chasing one.
Turning on Abelian symmetry (this session's actual goal) can't be
fully validated as "done" until this is fixed too, since even the
sym=off H2O baseline is now known-broken.

**Cross-perturbation identity check done** (2026-08-16, continued):
added a targeted diagnostic in `scrnc_jj_cc_lrspn.F` (guarded
`_DEBUG_LVL0`, fires only for the `(Iao_pair,Jao_pair)` pairs
`(9,41)`/`(41,9)`) splitting the final `Ddot(373(p),374(q))` into its
VO part (first 20 elements, the one-particle/singles-like `MAIBAR`
block) and doubles part (remaining 120 elements), plus raw-value and
checksum dumps of 373/374 themselves. Result: `X(9,41)=0.36434921829`
vs `X(41,9)=0.36486641546` (matches the matrix-parsed value exactly).
Decomposition: `Xvo(9,41)=0.36229369334`, `Xvo(41,9)=0.36370391709`
(diff +0.00141) and `Xdbl(9,41)=0.00205552495`,
`Xdbl(41,9)=0.00116249837` (diff -0.00089303) -- the two partially
cancel to give the observed +0.00051720 total. Both the VO
(one-particle) part AND the doubles part are individually MORE
asymmetric than the final combined total -- the bug is not confined to
the doubles machinery, and is present even in the simplest (VO-only)
piece of the construction.

**IEOMPROP=2 hypothesis tested and REFUTED**: the earlier
"[[project_scrnc_gw_method]]"-era note that "the EOM-effective-
Hamiltonian/CI-like route is known to differ from genuine linear
response... IEOMPROP=2 (Calcquad1) is the documented fix" looked like
a strong candidate explanation, especially since `grep -rl
h_IFLAGS2_eomprop` finds ONLY the read site in
`scrnc_set_vars.F` -- IFLAGS2(106) is never written anywhere in this
codebase, so IEOMPROP=2 is currently UNREACHABLE from any ZMAT
keyword. Tested directly: temporarily hardcoded `Ieomprop=2` in
`scrnc_set_vars.F` (clearly marked as a throwaway diagnostic, reverted
immediately after, never committed), rebuilt, reran the same H2O
cross-check. Result: asymmetry got WORSE, not better --
`X(9,41)=0.41718836125` vs `X(41,9)=0.42077131578`, diff 0.00358 (7x
larger than the IEOMPROP=0 default's 0.00052). This rules out "missing
quadratic correction" as the explanation.

**Scrnc_load_vec.F checked against its trusted ancestor**
(`tdee/tdee_load_vec.F`, found via the routine's own debug label
"@-tdee_load_vec" still present in scrnc's copy) -- diff is 100% clean,
literally identical modulo the renamed subroutine and one checksum-
call-signature difference (an extra unused `S` argument in the
original). List-number usage (`Mu_s=482`/`Mu_d=453`, offsets `1`/`2`
per spin for singles, `IRREPX`-indexed `+1`/`+2`/`+3` for
AAAA/BBBB/ABAB doubles) is self-consistent with how
`Scrnc_form_mubar_s.F`/`_d_uhf.F` write those same lists. No bug found.

**Where this leaves things, honestly**: FIVE scrnc-original or
scrnc-adapted routines across the entire L-side/R-side construction
chain (`Scrnc_form_mubar_00.F`, `Scrnc_form_mubar_s.F`,
`Scrnc_form_mubar_d_uhf.F`, `Scrnc_load_vec.F`, and
`Scrnc_prep_culomb_ints.F`'s raw-integral output) have all been
checked -- either against a decades-trusted vcceh/tdee reference
(clean diffs) or directly and empirically (MIJ symmetry, VO/doubles
decomposition) -- and NONE show an identifiable bug. The remaining
untested pieces are: (1) `Lineqy`/`Plineqy` themselves (general,
widely-used ACES2 infrastructure, so a NEW bug there is less likely a
priori, but scrnc's specific INVOCATION of them -- not the general
routines -- has not been comprehensively re-audited beyond the
Ioffset=Npert_count fix from 2026-08-14); (2)
`Scrnc_form_lrspn_cc_lists.F`'s list-size allocation for 373/374/446/
456 etc.; (3) something about how `Nsizec`/`Scrnc_newsize.F` computes
the JJ-assembly Ddot length, since a size MISMATCH between what 373
and 374 actually contain vs what `Nsizec` assumes could produce
exactly this kind of "close but wrong" cross-term without either side
being individually "wrong" in isolation. (3) is a genuinely new,
not-yet-considered angle and probably the best next lead.

**Scrnc_newsize.F checked too, also clean**: added a debug print
comparing its hardcoded `Idsymsz(Irrepx,1,3)/(2,4)/(13,14)` (AA/BB/AB
doubles sizes) against the same sizes computed via
`Idsymsz(Irrepx,ISYTYP(1,44),ISYTYP(2,44))` etc. (the indirect,
list-number-driven lookup `Scrnc_load_vec.F` uses for the lists that
actually store this data) -- exact match (10/10/100 both ways,
Nsizec=140 total) for H2O/STO-3G. Rules out a block-boundary/size
mismatch between how `Nsizec` is computed for the final Ddot and how
the lists were actually built.

**Six independent checks, all clean, root cause still not found**:
`Scrnc_form_mubar_00.F`, `Scrnc_form_mubar_s.F`,
`Scrnc_form_mubar_d_uhf.F` (all three vs their exact vcceh ancestors),
`Scrnc_load_vec.F` (vs its exact tdee ancestor), raw `MIJ` symmetry
(direct), and `Scrnc_newsize.F`'s sizing (direct, cross-checked two
ways). The IEOMPROP=2 alternative-physics hypothesis was also tested
and refuted (makes the asymmetry 7x worse, not better). At this point,
continuing to check individual routines in isolation has a low hit
rate -- six for six clean. The two remaining candidates worth trying
next time are qualitatively different from what's been tried:
1. `Lineqy`/`Plineqy`'s specific INVOCATION by
   `scrnc_lrspn_cc_main.F` (not the general routines themselves, which
   are trusted, widely-used ACES2 infrastructure) -- has not been
   comprehensively re-audited beyond the 2026-08-14 Ioffset fix.
2. A genuinely from-scratch independent numerical reference (e.g. a
   small Python/NumPy CCSD linear-response calculation for H2O/STO-3G)
   to compare scrnc's raw 373/374 vectors against element-by-element --
   qualitatively different from "read scrnc's own code more carefully,"
   which has now been tried six times without success. Given six clean
   code-reading passes, the bug may be a genuine physics/formula gap
   only visible against an independent calculation, not a local
   indexing slip visible from reading scrnc's code alone.

Once fixed, recheck the sym=on C2v case to see whether it's then
clean, or whether a third, genuinely symmetry-specific issue remains
on top of this one.

**Left in the tree, uncommitted, all `_DEBUG_LVL0`-guarded and
harmless**: the `Npert_count.EQ.9` MIJ-symmetry check in
`scrnc_lrspn_cc_main.F`, the `(9,41)`/`(41,9)` VO/doubles
cross-check in `scrnc_jj_cc_lrspn.F`, and the `Idsymsz` vs `ISYTYP`
comparison in `scrnc_newsize.F` -- all confirmed-clean diagnostics,
kept as templates/regression checks for whoever picks this up next.

---

## SESSION 2026-08-16 (evening, continued) -- independent numerical pilot CONFIRMS the bug is real (not missed physics)

Ajith's direction: build an independent numerical reference, a small
pilot code, no symmetry. Note on framing: PySCF is used here only
because it is a genuinely SEPARATE codebase from ACES II/scrnc, so it
cannot share whatever bug is specific to scrnc's own Fortran -- this is
NOT a claim that PySCF is more authoritative than ACES II. ACES II's
own SCF/CCSD/Lambda are the established, decades-validated reference
this whole project already trusts (and PySCF's own CCSD/Lambda code
was historically itself checked against codes in the ACES/CFOUR
lineage, not the other way around -- there was no PySCF in 1990).
Confirmed first: PySCF UHF energy (-74.96290541582) and CCSD
correlation energy (-0.04933703248) both match ACES II's own
`h2o_scrnc_ccsd_uhf_sto3g` output (-74.9629053675,
-0.049337003434) to 6-7 significant figures -- same geometry/basis,
good starting agreement.

**Deliberately avoided re-deriving scrnc's own response equations by
hand** (Hbar linear solve + Lambda-dressed property-gradient vector) --
doing so would risk introducing an independent bug of its own, which
would defeat the entire point of an independent check. Used a
lower-risk, physically equivalent route instead: the CCSD linear
response reciprocity theorem `<<A;B>>(0) == <<B;A>>(0)` (a direct
consequence of the CCSD energy Lagrangian being a smooth function of
perturbation strengths -- Clairaut/Schwarz symmetry of mixed partial
derivatives, standard 2n+1-rule Lagrangian derivation, Koch & Jorgensen
JCP 93, 3333 (1990)) can be tested via FINITE DIFFERENCE of the
CCSD+Lambda relaxed one-particle density matrix, using ONLY PySCF's own
trusted `cc.UCCSD`/`solve_lambda`/`make_rdm1` end-to-end -- no hand-coded
response formulas anywhere. Frozen HF orbitals throughout (core
Hamiltonian shifted by `eps*B`, orbitals NOT reoptimized), matching
scrnc's own "pure T/Lambda amplitude response, no orbital relaxation"
formulation (Lineqy solves only for T; there is no CPHF/Z-vector piece
anywhere in scrnc).

**Operators used**: diagonal AO-pair "self" operators for the SAME
physical AOs behind scrnc's own worst-offender pair (ACES2 AO-pair
index 9 = AO(2,2) = O 2s self-pair, index 41 = AO(6,6) = one H's 1s
self-pair) -- confirmed via `mol.ao_labels()` that PySCF's own AO
ordering for this exact STO-3G/H2O system places these same physical
orbitals at 0-based indices 1 and 5, i.e. an exact match, no
cross-code AO-ordering translation needed.

**Result** (`aces2/scrnc/independent_pilot/reciprocity_pilot.py`,
committed): at `eps=1e-3`, `d<A>/dlambda_B=-0.39660964`,
`d<B>/dlambda_A=-0.39660751`, relative difference `5.4e-6`. Halving
`eps` to `5e-4` shrinks the relative difference to `2.0e-6` --
consistent with `O(eps^2)` finite-difference truncation error, not a
real effect (a genuine physical asymmetry would NOT shrink when
tightening the finite-difference step). scrnc's own
`JJ(9,41)`/`JJ(41,9)` relative discrepancy is `~1.4e-3` -- roughly
1000x larger than anything finite-difference truncation could produce
here.

**Conclusion, decisively**: the reciprocity theorem genuinely holds
for this exact system and these exact physical operators, confirmed
independently via a completely different computational route
(amplitude-response finite difference vs analytic response-vector dot
product). scrnc's asymmetry is a REAL IMPLEMENTATION BUG, not a missed
physics effect or an inherent limitation of the CI-like/
EOM-effective-Hamiltonian construction. This also independently
corroborates the earlier finding that IEOMPROP=2 ("the quadratic
correction is the documented fix for a gap like this") was the wrong
explanation -- there is no genuine physics gap here to correct.

**Not yet done**: this pilot confirms the SYMMETRY expectation, not
scrnc's exact numerical value -- the raw magnitudes differ in sign and
~10% in scale from scrnc's own numbers, almost certainly just a
perturbation-strength/`Fact=2`/sign-convention difference between "
`d<A>/dlambda_B` for `H=H0+lambda*B`" and whatever exact convention
scrnc's `Ddot(373,374)*Fact` construction uses -- not yet reconciled,
and not needed for the symmetry conclusion above. Reconciling the
convention would let this pilot's `X^A`/`X^B`-equivalent quantities be
compared ELEMENT-BY-ELEMENT against scrnc's raw list-373/374 dumps
(already captured in the `Npert_count.EQ.9`/`(9,41)` diagnostics left
in the tree above), which would be a much sharper localization tool
than the reciprocity check alone -- the natural next step if this
picks back up.

**Committed**: `aces2/scrnc/independent_pilot/reciprocity_pilot.py`,
this log entry.

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
