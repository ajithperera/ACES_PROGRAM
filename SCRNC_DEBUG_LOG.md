# scrnc / GW manuscript debug log

Plain-text resume document, kept in the repo so it survives terminal
freezes / HPG outages and gets backed up on every `git push`. Read the
**RESUME POINT** section first on any restart.

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
