













































































































































































































      SUBROUTINE RCC_LAGRANGIAN_FIX(ICORE,MAXCOR,IUHF,RCCD,
     &                                    DRCCD)
C
C Ports the DCCSD Lagrangian-derived Hbar correction pattern (see
C aces2/vcc/post_dcc_mods.F, dcc_fixfbar.F, pdcc_quad2/3.F) to scrnc's
C RHF closed-shell RCCD/DRCCD pipeline. DCCSD has the same structural
C problem as (d)rCCD -- no compact e^T exponential wavefunction, so the
C naive Hbar=exp(-T)Hexp(T) formula does not apply and the correct
C expressions come from the Lagrangian equations instead. This is a
C first, empirically-tuned port, not an independently re-derived rCCD
C Lagrangian -- see SCRNC_DEBUG_LOG.md for the derivation this is based
C on and what remains uncertain.
C
C Two pieces, applied once, called from RCC_DRIVER right after the
C amplitude iterations converge (same process/JOBARC-list-namespace
C RCC_GENINT/RCC_DWMBEJ themselves ran in -- an earlier attempt calling
C this from Scrnc_driver instead, in the separate xscrnc process/
C module, got garbage dimensions for list 18/37 because scrnc's own
C pipeline steps reuse those same list numbers for unrelated purposes.
C See SCRNC_DEBUG_LOG.md.
C
C  1) Two-particle: subtract a self-referential
C     I(mb,ej) = t(mn,ef)<nj||fb> contraction from the already-built
C     W(mb,ej) sitting on lists 56 (singlet) / 54 (triplet) -- built
C     via the SAME CNTRCT machinery RCC_DWMBEJ itself uses (copied,
C     not shared, to avoid touching RCC_DWMBEJ's own iteration-time
C     use), with a unit (+1.0) coefficient instead of RCC_DWMBEJ's
C     own -C4-derived (+0.5) one, and SUBTRACTED (not added).
C
C  2) One-particle: fill scrnc_form_rcc_hbar.F's currently-all-zero
C     (see its own NONHF/HEFF=.FALSE. branch) Fmi/Fae quadratic blocks
C     via pdcc_quad2.F/pdcc_quad3.F -- confirmed GENERIC vcc utility
C     routines (not DCC-specific), reading RCC's own converged T2 off
C     list 46 (confirmed populated -- see RCC_DRIVER.F) and the
C     standard MO integral list 16. Stored at list 91/92 sub-index 9
C     (IOFF=8), the same slot DCC itself uses, to stay clear of
C     whatever (if anything) uses sub-index 1/2. List 91/92 sub-index 9
C     is not registered by any pure ring-CC pipeline step by default --
C     UPDMOI'd here first (see below), mirroring post_dcc_mods.F.
C
C STATUS (2026-08-20, H2O/STO-3G RHF vs UHF "RPA JJ matrix" ground
C truth, see SCRNC_DEBUG_LOG.md for full numbers, findings below and
C the CORR_2P scan/JJ(4,4)-blowup investigation are from a later
C session than the original -1.5 derivation): Part 1 ALONE is a real,
C if partial and INCOMPLETE, improvement on JJ(1,1) specifically
C (9.58% RHF/UHF gap -> 7.45% at CORR_2P=-1.5) and only mildly worsens
C the pre-existing JJ(4,4) anomaly. Part 1+Part 2 COMBINED improves
C JJ(1,1) further but makes JJ(4,4) blow up catastrophically --
C DIAGNOSED (not just observed): this is a genuine near-pole crossing,
C not a Part-2 coding bug -- a FACT2 strength scan (0/0.5/1.0) showed
C the ring-CC Hbar's 2 smallest eigenvalues (the SAME pre-existing
C "2 negative eigenvalues" bug documented earlier in this file, see the
C "SESSION 2026-08-20 (continued)" entries in SCRNC_DEBUG_LOG.md) cross
C zero LINEARLY and monotonically as Part 2 strength increases
C (-0.633/-0.540 -> -0.255/-0.142 -> +0.122/+0.256), landing very close
C to (not comfortably past) zero at full strength -- exactly the
C signature of a static response function passing near an excitation
C pole, not an implementation defect. Both corrections move these
C eigenvalues the SAME direction UHF's own all-positive spectrum
C implies is correct.
C
C A separate CORR_2P scan (Part 2 OFF, to stay away from that pole) at
C -1,-2,-3,-6,-10,-12 found JJ(1,1)'s gap shrinks monotonically toward
C 0 as CORR_2P becomes more negative (crossing zero near ~-10.5) --
C BUT JJ(2,2) and JJ(3,3)'s gaps get monotonically WORSE over the exact
C same range (baseline, CORR_2P=0, is their own best value). No single
C scalar CORR_2P improves all 3 "well-behaved" JJ components at once --
C this means the correction's FUNCTIONAL FORM (a single uniformly-
C scaled self-referential T*W contraction), not just its magnitude, is
C structurally incomplete; further coefficient scanning alone will not
C fix this. CORR_2P is left at -1.5 (its original derivation guess) --
C not because it is confirmed correct, but because no single scalar
C value is a strict improvement across all components, so there is no
C better single number to prefer. NEXT: needs either an actual
C Lagrangian re-derivation (not a DCC-analogy guess) or a genuinely
C different functional form (e.g. matrix-element/symmetry-block-
C dependent weighting) rather than more scalar tuning.
C
CEND
C 2026-08-20: this whole class of ringcc routine (rcc_dwmbej.F,
C cntrct.F) uses IMPLICIT INTEGER (A-Z) -- dimension/list-index
C variables (DSSYWA, NMSYWA, IRREP, NSIZE, ...) outnumber the float
C locals here, and matching that convention (rather than the more
C common ACES2 IMPLICIT DOUBLE PRECISION) avoids the classic Fortran
C implicit-typing trap. A first version of this file used IMPLICIT
C DOUBLE PRECISION and only explicitly declared 3 of its own index
C variables -- every OTHER index var (DSSYWA, DSSYTA, NMSYWA, NMSYTA,
C NUMSYQ, IRREP, ...) silently fell through to DOUBLE PRECISION,
C corrupting every list-dimension lookup and the CNTRCT/DGEMM calls
C downstream (confirmed via debug prints: ISYTYP itself was fine,
C 9/10 for lists 18/37/56 alike; IRPDPD(irrep,9) computed correctly
C into an INTEGER local (DISSYQ) but came out as reinterpreted-bit
C garbage into a same-valued but DOUBLE-PRECISION-typed local,
C DSSYWA). See SCRNC_DEBUG_LOG.md and [[feedback_aces_fortran_implicit_typing]].
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ZERO,ONE,ONEM,CORR_2P,FACT2
      LOGICAL RCCD,DRCCD
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT1AA,NT1BB,NF1AA,NF1BB,NF2AA,
     &             NF2BB
      COMMON /FLAGS/ IFLAGS(100)
C
      DATA ZERO,ONE,ONEM /0.0D0,1.0D0,-1.0D0/
C
C Empirically-tuned overall coefficient for the two-particle piece --
C see header. Applied as: Hbar_new = Hbar_old + CORR_2P * I_unit, where
C I_unit is the SAME T*W contraction RCC_DWMBEJ builds with a unit
C (+1.0) scale. -1.5 makes the net two-particle term go from RCC's
C current "+0.5*contraction" to DCC's analogous "-1.0*contraction"
C (0.5 + (-1.5) = -1.0) -- see SCRNC_DEBUG_LOG.md for the derivation.
C NOT yet empirically confirmed against the JJ-matrix ground truth.
C
C 2026-08-20: the ORIGINAL "SUPERSEDED" note here (see git history)
C referred to the BROKEN rebuild-based pipeline fix (e48b5cc0), since
C corrected (c4169bec, snapshot-based, validated bit-identical against
C ringcl's independent Hbar diagonalization). This routine still runs
C INSIDE xringcc, BEFORE the snapshot is taken in xringcl -- so unlike
C the broken rebuild-based fix, THIS correction's edits to list 56 DO
C survive into the snapshot scrnc reads. Re-tested against the
C validated baseline at CORR_2P=-1.5 and -10: NEITHER helps. JJ11 is
C nearly insensitive to CORR_2P now (-6.16% -> -6.22% from 0 to -10),
C while JJ22/JJ33/JJ44 all get WORSE as CORR_2P goes more negative
C (JJ22: -12.91% -> -14.71%; JJ33: -6.49% -> -7.21%; JJ44: -13.20% ->
C -13.63%) -- a DIFFERENT pattern than the old pre-pipeline-fix scan
C (where JJ11 was the dominant, sensitive mover). CORR_2P=0 (no
C correction at all) is the best available value against the validated
C baseline -- left at 0. The residual ~6-13% undershoot needs a
C genuinely different fix (real Lagrangian re-derivation, not scalar
C tuning on this formula) -- see SCRNC_DEBUG_LOG.md.
      PARAMETER (CORR_2P = 0.0D0)
C
      IF (IUHF .NE. 0) RETURN
      IF (.NOT. (RCCD .OR. DRCCD)) RETURN
C
C ============================================================
C  Part 1: two-particle I(mb,ej) correction, singlet (list 56)
C ============================================================
C
      LISTTA = 37
      LISTWA = 18
      LISTWB = 21
      LISTQ  = 56
C
      DO IRREP = 1, NIRREP
        DSSYWA = IRPDPD(IRREP,ISYTYP(1,LISTWA))
        DSSYTA = IRPDPD(IRREP,ISYTYP(1,LISTTA))
        NMSYWA = IRPDPD(IRREP,ISYTYP(2,LISTWA))
        NMSYTA = IRPDPD(IRREP,ISYTYP(2,LISTTA))
        DISSYQ = IRPDPD(IRREP,ISYTYP(1,LISTQ))
        NUMSYQ = IRPDPD(IRREP,ISYTYP(2,LISTQ))
C
        I000 = 1
        I010 = I000 + IINTFP*DISSYQ*NUMSYQ
        I020 = I010 + IINTFP*NMSYWA*DSSYWA
        I030 = I020 + IINTFP*NMSYTA*DSSYTA
        I040 = I030 + IINTFP*DISSYQ*NUMSYQ
        IF (I040 .GT. MAXCOR) CALL INSMEM('SCRNC_RCC_LFIX',I040,
     &                                     MAXCOR)
C
        CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
        CALL GETLST(ICORE(I020),1,NMSYTA,1,IRREP,LISTTA)
        CALL CNTRCT(ICORE(I000),ICORE(I010),ICORE(I020),
     &              DISSYQ,NUMSYQ,DSSYTA,NMSYTA,DSSYWA,NMSYWA,
     &              IRREP,LISTWA,LISTWB,LISTTA,'WxT',ZERO,ONE,1)
C
        NSIZE = DISSYQ*NUMSYQ
        CALL GETALL(ICORE(I030),NSIZE,IRREP,LISTQ)
        CALL SAXPY(NSIZE,CORR_2P,ICORE(I000),1,ICORE(I030),1)
        CALL PUTALL(ICORE(I030),NSIZE,IRREP,LISTQ)
      ENDDO
C
C ============================================================
C  Part 1b: two-particle I(mb,ej) correction, triplet (list 54)
C ============================================================
C
      LISTTB = 34
      LISTQ  = 54
C
      DO IRREP = 1, NIRREP
        DSSYWA = IRPDPD(IRREP,ISYTYP(1,LISTWA))
        DSSYTB = IRPDPD(IRREP,ISYTYP(1,LISTTB))
        NMSYWA = IRPDPD(IRREP,ISYTYP(2,LISTWA))
        NMSYTB = IRPDPD(IRREP,ISYTYP(2,LISTTB))
        DISSYQ = IRPDPD(IRREP,ISYTYP(1,LISTQ))
        NUMSYQ = IRPDPD(IRREP,ISYTYP(2,LISTQ))
C
        I000 = 1
        I010 = I000 + IINTFP*DISSYQ*NUMSYQ
        I020 = I010 + IINTFP*NMSYWA*DSSYWA
        I030 = I020 + IINTFP*NMSYTB*DSSYTB
        I040 = I030 + IINTFP*DISSYQ*NUMSYQ
        IF (I040 .GT. MAXCOR) CALL INSMEM('SCRNC_RCC_LFIX',I040,
     &                                     MAXCOR)
C
        CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
        CALL GETLST(ICORE(I020),1,NMSYTB,1,IRREP,LISTTB)
        CALL RCC_TRP_CNTRCT(ICORE(I000),ICORE(I010),ICORE(I020),
     &              DISSYQ,NUMSYQ,DSSYTB,NMSYTB,DSSYWA,NMSYWA,
     &              IRREP,LISTWA,LISTWB,LISTTB,'WxT',ZERO,ONE,1)
C
        NSIZE = DISSYQ*NUMSYQ
        CALL GETALL(ICORE(I030),NSIZE,IRREP,LISTQ)
        CALL SAXPY(NSIZE,CORR_2P,ICORE(I000),1,ICORE(I030),1)
        CALL PUTALL(ICORE(I030),NSIZE,IRREP,LISTQ)
      ENDDO
C
C ============================================================
C  Part 2: one-particle Fmi/Fae quadratic correction
C ============================================================
C
C List 91/92 sub-index 9 (the IOFF=8 slot PDCC_QUAD2/3 target) is never
C registered by any pure ring-CC pipeline step (unlike a plain CCSD run,
C where xvcc's own setup does this) -- PUTLST asserts "list does not
C exist" without it. Register via UPDMOI, mirroring post_dcc_mods.F's
C own `Call Updmoi(1,IHHA_LENGTH,9,91,0,0)` pattern exactly.
      IHHA_LENGTH = IRPDPD(1,21)
      IPPA_LENGTH = IRPDPD(1,19)
      CALL UPDMOI(1,IHHA_LENGTH,9,91,0,0)
      CALL UPDMOI(1,IPPA_LENGTH,9,92,0,0)
C
C 2026-08-20: original FACT2 scan (0, 0.5, 1.0) against the OLD,
C corrupted (pre-pipeline-fix) baseline found the 2 small Hbar
C eigenvalues crossing zero LINEARLY as Part-2 strength increases, with
C JJ(4,4) blowing up catastrophically near the crossing (1.34 -> 33.2).
C
C 2026-08-21: re-tested FACT2=1.0 against the NOW-VALIDATED (c4169bec)
C baseline -- and discovered the FIRST attempt at this re-test (same
C day) was itself a false negative: scrnc_form_rcc_hbar.F's read of
C list 91/92 sub-index 9 (Part 2's actual output) was commented out
C (dead wire, see that file's own history), so EVERY prior FACT2 test,
C including the ORIGINAL pre-pipeline-fix blowup scan above, may never
C have actually exercised Part 2's effect on Hbar(i,j)/Hbar(a,b) at
C all -- unresolved which of PDCC_QUAD2/3's other side effects (if any)
C produced that original blowup instead; flagged, not chased further.
C Re-enabled the read (this session) and re-ran FACT2=1.0 for real: NO
C blowup, no near-zero eigenvalue crossing this time -- instead a
C clean, coherent WORSENING of all 4 JJ components (JJ11 -6.16%->
C -6.62%, JJ22 -12.91%->-13.32%, JJ33 -6.49%->-6.75%, JJ44 -13.20%->
C -14.24%). Unambiguous negative result, same conclusion as CORR_2P's
C own scan: Part 2 in its current form does not help. FACT2 left at 0.
C See SCRNC_DEBUG_LOG.md.
      FACT2 = 0.0D0
      IF (FACT2 .NE. 0.0D0) THEN
        CALL PDCC_QUAD2(ICORE,MAXCOR,IUHF,FACT2,8,0)
        CALL PDCC_QUAD3(ICORE,MAXCOR,IUHF,FACT2,8,0)
      ENDIF
C
      RETURN
      END
