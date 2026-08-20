













































































































































































































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
C truth, see SCRNC_DEBUG_LOG.md for full numbers): Part 1 ALONE is a
C real, if partial, improvement on the 3 "mixed-symmetry" JJ components
C (9.58% RHF/UHF gap on JJ(1,1) -> 7.45%) and only mildly worsens the
C pre-existing JJ(4,4) anomaly (0.71 -> 1.34 vs UHF's 0.030, was already
C 23.7x off at baseline). Part 1+Part 2 COMBINED improves the mixed
C components further (JJ(1,1) gap -> 3.76%) but makes JJ(4,4) blow up
C catastrophically (33.2, vs 1.34 with Part 1 alone) -- Part 2 is
C therefore DISABLED below pending investigation into why (leading
C guess: JJ(4,4)'s AO-pair symmetry sector already sits close to some
C near-degenerate/small-denominator structure at baseline -- see the
C dead-ends section's note that this component is exactly symmetry-
C decoupled from 1-3 in both references -- and Part 2's one-particle
C shift pushes an RPA-type response denominator through a near-zero
C crossing for that sector specifically). CORR_2P's -1.5 value is
C UNCHANGED from its original derivation guess -- not yet tuned; given
C Part 1 alone only partially closes the gap, it is likely still not
C the exactly correct coefficient, just directionally right.
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
      DOUBLE PRECISION ZERO,ONE,ONEM,CORR_2P
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
      PARAMETER (CORR_2P = -1.5D0)
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
C TEMPORARILY DISABLED to isolate which piece (this one-particle Fmi/
C Fae correction, or the two-particle I(mb,ej) piece above) is
C responsible for JJ(4,4)'s blowup (0.71 -> 33.2) after the first full
C test. See SCRNC_DEBUG_LOG.md.
C     CALL PDCC_QUAD2(ICORE,MAXCOR,IUHF,ONE,8,0)
C     CALL PDCC_QUAD3(ICORE,MAXCOR,IUHF,ONE,8,0)
C
      RETURN
      END
