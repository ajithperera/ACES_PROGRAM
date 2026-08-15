













































































































































































































      SUBROUTINE SCRNC_FORM_RCC_HBAR_INV(WHBAR,W,MAXCOR,IUHF,ILOC)
C
C THIS ROUTINE DRIVES THE drCCD Hbar diagonalizations
C
C 2026-08-15: reworked to take the pristine raw Hbar (WHBAR, indexed
C per irrep via ILOC, read ONCE by the new Scrnc_rcc_hbar_read.F) as
C input, instead of re-reading it via GETLST itself. This routine only
C COPIES from WHBAR into its own output buffer W (also ILOC-indexed,
C a separate memory region) before shifting/inverting -- WHBAR is never
C written to, so it stays a valid pristine snapshot across repeated
C calls (one per +/-Freq pass of Scrnc_driver.F's dynamic Icount loop).
C Necessary because repeated GETLST re-reads of the underlying lists
C were found to add a spurious additive increment on every call (see
C Scrnc_rcc_hbar_read.F's header and SCRNC_DEBUG_LOG.md).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSIZ
      LOGICAL CIS
      LOGICAL REALFREQ
      DOUBLE PRECISION WHBAR(MAXCOR),W(MAXCOR)
      DIMENSION ILOC(8)

      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREP0(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/FREQINFO/FREQ,REALFREQ,XFREQ(200),NO_FREQS
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end

      DATA TWO, TOL /2.0D0,1.0D-08/
      Write(6,"(1x,a)") "-----Entering scrnc_form_rcc_hbar_inv------"
C
      CIS = .FALSE.

      DO IRREP=1,NIRREP

         IF (IUHF.EQ.0)THEN
            NUMDIS=IRPDPD(IRREP,ISYTYP(2,56))
            MATDIM = NUMDIS
         ELSE
            NUMAA=IRPDPD(IRREP,ISYTYP(1,19))
            NUMBB=IRPDPD(IRREP,ISYTYP(1,20))
            MATDIM=NUMAA+NUMBB
         ENDIF

         I000 = ILOC(IRREP)
         I010 = I000 + MATDIM * MATDIM
         I020 = I010 + MATDIM
         IEND = I020
         IF (IEND .GT. MAXCOR) CALL INSMEM("SCRNC_FORM_RCC_HBAR_INV",
     +                                     IEND,MAXCOR)

         CALL DCOPY(MATDIM*MATDIM,WHBAR(I000),1,W(I000),1)

        Write(6,"(a)") " The drCCD Hbar matrix"
        Iprint=Iloc(Irrep)
        call output(W(Iprint),1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
C Do a full non-symmetric diagonalization

C Independent Hbar diagonalization (Ajith's request 2026-08-14), to
C compare against ringcc's own eigenvalues from Rcl_check_rcchbar (see
C SCRNC_DEBUG_LOG.md). DGEEV destroys its input matrix in-place, so
C snapshot Hbar into scratch first -- this is a read-only diagnostic and
C must not perturb the production Hbar^(-1) that MINV builds below.
         Ihcopy = IEND
         Ievr   = Ihcopy + MATDIM*MATDIM
         Ievi   = Ievr   + MATDIM
         Ivl    = Ievi   + MATDIM
         Ivr    = Ivl    + MATDIM*MATDIM
         Idgwrk = Ivr    + MATDIM*MATDIM
         Idgend = Idgwrk + 4*MATDIM
         IF (Idgend .GT. MAXCOR) CALL INSMEM(
     +                        "SCRNC_FORM_RCC_HBAR_INV-diag",
     +                        Idgend,MAXCOR)

         CALL DCOPY(MATDIM*MATDIM,W(I000),1,W(Ihcopy),1)
         CALL DGEEV("V","V",MATDIM,W(Ihcopy),MATDIM,W(Ievr),W(Ievi),
     +              W(Ivl),MATDIM,W(Ivr),MATDIM,W(Idgwrk),4*MATDIM,
     +              Ierror)

         Write(6,*)
         Write(6,"(a)") " scrnc's own Hbar(ai,bj) is fully diagonalized"
         Write(6,*)
         Bohrtoev = 27.2113961D0
         Write(6,"(a)") " The real eigenvalues"
         Write(6,"(5(1x,F15.8))") (W(Ievr+I-1)*Bohrtoev,I=1,MATDIM)
         Write(6,"(a)") " The imaginary eigenvalues"
         Write(6,"(5(1x,F15.8))") (W(Ievi+I-1)*Bohrtoev,I=1,MATDIM)

C Independent cross-check of the (Hbar^-1)^T shortcut used by the L-side
C in scrnc_jj_rcc_lrspn.F: (A^-1)^T = (A^T)^-1 is a pure linear-algebra
C identity, so it must hold to numerical precision for ANY invertible A,
C symmetric or not. Verify it directly here rather than trusting the
C algebra -- build Hbar^T explicitly and invert IT independently, then
C compare against the transpose of the (already-computed-elsewhere)
C Hbar^-1. Copy from W(I000), which still holds pristine Hbar here (the
C DGEEV diagnostic above destroyed its OWN copy at Ihcopy, not this one;
C the main MINV call below is what destroys W(I000), and that hasn't
C run yet).
         Ihbt  = Idgend
         Iscr2 = Ihbt + MATDIM*MATDIM
         Idend = Iscr2 + MATDIM
         IF (Idend .GT. MAXCOR) CALL INSMEM(
     +                        "SCRNC_FORM_RCC_HBAR_INV-xpchk",
     +                        Idend,MAXCOR)
         DO Jc = 1, MATDIM
            DO Ic = 1, MATDIM
               W(Ihbt-1+(Jc-1)*MATDIM+Ic) =
     +            W(I000-1+(Ic-1)*MATDIM+Jc)
            ENDDO
         ENDDO
         CALL MINV(W(Ihbt),MATDIM,MATDIM,W(Iscr2),Dett,TOL,0)

C Dynamic (frequency-dependent) response, 2026-08-15: build
C (Hbar + Freq*1)^(-1) instead of Hbar^(-1) when Freq is nonzero (sign
C confirmed empirically: gives the same increasing-toward-resonance
C trend as the CCSD path's Lineqy/Plineqy-driven response). The caller
C (Scrnc_driver.F) rebuilds this fresh from the pristine WHBAR for each
C +/-Freq pass of its dynamic Icount loop and shares the resulting
C shifted inverse with both the R-side (Scrnc_lrspn_rcc_main.F) and
C L-side (Scrnc_jj_rcc_lrspn.F) routines for that pass, same as the
C static Hbar^(-1) always was.

         IF (FREQ .NE. 0.0D0) THEN
            DO ID1 = 1, MATDIM
               W(I000-1+(ID1-1)*MATDIM+ID1) =
     +            W(I000-1+(ID1-1)*MATDIM+ID1) + FREQ
            ENDDO
         ENDIF

         CALL MINV(W(I000),MATDIM,MATDIM,W(I010),DET,TOL,0)
        Iprint=Iloc(Irrep)
        Write(6,"(a)") " The drCCD Hbar inverse matrix"
        call output(W(Iprint),1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)

        Xpmax = 0.0D0
        DO Jc = 1, MATDIM
           DO Ic = 1, MATDIM
C (Hbarinv)^T(Ic,Jc) = Hbarinv(Jc,Ic); compare to (Hbar^T)^-1(Ic,Jc)
              Xp1 = W(I000-1+(Ic-1)*MATDIM+Jc)
              Xp2 = W(Ihbt-1+(Jc-1)*MATDIM+Ic)
              Xpmax = MAX(Xpmax,ABS(Xp1-Xp2))
           ENDDO
        ENDDO
        Write(6,"(a,D15.6)") " max|(Hbarinv)^T - (Hbar^T)^-1| = ",
     +                        Xpmax

      ENDDO
      RETURN
      END
