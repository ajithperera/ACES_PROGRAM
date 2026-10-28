      SUBROUTINE SCRNC_RCC_HBAR_READ(W,MAXCOR,IUHF,ILOC,IMULT)
C
C Reads the raw drCCD Hbar matrix (RHF, or UHF AA/BB/AABB/BBAA blocks)
C into W, indexed per irrep via ILOC. Split out of the old
C Scrnc_form_rcc_hbar_inv.F (2026-08-15) so it can run exactly ONCE,
C before Scrnc_driver.F's dynamic +/-Freq Icount loop, instead of being
C re-read via GETLST on every pass. Repeated GETLST calls to these
C lists were found empirically to add a spurious, constant, orbital-
C energy-gap-sized increment to the diagonal on every call (root cause
C not fully traced to its source -- confirmed additive and exactly
C reproducible across repeated calls even with Freq=0 throughout, see
C SCRNC_DEBUG_LOG.md). Reading once into a pristine in-memory copy that
C Scrnc_form_rcc_hbar_inv.F then only COPIES from (never re-reads)
C avoids the issue entirely, whatever its root cause.
C
C 2026-08-19: added IMULT (unused before -- this routine always read
C list 56, even during Scrnc_driver.F's Imult=2/triplet pass of its
C RHF RCCD two-channel loop). Scrnc_form_rcc_hbar.F already writes the
C singlet Hbar to list 56 (Imult=1) and the triplet Hbar to list 54
C (Imult=2) -- this reader never picked up that distinction, so the
C "triplet channel" response was actually being built from a second,
C redundant read of the SINGLET Hbar, mislabeled as triplet. Confirmed
C via a from-scratch Hbar-eigenvalue check (RHF's two separate
C singlet/triplet diagonalizations match UHF's combined one to 7
C decimal places -- see SCRNC_DEBUG_LOG.md), which proves list 54
C itself is correct; this routine just never read it.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSIZ
      LOGICAL CIS
      DOUBLE PRECISION W(MAXCOR)
      DIMENSION ILOC(8)

      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREP0(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
#include "flags.com"

      CIS = .FALSE.

      IF (CIS) THEN
         LISTDIAG = 42
      ELSE IF (IUHF.EQ.0 .AND. IMULT.EQ.2) THEN
         LISTDIAG = 54
      ELSE
         LISTDIAG = 56
      ENDIF

      DO IRREP=1,NIRREP

         IF (IUHF.EQ.0)THEN
            NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTDIAG))
            I000=ILOC(IRREP)
            CALL GETLST(W(I000),1,NUMDIS,1,IRREP,LISTDIAG)
         ELSE

            NUMAA=IRPDPD(IRREP,ISYTYP(1,19))
            NUMBB=IRPDPD(IRREP,ISYTYP(1,20))
            MATDIM=NUMAA+NUMBB

            IOFFAAAA=ILOC(IRREP)
            IOFFBBBB=IOFFAAAA+MATDIM*NUMAA+NUMAA
            IOFFBBAA=IOFFAAAA+NUMAA
            IOFFAABB=IOFFAAAA+MATDIM*NUMAA

            IF (CIS) THEN
               LISTW=23
            ELSE
               LISTW=54
            ENDIF

            NUMDIS=IRPDPD(IRREP,ISYTYP(1,LISTW))
            IOFF=IOFFAAAA
            DO 11 IDIS=1,NUMDIS
               CALL GETLST(W(IOFF),IDIS,1,1,IRREP,LISTW)
               IOFF=IOFF+MATDIM
11          CONTINUE

            IF (CIS) THEN
               LISTW=24
            ELSE
               LISTW=55
            ENDIF

            NUMDIS=IRPDPD(IRREP,ISYTYP(1,LISTW))
            IOFF=IOFFBBBB
            DO 12 IDIS=1,NUMDIS
               CALL GETLST(W(IOFF),IDIS,1,1,IRREP,LISTW)
               IOFF=IOFF+MATDIM
12          CONTINUE

            IF (CIS) THEN
               LISTW1=17
               LISTW2=17
            ELSE
               LISTW1=56
               LISTW2=57
            ENDIF

            NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
            IOFF1=IOFFBBAA
            IOFF2=IOFFAABB

            DO 13 IDIS=1,NUMDIS
               CALL GETLST(W(IOFF1),IDIS,1,1,IRREP,LISTW1)
               CALL GETLST(W(IOFF2),IDIS,1,1,IRREP,LISTW2)
               IOFF1=IOFF1+MATDIM
               IOFF2=IOFF2+MATDIM
13          CONTINUE

         ENDIF

      ENDDO
      RETURN
      END
