      SUBROUTINE RCL_DRIVER_HBARSAFE(WORK,MAXCOR,IUHF,SOLVE_4LAMBDA)
C
C 2026-08-20 (see SCRNC_DEBUG_LOG.md): RCL_DRIVER (the real iterative
C Lambda-equation solver, always invoked from Rcl_main whenever
C DENSITY=RESPONSE is set, regardless of the excite= keyword) reads
C AND OVERWRITES lists 54(-59) -- the same lists Rcl_hbar/RCC_GENINT
C built as the plain ground-state Hbar, and the same lists scrnc's own
C Scrnc_form_rcc_hbar.F reads downstream, in the SEPARATE xscrnc
C process, after xringcl has fully exited.
C
C A save-before/restore-after wrapper around RCL_DRIVER (first attempt,
C see git history) does NOT work: RCL_DRIVER itself calls
C `CALL ACES_FIN` + `CALL ACES_EXIT(0)` internally upon completing the
C Lambda solve (same pattern RCC_DRIVER.F uses in xringcc) -- it never
C returns to its caller, so a restore-after-the-call is architecturally
C impossible; confirmed empirically (a checksum placed immediately
C after the CALL RCL_DRIVER line never printed).
C
C Fix instead: snapshot lists 54(-59) to NEW, dedicated list numbers
C (154/156, offset +100, mirroring the DCC precedent's own +200 scratch-
C list convention elsewhere in this codebase) BEFORE calling
C RCL_DRIVER, so the ORIGINAL ground-state Hbar survives under a
C different name regardless of what RCL_DRIVER does to 54(-59)
C afterward or whether it ever returns. Scrnc_form_rcc_hbar.F reads
C from 154/156 instead of 54/56 for RHF RCCD/DRCCD. Called from
C Rcl_main.F in place of a direct `Call Rcl_driver(...)`.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SOLVE_4LAMBDA
      DIMENSION WORK(MAXCOR)
      INTEGER IRREPX,IMODE
C
#include "machsp.com"
#include "sympop.com"
C
      IRREPX = 1
      IMODE  = 0
C
      IF (IUHF .NE. 0) THEN
         CALL INIPCK(1,ISYTYP(1,54),ISYTYP(2,54),154,IMODE,0,1)
         CALL INIPCK(1,ISYTYP(1,55),ISYTYP(2,55),155,IMODE,0,1)
         CALL INIPCK(1,ISYTYP(1,56),ISYTYP(2,56),156,IMODE,0,1)
         CALL INIPCK(1,ISYTYP(1,57),ISYTYP(2,57),157,IMODE,0,1)
         CALL INIPCK(1,ISYTYP(1,58),ISYTYP(2,58),158,IMODE,0,1)
         CALL INIPCK(1,ISYTYP(1,59),ISYTYP(2,59),159,IMODE,0,1)
      ELSE
         CALL INIPCK(1,ISYTYP(1,54),ISYTYP(2,54),154,IMODE,0,1)
         CALL INIPCK(1,ISYTYP(1,56),ISYTYP(2,56),156,IMODE,0,1)
         CALL INIPCK(1,ISYTYP(1,58),ISYTYP(2,58),158,IMODE,0,1)
      ENDIF
C
      IF (IUHF .NE. 0) THEN
         ISZ54 = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         ISZ55 = IDSYMSZ(IRREPX,ISYTYP(1,55),ISYTYP(2,55))
         ISZ56 = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         ISZ57 = IDSYMSZ(IRREPX,ISYTYP(1,57),ISYTYP(2,57))
         ISZ58 = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         ISZ59 = IDSYMSZ(IRREPX,ISYTYP(1,59),ISYTYP(2,59))
      ELSE
         ISZ54 = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         ISZ56 = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         ISZ58 = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
      ENDIF
C
      I900 = 1
      I910 = I900 + ISZ54*IINTFP
      IF (IUHF .NE. 0) THEN
         I920 = I910 + ISZ55*IINTFP
         I930 = I920 + ISZ56*IINTFP
         I940 = I930 + ISZ57*IINTFP
         I950 = I940 + ISZ58*IINTFP
         IEND = I950 + ISZ59*IINTFP
      ELSE
         I920 = I910 + ISZ56*IINTFP
         IEND = I920 + ISZ58*IINTFP
      ENDIF
      IF (IEND .GT. MAXCOR) CALL INSMEM('RCL_DRIVER_HBARSAFE',IEND,
     &                                   MAXCOR)
C
      CALL GETALL(WORK(I900),ISZ54,IRREPX,54)
      CALL PUTALL(WORK(I900),ISZ54,IRREPX,154)
      IF (IUHF .NE. 0) THEN
         CALL GETALL(WORK(I910),ISZ55,IRREPX,55)
         CALL PUTALL(WORK(I910),ISZ55,IRREPX,155)
         CALL GETALL(WORK(I920),ISZ56,IRREPX,56)
         CALL PUTALL(WORK(I920),ISZ56,IRREPX,156)
         CALL GETALL(WORK(I930),ISZ57,IRREPX,57)
         CALL PUTALL(WORK(I930),ISZ57,IRREPX,157)
         CALL GETALL(WORK(I940),ISZ58,IRREPX,58)
         CALL PUTALL(WORK(I940),ISZ58,IRREPX,158)
         CALL GETALL(WORK(I950),ISZ59,IRREPX,59)
         CALL PUTALL(WORK(I950),ISZ59,IRREPX,159)
      ELSE
         CALL GETALL(WORK(I910),ISZ56,IRREPX,56)
         CALL PUTALL(WORK(I910),ISZ56,IRREPX,156)
         CALL GETALL(WORK(I920),ISZ58,IRREPX,58)
         CALL PUTALL(WORK(I920),ISZ58,IRREPX,158)
      ENDIF
C
      CALL RCL_DRIVER(WORK,MAXCOR,IUHF,SOLVE_4LAMBDA)
C
C RCL_DRIVER calls ACES_FIN+ACES_EXIT(0) internally on completion and
C never returns -- any code here is unreachable in practice, kept only
C in case some future RCL_DRIVER code path DOES return normally.
      RETURN
      END
