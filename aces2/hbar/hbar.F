#include "flags.h"

      SUBROUTINE HBAR_(ICORE,ICRSIZ,IUHF)
      implicit none

      integer icore(icrsiz),icrsiz,iuhf,i0

c COMMON BLOCKS
#include "flags.com"
#include "flags2.com"
#include "machsp.com"

c ----------------------------------------------------------------------
      I0 = 1
C
      Call Drive_Hbar(Icore(i0), icrsiz/iintfp, Iuhf)

C
c ----------------------------------------------------------------------
      RETURN
      end


C ----------------------------------------------------------------------
C Thin PROGRAM wrapper restoring the standalone xHBAR entry point
C (MAIN__), lost when this module was converted to a SUBROUTINE for
C pyaces (see project memory: ACES_PROGRAM unification, 2026-08). This
C calls the shared aces_init/aces_fin pair so the standalone binary
C manages its own memory exactly as the classic driver's own separate
C process used to, before the pyaces conversion removed that from
C inside the subroutine itself.
C ----------------------------------------------------------------------
      PROGRAM XHBAR
      IMPLICIT NONE
#include "icore.com"
#include "istart.com"
      INTEGER IUHF
      CALL ACES_INIT(ICORE, I0, ICRSIZ, IUHF, .TRUE.)
      CALL HBAR_(ICORE(I0), ICRSIZ, IUHF)
      CALL ACES_FIN
      STOP
      END
