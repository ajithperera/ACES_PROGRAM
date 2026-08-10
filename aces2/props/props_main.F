#include "flags.h"

      SUBROUTINE PROPS_(IUHF)
      Implicit none
      Integer iuhf

c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "sympop.com"
#include "machsp.com"
       
c ----------------------------------------------------------------------
      Call Aces_init(icore,i0,icrsiz,iuhf,.true.)

      Call Props(icore(i0),Icrsiz,iuhf)

c ----------------------------------------------------------------------

      RETURN
      End


C ----------------------------------------------------------------------
C Thin PROGRAM wrapper restoring the standalone xprops entry point
C (MAIN__). PROPS_ is self-contained (calls Aces_init itself), so no
C aces_init call is needed here, only aces_fin to close out cleanly.
C ----------------------------------------------------------------------
      PROGRAM XPROPS
      IMPLICIT NONE
      INTEGER IUHF
      CALL PROPS_(IUHF)
      CALL ACES_FIN
      STOP
      END
