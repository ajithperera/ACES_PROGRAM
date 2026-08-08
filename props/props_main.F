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

