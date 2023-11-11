
      Subroutine Props_()
      Implicit none

#include "icore.com"
#include "istart.com"
#include "machsp.com"

      Integer iuhf

      Call Aces_init(icore,i0,icrsiz, iuhf, .true.)
      Call Propints(icore(i0),icrsiz/iintfp,icore(i0),icrsiz)
      Call Aces_fin
      Call C_free(Icore)

      Return
      End

