      Subroutine Ee_()
CSSS      Program Ee
      Implicit None
      Integer Iuhf

c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"

c ----------------------------------------------------------------------
       Call aces_init(icore,i0,icrsiz, iuhf, .true.)
       Call init_legacies()
       Call Vee(icore(i0),Icrsiz,Iuhf)
       Call aces_fin
       Call C_free(Icore)
c ----------------------------------------------------------------------
      Return
CSSS      Stop
      End

