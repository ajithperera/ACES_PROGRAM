      Subroutine Tran_()

      Implicit None
      Integer Iuhf

c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"
#include "jobarc.com"

c ----------------------------------------------------------------------
       Call aces_init(icore,i0,icrsiz, iuhf, .true.)
       Call Init_legacies()
       Call Vtran(icore(i0),Icrsiz,Iuhf)
       Call C_free(Icore)
       Call aces_fin
c ----------------------------------------------------------------------
      Return
      End

