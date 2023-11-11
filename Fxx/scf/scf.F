      Subroutine Scf_()

      Implicit None
      Integer Iuhf

c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"

c ----------------------------------------------------------------------
       Call Aces_init(Icore,I0,Icrsiz,Iuhf,.True.)
       Call Init_legacies()
       Call Vscf(Icore(I0),Icrsiz,Iuhf)
       Call Aces_fin
c ----------------------------------------------------------------------
      Return
      End

