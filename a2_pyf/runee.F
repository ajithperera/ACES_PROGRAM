      Subroutine Runee(Iuhf)

      Implicit None
      Integer Iuhf

#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"
#include "jobarc.com"

      Call Joda_()
      Call Aces_ja_fin()
      Call Ints_()

      Call aces_init(icore,i0,icrsiz, iuhf, .true.)

      Call V2ja(icore(i0),Icrsiz,Iuhf)
      Call Vscf(Icore(i0),Icrsiz,Iuhf)
      Call aces_init_chemsys()
      Call Vtran(icore(i0),Icrsiz,Iuhf)
      Call aces_init_chemsys()
      Call Init_legacies()
      Call Intproc(icore(i0),Icrsiz,Iuhf)
      Call Vcc(icore(i0),Icrsiz,Iuhf)
      Call Vlambda(icore(i0),Icrsiz,Iuhf)
      Call Vee(icore(i0),Icrsiz,Iuhf)

      Call aces_fin
      Call C_free(Icore)
 
      Return
      End

