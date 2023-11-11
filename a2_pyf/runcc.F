      Subroutine Runcc(Iuhf)

      Implicit None
      Integer Iuhf 
      Integer Length 
      Double Precision Tmp
      Logical Scf_done

#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"
#include "jobarc.com"


      Call Getrec(0,"JOBARC", "SCFEVCA0",Length,Tmp)
      Scf_done = .False. 
      If (Length .Gt. 0) Scf_done = .True.

      If (Scf_done) Then
         Call aces_init_chemsys()
         Call Vtran(icore(i0),Icrsiz,Iuhf)
         Call aces_init_chemsys()
         Call Init_legacies()
         Call Intproc(icore(i0),Icrsiz,Iuhf)
         Call Vcc(icore(i0),Icrsiz,Iuhf)
      Else
         Call Ints_(Icore(i0),Icrsiz,Iuhf)
         Call V2ja(icore(i0),Icrsiz,Iuhf)
         Call Vscf(Icore(i0),Icrsiz,Iuhf)
         Call aces_init_chemsys()
         Call Vtran(icore(i0),Icrsiz,Iuhf)
         Call aces_init_chemsys()
         Call Init_legacies()
         Call Intproc(icore(i0),Icrsiz,Iuhf)
         Call Vcc(icore(i0),Icrsiz,Iuhf)
      Endif 

      Call aces_fin
      Call C_free(Icore)

      Return
      End

