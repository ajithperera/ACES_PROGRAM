#include "flags.h"
      Subroutine Run1props(Iuhf)

      Implicit None
      Integer Iuhf 
      Logical FileExist

#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"
#include "jobarc.com"

      FileExist=.false.
      inquire(file='IIII',exist=FileExist)
 
      If (FileExist) Then 
         Call Vscf(Icore(i0),Icrsiz,Iuhf)
         Call Propints(icore(i0),icrsiz/iintfp,icore(i0),icrsiz)
      Else 
         Call Ints_(Icore(i0),Icrsiz,Iuhf)
         Call V2ja(icore(i0),Icrsiz,Iuhf)
         Call Vscf(Icore(i0),Icrsiz,Iuhf)
         Call Propints(icore(i0),icrsiz/iintfp,icore(i0),icrsiz)
      Endif 

      If (Iflags(h_IFLAGS_calc) .Gt. 0) Then
         Call aces_init_chemsys()
         Call Vtran(icore(i0),Icrsiz,Iuhf)
         Call aces_init_chemsys()
         Call Init_legacies()
         Call Intproc(icore(i0),Icrsiz,Iuhf)
         Call Vcc(icore(i0),Icrsiz,Iuhf)
         Call Vlambda(icore(i0),Icrsiz,Iuhf)
         Call Vdens(icore(i0),Icrsiz,Iuhf)
         Call Cmp_props(Icore(I0),Icrsiz,Iuhf)
      Else
         Call Cmp_props(Icore(I0),Icrsiz,Iuhf)
      Endif 

      Call aces_fin
      Call C_free(Icore)

      Return
      End

