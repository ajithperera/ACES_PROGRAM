      Subroutine Runscf(Iuhf)

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
      Else 
         Call Ints_(Icore(i0),Icrsiz,Iuhf)
         Call V2ja(icore(i0),Icrsiz,Iuhf)
         Call Vscf(Icore(i0),Icrsiz,Iuhf)
      Endif 

C Do not call aces_fin (use might be intending to run something
C else after external manupilations). 

CSSS      Call aces_fin
CSSS      Call C_free(Icore)

      Return
      End

