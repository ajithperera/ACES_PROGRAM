      Subroutine Runints(Iuhf)

      Implicit None
      Integer*8 Iuhf

#include "icore.com"
#include "istart.com"
    
      Call Ints_(Icore(i0),Icrsiz,Iuhf)
      Call V2ja(icore(i0),Icrsiz,Iuhf)

      Return
      End

