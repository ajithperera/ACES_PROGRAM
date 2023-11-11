#include "flags.h"

      Subroutine aces2_main()
      implicit none

      integer iuhf 
c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"

c ----------------------------------------------------------------------
         Call Joda_()
         Call Aces_ja_fin
         Call aces_init(icore,i0,icrsiz, iuhf, .true.)
CSSS         Call Test_py(icore(i0),Icrsiz,Scr,Ndim,Iuhf)
c ----------------------------------------------------------------------
      Return
      End 

