#include "flags.h"

      SUBROUTINE HBAR_(ICORE,ICRSIZ,IUHF)
      implicit none

      integer icore(icrsiz),icrsiz,iuhf,i0

c COMMON BLOCKS
#include "flags.com"
#include "flags2.com"
#include "machsp.com"

c ----------------------------------------------------------------------
      I0 = 1
C
      Call Drive_Hbar(Icore(i0), icrsiz/iintfp, Iuhf)

C
c ----------------------------------------------------------------------
      RETURN
      end

