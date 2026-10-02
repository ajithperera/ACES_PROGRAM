      subroutine oep_mem_init(maxmem, irppmem, imemget, icrsiz)
c 
      implicit double precision (a-h, o-z)
c
#include <machsp.com>
#include <flags.com>
#include <flags2.com>
      integer icrsiz
c     Never ask for more than the core still available (ICRSIZ); round
c     down to even so IMEMGET stays <= ICRSIZ.
      maxmem  = min(iflags (36), icrsiz) - irppmem
     &        - iflags(44)*iflags(37)
      maxmem  = maxmem      - mod (maxmem, 2)
      imemget = maxmem      + irppmem + iflags(44)*iflags(37)
      imemget = imemget     - mod (imemget, 2)
     
      Return
      End

