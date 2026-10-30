c Short-range Becke exchange potential (see becke_sr_spin)
c
c pota is de/d(roa) and grcompa is de/d|grad roa|, as in pot_exch_becke.

      subroutine pot_exch_becke_sr(pota,grcompa,potb,grcompb)

      implicit none
#include <aces.h>
#include <iuhf.com>
#include <dftfunc.com>
#include <callstack.com>

      M_REAL
     &    pota,grcompa,potb,grcompb

      M_REAL
     &    e

      callstack_curr='POT_EXCH_BECKE_SR'

      pota=0.d0
      grcompa=0.d0
      if (roa.ge.thresh)
     &   call becke_sr_spin(roa,gro2a,e,pota,grcompa,.true.)

      if (iuhf.eq.0) then
        potb=pota
        grcompb=grcompa
      else
        potb=0.d0
        grcompb=0.d0
        if (rob.ge.thresh)
     &     call becke_sr_spin(rob,gro2b,e,potb,grcompb,.true.)
      end if

      return
      end
