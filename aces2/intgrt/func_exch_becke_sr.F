c Short-range Becke exchange functional (see becke_sr_spin)

      subroutine func_exch_becke_sr(func)

      implicit none
#include <aces.h>
#include <dftfunc.com>
#include <iuhf.com>
#include <callstack.com>

      M_REAL
     &    func

      M_REAL
     &    funa,funb,dum1,dum2

      callstack_curr='FUNC_EXCH_BECKE_SR'

      funa=0.d0
      if (roa.ge.thresh)
     &   call becke_sr_spin(roa,gro2a,funa,dum1,dum2,.false.)

      if (iuhf.eq.0) then
        funb=funa
      else
        funb=0.d0
        if (rob.ge.thresh)
     &     call becke_sr_spin(rob,gro2b,funb,dum1,dum2,.false.)
      end if

      func=funa+funb

      return
      end
