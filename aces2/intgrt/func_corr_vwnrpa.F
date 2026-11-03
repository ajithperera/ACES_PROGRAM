c Vosko, Wilk, Nusair correlation energy fitted to the RPA.
c
c    Ec(r,zeta) = e(p) + f(zeta) [e(f) - e(p)]
c    f(zeta)    = [(1+zeta)**4/3 + (1-zeta)**4/3 - 2]/(2**4/3 - 2)
c
c with e(p) and e(f) of the form [4.4] and the RPA parameters of the
c paper (Table 5, energies divided by 2 for Hartrees).  This is the
c VWN_RPA of libxc, which qtp17 uses; vwn (func_corr_vwn) is VWN5.
c
c S. H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys., 58, 1200, 1980.

      subroutine func_corr_vwnrpa(func)

      implicit none
#include <aces.h>
#include <dftfunc.com>
#include <callstack.com>

      M_REAL
     &    func

      M_REAL
     &    con43,con13,
     &    x0p,x0f,cp,cf,bp,bf,ap,af

      parameter(con43=4.d0/3.d0,con13=1.d0/3.d0,
     &    x0p=-0.409286d0,ap=0.03109070d0,  bp=13.0720d0,cp=42.7198d0,
     &    x0f=-0.743294d0,af=0.5d0*ap,      bf=20.1231d0,cf=101.578d0)

      M_REAL
     &    qi,xxi,b,c,x,x0,
     &    ei,a,q,xx,xx0,qp,qf,xs,xxp,xxf,xx0p,xx0f,ep,ef,fz

      qi(b,c)=dsqrt(4.d0*c-b*b)
      xxi(x,b,c)=x*x+b*x+c
      ei(a,x,x0,b,q,xx,xx0)=a*(dlog(x*x/xx) +
     &    (b+b)/q*datan(q/(x+x+b)) -
     &    b*x0/xx0* (dlog((x-x0)**2/xx)+
     &    (2.d0*(b+x0+x0)/q)*datan(q/(x+x+b))  )  )

      callstack_curr='FUNC_CORR_VWNRPA'

      func=0.d0
      if (ro.ge.thresh) then
         qp=qi(bp,cp)
         qf=qi(bf,cf)
         xs=dsqrt(rs)
         xxp=xxi(xs,bp,cp)
         xxf=xxi(xs,bf,cf)
         xx0p=xxi(x0p,bp,cp)
         xx0f=xxi(x0f,bf,cf)
         ep=ei(ap,xs,x0p,bp,qp,xxp,xx0p)
         ef=ei(af,xs,x0f,bf,qf,xxf,xx0f)

         fz=((1.d0+zeta)**con43+(1.d0-zeta)**con43-2.d0)/
     &      (2.d0**con43-2.d0)
         func=ro*(ep+fz*(ef-ep))
      end if

      return
      end
