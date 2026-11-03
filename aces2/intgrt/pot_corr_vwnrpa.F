c Vosko, Wilk, Nusair RPA correlation potential (see func_corr_vwnrpa).

      subroutine pot_corr_vwnrpa(pota,potb)

      implicit none
#include <aces.h>
#include <dftfunc.com>
#include <callstack.com>

      M_REAL
     &    pota,potb

      M_REAL
     &    con43,con13,
     &    x0p,x0f,cp,cf,bp,bf,ap,af

      parameter(con43=4.d0/3.d0,con13=1.d0/3.d0,
     &    x0p=-0.409286d0,ap=0.03109070d0,  bp=13.0720d0,cp=42.7198d0,
     &    x0f=-0.743294d0,af=0.5d0*ap,      bf=20.1231d0,cf=101.578d0)

      M_REAL
     &    qi,b,c,xxi,x,ei,a,x0,q,xx,xx0,delei,qp,qf,xxp,xxf,
     &    xx0p,xx0f,ep,ef,delep,delef,fz,delfz,delx,delzeta,delzetb,
     &    fact1,fact2,ec

      qi(b,c)=dsqrt(4.d0*c-b*b)
      xxi(x,b,c)=x*x+b*x+c
      ei(a,x,x0,b,q,xx,xx0)=a*(dlog(x*x/xx) +
     &    (b+b)/q*datan(q/(x+x+b)) -
     &    b*x0/xx0* (dlog((x-x0)**2/xx)+
     &    (2.d0*(b+x0+x0)/q)*datan(q/(x+x+b))  )  )
      delei(a,x,x0,b,q,xx,xx0)=a*(2/x-(2*x+b)/xx-4*b/((2*x+b)**2+q**2)-
     &    b*x0/xx0*(2.d0/(x-x0)-(2*x+b)/xx-
     &    4.d0*(2.d0*x0+b)/((2*x+b)**2+q**2)))

      callstack_curr='POT_CORR_VWNRPA'

      if (ro.ge.thresh) then

        qp=qi(bp,cp)
        qf=qi(bf,cf)
        x=sqrt(rs)
        xxp=xxi(x,bp,cp)
        xxf=xxi(x,bf,cf)
        xx0p=xxi(x0p,bp,cp)
        xx0f=xxi(x0f,bf,cf)
        ep=ei(ap,x,x0p,bp,qp,xxp,xx0p)
        ef=ei(af,x,x0f,bf,qf,xxf,xx0f)
        delep=delei(ap,x,x0p,bp,qp,xxp,xx0p)
        delef=delei(af,x,x0f,bf,qf,xxf,xx0f)

        fz=((1.d0+zeta)**con43+(1.d0-zeta)**con43-2.d0)/
     &     (2.d0**con43-2.d0)
        delfz=con43*((1.d0+zeta)**con13-(1.d0-zeta)**con13)/
     &        (2.d0**con43-2.d0)
c       d(x)/d(ro) and d(zeta)/d(roa), d(zeta)/d(rob)
        delx=-(x/ro)/6.d0
        delzeta=(1.d0-zeta)/ro
        delzetb=-(1.d0+zeta)/ro

        fact1=delx*(delep+fz*(delef-delep))
        fact2=delfz*(ef-ep)
        ec=ep+fz*(ef-ep)

        pota=ec+ro*(fact1+delzeta*fact2)
        potb=ec+ro*(fact1+delzetb*fact2)
      else
        pota=0.d0
        potb=0.d0
      end if

      return
      end
