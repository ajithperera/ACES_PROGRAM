c Short-range Becke exchange of one spin density
c
c The Becke exchange energy density of spin density ro,
c
c    eb = - ro**(4/3) K(x),  K(x) = Cx + beta x**2/(1+6 beta x asinh(x))
c
c with x = |grad ro|/ro**(4/3) and Cx the LDA constant, is attenuated as
c in Iikura, Tsuneda, Yanai and Hirao, JCP 115, 3540 (2001): eb is
c taken as the exchange of a uniform gas with local Fermi wavevector
c
c    k = sqrt(9 pi/(2 K)) ro**(1/3)
c
c (k is the Fermi wavevector itself when K=Cx), giving
c
c    e = eb F(a),  a = omega/(2k)
c
c with F from ityh_attn and omega = int_rs_omega.
c
c    ro   : the spin density (ro >= thresh)
c    gro2 : |grad ro|**2
c    e    : the energy density
c    vro  : de/dro               (only if doder)
c    vgro : de/d|grad ro|        (only if doder)

      subroutine becke_sr_spin(ro1,gro21,e,vro,vgro,doder)

      implicit none
#include <aces.h>
#include <intgrtflags.com>

      M_REAL
     &    ro1,gro21,e,vro,vgro
      logical doder

      M_REAL
     &    con13,con43,con98,alpha,beta
      parameter (con13=1.d0/3.d0,con43=4.d0/3.d0,con98=9.d0/8.d0)
      parameter (alpha=2.d0/3.d0)
      parameter (beta=0.0042d0)

      M_REAL
     &    pi,conlda,x,x2,ash,d,ck,dck,ro13,ro43,a,f,df,eb,
     &    dxdro,dxdgro,dadro,dadgro

      pi=acos(-1.d0)
      conlda=alpha*(2.d0**con13)*con98*((3.d0/pi)**con13)

      ro13=ro1**con13
      ro43=ro1*ro13
      x2=gro21/(ro43*ro43)
      x=sqrt(x2)
      ash=log(x+sqrt(x2+1.d0))
      d=1.d0+6.d0*beta*x*ash
      ck=conlda+beta*x2/d
      eb=-ro43*ck

      a=0.5d0*int_rs_omega*sqrt(2.d0*ck/(9.d0*pi))/ro13
      call ityh_attn(a,f,df)
      e=eb*f

      if (doder) then
c        dK/dx
         dck=beta*(2.d0*x*d-6.d0*beta*x2*(ash+x/sqrt(x2+1.d0)))/(d*d)
         dxdro=-con43*x/ro1
         dxdgro=1.d0/ro43
         dadro=a*(0.5d0*dck/ck*dxdro-con13/ro1)
         dadgro=a*0.5d0*dck/ck*dxdgro
         vro=-con43*ro13*(ck-x*dck)*f+eb*df*dadro
         vgro=-dck*f+eb*df*dadgro
      end if

      return
      end
