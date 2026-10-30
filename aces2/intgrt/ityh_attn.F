c Attenuation function of short-range LDA-type exchange
c
c The exchange energy density of a uniform gas with Fermi wavevector k
c changes by the factor F(a), a = omega/(2k), when the interaction
c 1/r12 is replaced by erfc(omega*r12)/r12:
c
c    F(a) = 1 - (8/3) a [ sqrt(pi) erf(1/(2a)) + (2a - 4a**3) E
c                         - 3a + 4a**3 ],     E = exp(-1/(4a**2))
c
c    dF/da = -(8/3) [ B(a) + a (12 a**2 (1-E) - 3) ]
c
c where B(a) is the bracket above.  For large a the bracket is a
c difference of large numbers, so there F is summed from its expansion
c in u = 1/a**2 instead,
c
c    F = u/36 - u**2/960 + u**3/26880 - u**4/829440 + ...
c
c Used for becke_sr in the form of Iikura, Tsuneda, Yanai and Hirao,
c JCP 115, 3540 (2001).
c
c    a  : the argument (a > 0)
c    f  : F(a)
c    df : dF/da

      subroutine ityh_attn(a,f,df)

      implicit none
#include <aces.h>

      M_REAL
     &    a,f,df

      integer nser,i
      parameter (nser=8)
      M_REAL
     &    con83,aser,pi,e,b,u,up,cser(nser)
      parameter (con83=8.d0/3.d0)
      parameter (aser=2.d0)
c     cser(i) = 1/36, -1/960, 1/26880, -1/829440, 1/28385280,
c               -1/1073479680, 1/44590694400, -1/2021444812800
      data cser /2.7777777777777778d-02,-1.0416666666666667d-03,
     &           3.7202380952380952d-05,-1.2056327160493827d-06,
     &           3.5229527417027417d-08,-9.3155000381562881d-10,
     &           2.2426203795561436d-11,-4.9469567196091402d-13/

      if (a.gt.aser) then
         u=1.d0/(a*a)
         f=0.d0
         df=0.d0
         up=1.d0
         do i=1,nser
c           up = u**(i-1)
            f=f+cser(i)*up*u
            df=df+dble(i)*cser(i)*up
            up=up*u
         end do
c        dF/da = dF/du * du/da, du/da = -2/a**3
         df=-2.d0*df*u/a
      else
         pi=acos(-1.d0)
         e=exp(-1.d0/(4.d0*a*a))
         b=sqrt(pi)*erf(1.d0/(2.d0*a))+(2.d0*a-4.d0*a**3)*e
     &     -3.d0*a+4.d0*a**3
         f=1.d0-con83*a*b
         df=-con83*(b+a*(12.d0*a*a*(1.d0-e)-3.d0))
      end if

      return
      end
