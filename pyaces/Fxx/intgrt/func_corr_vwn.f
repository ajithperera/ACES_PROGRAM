










c Vosko, Wilk, Nusair correlation energy.
c
c Written by S. Beck; Modified by S. Ivanov
c
c
c    [4.4] defines e(p), e(f), and e(a)
c    [2.4] Ec(r,zeta)     = e(p) + del ec(r,zeta)
c    [4.7] del ec(r,zeta) = ...
c          del ec(r,1)/del eRPA(r,1) ~ 1   (see table 5)
c    [3.2] del eRPA(r,zeta) = ...
c          alphaRPA == e(a)
c    [2.3] f(zeta) =
c          f''(0) = 8/9 [2**4/3 - 2]**-1
c          betaRPA(r) = beta1(r)
c    [4.8] beta1(r) =
c          del ec(r,1) = e(f) - e(p)   (see table 5)
c
c Note that the energy units in the paper are in Rydberts.  To switch
c to Hartrees, the values of A are all divided by 2.
c
c S. H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys., 58, 1200, 1980.
c
      
           subroutine func_corr_vwn(func)

      implicit none







c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler























cYAU - ACES3 stuff . . . we hope - #include <aces.par>






c The following values are used in determining the values of the various
c functionals:
c
c    roa    : The alpha density
c    rob    :
c    ro     : roa+rob
c    roinv  : 1/ro
c    rom    : roa-rob
c    rs     : (3 / 4 pi ro)**(1/3)
c    zeta   : (roa-rob)/ro
c
c    gradx  : d(roa)/dx + d(rob)/dx
c    grady  :
c    gradz  :
c    gradxm : d(roa)/dx - d(rob)/dx
c    gradym :
c    gradzm :
c    gro2   : gradx**2 + grady**2 + gradz**2
c    gro    : sqrt(gro2)
c    gro2a  : [d(roa)/dx]**2 + [d(roa)/dy]**2 + [d(roa)/dz]**2
c    gro2b  :
c
c    hesxx  : d2(roa)/dx2 + d2(rob)/dx2
c    hesyy  :
c    heszz  :
c    hesxy  :
c    hesxz  :
c    hesyz  :
c    hesxxm : d2(roa)/dx2 - d2(rob)/dx2
c    hesyym :
c    heszzm :
c    hesxym :
c    hesxzm :
c    hesyzm :
c    xlap   : hesxx + hesyy + heszz
c    xlapm  : hesxxm + hesyym + heszzm
c    trm1   :
c    trm1m  :
c    trm2   :
c    trm2m  :
c    grdaa  : gradroa.gradroa
c    grdbb  :
c    grdab  :
c    thresh : Parameter which tells what is the cutoff value
c           :  for is for the density, functionals and potetentials.
c           : The flag is cutoff and the default value is 10**(-12).
      double precision
     &    roa,rob,ro,rom,rs,zeta,roinv,
     &    gradx,grady,gradz,gradxm,gradym,gradzm,gro,gro2,gro2a,gro2b,
     &    hesxx,hesyy,heszz,hesxy,hesxz,hesyz,
     &    hesxxm,hesyym,heszzm,hesxym,hesxzm,hesyzm,
     &    xlap,xlapm,trm1,trm1m,trm2,trm2m,
     &    grdaa,grdbb,grdab,thresh

c There are several types of exchange and correlation energy functionals
c as well as hybrid schemes.  There are also corresponding potentials
c for use in Kohn-Sham.  Define the number of each type of energy functional
c and parameters to make calling them clearer.
c
c We want to have freedom to sellect potential and functional for
c better flexibility. Usually, the choice of functional determines 
c the poential  to be used in the  KS calculations.
c 
c
c fun_num_exch  : the number of exchange energy functionals (or potentials)
c fun_exch      : the total exchange energy contribution from the current
c                 atom
c tot_exch      : the total exchange energy contribution from all atoms
c coef_exch     : the coefficient of each exchange functional
c ene_exch      : the total energy with this functional
c nam_exch      : the name of the functional
c coef_pot_exch : the coefficient of each exchange component in the KS potential
c coef_pot_corr : the coefficient of each correlation component in the KS potential

      integer fun_num_exch,fun_num_corr,fun_num_hyb
      parameter (fun_num_exch = 5)
      parameter (fun_num_corr = 6)
      parameter (fun_num_hyb  = 1)

      double precision
     &    fun_exch(fun_num_exch),fun_corr(fun_num_corr),
     &    fun_hyb(fun_num_hyb),
     &    tot_exch(fun_num_exch),tot_corr(fun_num_corr),
     &    tot_hyb(fun_num_hyb),
     &    ene_exch(fun_num_exch),ene_corr(fun_num_corr),
     &    ene_hyb(fun_num_hyb),
     &    coef_exch(fun_num_exch),coef_corr(fun_num_corr),
     &    coef_pot_exch(fun_num_exch),coef_pot_corr(fun_num_corr)

      double precision vxc_ksalpha,vxc_ksbeta

      integer fun_exch_none,fun_corr_none,fun_special,
     &    fun_dft_none,
     &    fun_exch_lda,fun_exch_becke,fun_exch_pbe,fun_exch_pw91,
     &    fun_exch_hf,
     &    fun_corr_vwn,fun_corr_lyp,fun_corr_pbe,
     &    fun_corr_pw91, fun_corr_wl,fun_corr_wi,
     &    fun_hyb_b3lyp

      character*50
     &    nam_exch(fun_num_exch),nam_corr(fun_num_corr),
     &    nam_hyb(fun_num_hyb)
      character*80 nam_func,nam_kspot
      character*10
     &    abb_exch(fun_num_exch),abb_corr(fun_num_corr),
     &    abb_hyb(fun_num_hyb)

      parameter (fun_exch_none = 0)
      parameter (fun_exch_lda  = 1)
      parameter (fun_exch_becke= 2)
      parameter (fun_exch_pbe  = 3)
      parameter (fun_exch_pw91 = 4)
      parameter (fun_exch_hf   = 5)

      parameter (fun_corr_none = 0)
      parameter (fun_corr_vwn  = 1)
      parameter (fun_corr_lyp  = 2)
      parameter (fun_corr_pbe  = 3)
      parameter (fun_corr_pw91 = 4)
      parameter (fun_corr_wl   = 5)
      parameter (fun_corr_wi   = 6)

      parameter (fun_special   = 0)

      parameter (fun_dft_none  = -1)
      parameter (fun_hyb_b3lyp = 1)

      common /dftfunc/ fun_exch,fun_corr,fun_hyb,tot_exch,tot_corr,
     &    tot_hyb,ene_exch,ene_corr,ene_hyb,coef_exch,coef_corr,
     &    roa,rob,ro,rom,rs,zeta,roinv,
     &    gradx,grady,gradz,gradxm,gradym,gradzm,gro,gro2,gro2a,gro2b,
     &    hesxx,hesyy,heszz,hesxy,hesxz,hesyz,
     &    hesxxm,hesyym,heszzm,hesxym,hesxzm,hesyzm,
     &    xlap,xlapm,trm1,trm1m,trm2,trm2m,
     &    grdaa,grdbb,grdab,thresh,vxc_ksalpha,vxc_ksbeta,
     &    coef_pot_exch,coef_pot_corr

      common /dftfuncc/ nam_exch,nam_corr,nam_hyb,nam_func,
     &    abb_exch,abb_corr,abb_hyb,nam_kspot

      save /dftfunc/
      save /dftfuncc/

c Parameters for hybrid methods:
      double precision
     &    b3lypa,b3lypb,b3lypc
      parameter (b3lypa=0.20d0)
      parameter (b3lypb=0.72d0)
      parameter (b3lypc=0.81d0)




c This contains the global string for identifying the current subroutine
c or function (provided the programmer set it).  cf. tools/callstack.F
c BE GOOD AND RESET CURR ON EXIT!

      character*64                callstack_curr,callstack_prev
      common /callstack_curr_com/ callstack_curr,callstack_prev
      save   /callstack_curr_com/



      double precision
     &    func

      double precision
     &    con43,con13,con98,pi,
     &    x0p,x0f,x0a,cp,cf,ca,bp,bf,ba,ap,af,aa

      parameter(con43=4.d0/3.d0,con13=1.d0/3.d0,con98=9.d0/8.d0,
     &    x0p=-0.10498d0, ap=0.03109070d0,   bp=3.72744d0, cp=12.9352d0,
     &    x0f=-0.325d0,   af=0.5d0*ap,       bf=7.06042d0, cf=18.0578d0,
     &    x0a=-0.0047584d0,ba=1.13107d0,
     &    ca=13.0045d0)

      double precision
     &    qi,xxi,b,c,x,x0,g1,g2,
     &    ei,a,q,xx,xx0,qp,qf,qa,xs,xxp,xxf,
     &    xxa,xx0p,xx0f,xx0a,ep,ef,ea

      qi(b,c)=dsqrt(4.d0*c-b*b)
      xxi(x,b,c)=x*x+b*x+c
      ei(a,x,x0,b,q,xx,xx0)=a*(dlog(x*x/xx) +
     &    (b+b)/q*datan(q/(x+x+b)) -
     &    b*x0/xx0* (dlog((x-x0)**2/xx)+
     &    (2.d0*(b+x0+x0)/q)*datan(q/(x+x+b))  )  )

      callstack_curr='FUNC_CORR_VWN'

      pi = acos(-1.d0)
      aa = -1.d0/(6.d0*pi*pi)

      func=0.d0
      if (ro.ge.thresh) then 
         qp=qi(bp,cp)
         qf=qi(bf,cf)
         qa=qi(ba,ca)
         xs=dsqrt(rs)
         xxp=xxi(xs,bp,cp)
         xxf=xxi(xs,bf,cf)
         xxa=xxi(xs,ba,ca)
         xx0p=xxi(x0p,bp,cp)
         xx0f=xxi(x0f,bf,cf)
         xx0a=xxi(x0a,ba,ca)
         ep=ei(ap,xs,x0p,bp,qp,xxp,xx0p)
         ef=ei(af,xs,x0f,bf,qf,xxf,xx0f)
         ea=ei(aa,xs,x0a,ba,qa,xxa,xx0a)

         g1=4.d0/(9.d0*(2.d0**con13-1.d0))
         g2=con98*((1.d0+zeta)**con43+(1.d0-zeta)**con43-2.d0)
         func=1.d0+(g1*(ef-ep)/ea-1.d0)*zeta**4

         func=ro*(ep+ea*g2*func)
      end if

      return
      end

