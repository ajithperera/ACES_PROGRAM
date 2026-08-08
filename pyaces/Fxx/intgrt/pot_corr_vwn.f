











      subroutine pot_corr_vwn(pota,potb)

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
     &    pota,potb

      double precision
     &    con43,con13,con98,con49,con32,pi,
     &    x0p,x0f,x0a,cp,cf,ca,bp,bf,ba,ap,af,aa

      parameter(con43=4.d0/3.d0,con13=1.d0/3.d0,con98=9.d0/8.d0,
     &    con49=4.d0/9.d0,con32=3.d0/2.d0,
     &    x0p=-0.10498d0, ap=0.03109070d0,   bp=3.72744d0, cp=12.9352d0,
     &    x0f=-0.325d0,   af=0.5d0*ap,       bf=7.06042d0, cf=18.0578d0,
     &    x0a=-0.0047584d0,ba=1.13107d0,
     &    ca=13.0045d0)

      double precision
     &    qi,b,c,xxi,x,ei,a,x0,q,xx,xx0,delei,qp,qf,qa,xxp,xxf,xxa,
     &    xx0p,xx0f,xx0a,ep,ef,ea,delep,delef,delea,con,
     &    g,h,delg,delh,delx,delzeta,delzetb,deleca,delecb,fact1,fact2,
     &    zeta4,ec

      qi(b,c)=dsqrt(4.d0*c-b*b)
      xxi(x,b,c)=x*x+b*x+c
      ei(a,x,x0,b,q,xx,xx0)=a*(dlog(x*x/xx) +
     &    (b+b)/q*datan(q/(x+x+b)) -
     &    b*x0/xx0* (dlog((x-x0)**2/xx)+
     &    (2.d0*(b+x0+x0)/q)*datan(q/(x+x+b))  )  )
      delei(a,x,x0,b,q,xx,xx0)=a*(2/x-(2*x+b)/xx-4*b/((2*x+b)**2+q**2)-
     &    b*x0/xx0*(2.d0/(x-x0)-(2*x+b)/xx-
     &    4.d0*(2.d0*x0+b)/((2*x+b)**2+q**2)))

      callstack_curr='POT_CORR_VWN'

      pi = acos(-1.d0)
      aa = -1.d0/(6.d0*pi*pi)

      if (ro.ge.thresh) then

        qp=qi(bp,cp)
        qf=qi(bf,cf)
        qa=qi(ba,ca)
        x=sqrt(rs)
        xxp=xxi(x,bp,cp)
        xxf=xxi(x,bf,cf)
        xxa=xxi(x,ba,ca)
        xx0p=xxi(x0p,bp,cp)
        xx0f=xxi(x0f,bf,cf)
        xx0a=xxi(x0a,ba,ca)
        ep=ei(ap,x,x0p,bp,qp,xxp,xx0p)
        ef=ei(af,x,x0f,bf,qf,xxf,xx0f)
        ea=ei(aa,x,x0a,ba,qa,xxa,xx0a)
        delep=delei(ap,x,x0p,bp,qp,xxp,xx0p)
        delef=delei(af,x,x0f,bf,qf,xxf,xx0f)
        delea=delei(aa,x,x0a,ba,qa,xxa,xx0a)

        con=con49/(2.d0**con13-1.d0)
        g=con98*((1.d0+zeta)**con43+(1.d0-zeta)**con43-2.d0)
        delg=con32*((1.d0+zeta)**con13-(1.d0-zeta)**con13)
        h=con*(ef-ep)/ea-1.d0
        delh=con/ea*(delef-delep-(ef-ep)*delea/ea)
        delx=-(x/ro)/6.d0
        delzeta=(1.d0-zeta)/ro
        delzetb=-(1.d0+zeta)/ro
        zeta4=zeta**4

        fact1=delx*(delep+delea*g*(1.d0+h*zeta4)+ea*g*delh*zeta4)
        fact2=ea*(delg*(1.d0+h*zeta4)+4.d0*g*h*zeta**3)
        deleca=fact1+delzeta*fact2
        delecb=fact1+delzetb*fact2
        ec=ep+ea*g*(1.d0+h*zeta4)

        pota=ec+ro*deleca
        potb=ec+ro*delecb
      else
        pota=0.d0
        potb=0.d0
      end if

      return
      end
