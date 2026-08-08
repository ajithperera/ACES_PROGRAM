










c     Perdew-Burke-Ernzerhof Correlation Potential
c
c     By S. Ivanov using the original version by Kieron Burke.
c
c     Phys. Rev. Lett. 77, 3865 (1996); ibid 78, 1396 (1997)
c     Phys. Rev. B 54, 16 533 (1996).
c     Phys. Rev. B 45, 13 244 (1992).

      subroutine pot_corr_pbe(pota,grcompa,potb,grcompb)

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
     &    pota,potb,grcompa,grcompb

      double precision
     &   con13,con23,con16,beta,gamma,delta,eta,pi

      parameter (con13=1.d0/3.d0,con23=2.d0/3.d0,con16=1.d0/6.d0)
      parameter (gamma=0.03109069086965489503494086371273d0)
      parameter (beta=0.06672455060314922d0,delta=beta/gamma)
      parameter (eta=1.d-12)

      double precision
     &   confk,consk,
     &   g,g3,sk,twoksg,t,t2,t4,b,b2,ec,fk,pon,
     &   q4,q5,h,ecrs,eczeta

      double precision
     &   g4,t6,gz,fac,bg,bec,q8,q9,hb,hrs,qln,
     &   hz,ht,comm

c g=phi(zeta)
c rs=(3/(4pi*rho))^(1/3)=local Seitz radius
c sk=Ks=Thomas-Fermi screening wavevector=sqrt(4fk/pi)
c twoksg=2*Ks*phi
c t=correlation dimensionless gradient=|grad rho|/(2*Ks*phi*rho)
c ec=lsd correlation energy
c h=gradient correction to correlation energy

 
      callstack_curr='POT_CORR_PBE'  

      pi = acos(-1.d0)

      if (ro.ge.thresh) then
        
        confk=(3.d0*pi*pi)**con13
        consk=dsqrt(4.d0*confk/pi)

        g=((1.d0+zeta)**con23+(1.d0-zeta)**con23)*0.5d0
        fk=confk*(ro**con13)
        sk=consk*(ro**con16)
        twoksg=2.d0*sk*g
        t=gro/(twoksg*ro)
     

c find LSD energy contributions

        call pot_corr_lsd2(pota,potb,ec,ecrs,eczeta)

c PBE correlation potential

        g3=g*g*g
        pon=-ec/(g3*gamma)
        b=delta/(dexp(pon)-1.d0)
        b2=b*b
        t2=t*t
        t4=t2*t2
        q4=1.d0+b*t2
        q5=1.d0+b*t2+b2*t4
        qln=1.d0+delta*q4*t2/q5
        h=g3*(beta/delta)*dlog(qln)

        g4=g3*g
        t6=t4*t2
        gz=(((1.d0+zeta)**2+eta)**(-con16)-
     &    ((1.d0-zeta)**2+eta)**(-con16))*con13

        fac=delta/b+1.d0 
        bg=-3.d0*b2*ec*fac/(beta*g4) 
        bec=b2*fac/(beta*g3) 
        q8=q5*q5+delta*q4*q5*t2 
        q9=1.d0+2.d0*b*t2 
        hb=-beta*g3*B*t6*(2.d0+b*t2)/q8

        hrs=hb*bec*ecrs

        hz=3.d0*gz*h/g+hb*(bg*gz+bec*eczeta)

        ht=2.d0*beta*g3*q9/q8

        comm=h-con13*rs*hrs-hz*zeta
     &    +(t2*ht*gz*zeta/g)-7.d0*con16*t2*ht

        pota=pota+comm+hz-(t2*ht*gz/g)
        potb=potb+comm-hz+(t2*ht*gz/g)

        comm=ro*ht*(t2/gro)
        grcompa=comm
        grcompb=comm
      else
        pota=0.d0
        potb=0.d0
        grcompa=0.d0
        grcompb=0.d0
      end if
      return
      end
