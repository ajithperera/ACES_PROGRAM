










c     Perdew-Wang 91 Correlation Potential
c
c     By S. Ivanov using the original version by Kieron Burke.
c
c     Phys. Rev. B 54, 16 533 (1996).
c     Phys. Rev. B 45, 13 244 (1992).


      subroutine pot_corr_pw91(pota,grcompa,potb,grcompb)

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
     &   con13,con23,con16,pi

      parameter (con13=1.d0/3.d0,con23=2.d0/3.d0,con16=1.d0/6.d0)

      double precision
     &   xnu,cc0,cx,
     &   c1,c2,c3,c4,c5,c6,a4,
     &   alfa,beta,delta,eta

      parameter (xnu=15.75592d0,cc0=0.004235d0,cx=-0.001667212d0)
      parameter (c1=0.002568d0,c2=0.023266d0,c3=7.389d-6,c4=8.723d0,
     &           c5=0.472d0,c6=7.389d-2,a4=100.0d0)
      parameter(alfa=0.090d0,beta=xnu*cc0,delta=2.d0*alfa/beta)
      parameter (eta=1.d-12)

      double precision
     &   confk,consk,
     &   g,g3,g4,sk,twoksg,t,t2,t4,t6,b,b2,ec,fk,pon,
     &   rs2,rs3,q4,q5,q6,q7,r0,r1,r2,r3,r4,cc,coeff,h0,h1,h

      double precision
     &   ecrs,eczeta,gz,ccrs,rs13,bg,bec,fac,q8,q9,
     &   h0b,h0z,h0rs,h0t,h1z,h1rs,h1t,hz,hrs,ht,comm

c g=phi(zeta)
c rs=(3/(4pi*rho))^(1/3)=local Seitz radius
c sk=Ks=Thomas-Fermi screening wavevector=sqrt(4fk/pi)
c twoksg=2*Ks*phi
c t=correlation dimensionless gradient=|grad rho|/(2*Ks*phi*rho)
c ec=lsd correlation energy
c h=gradient correction to correlation energy

      callstack_curr='POT_CORR_PW91'  

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

c PW91 correlation potential

        g3=g*g*g
        g4=g3*g
        pon=-ec*delta/(g3*beta)

        b=delta/(dexp(pon)-1.d0)
        b2=b*b

        t2=t*t
        t4=t2*t2

        rs2=rs*rs
        rs3=rs2*rs

        q4=1.d0+b*t2
        q5=1.d0+b*t2+b2*t4
        q6=c1+c2*rs+c3*rs2
        q7=1.d0+c4*rs+c5*rs2+c6*rs3

        cc=-cx+q6/q7
        coeff=cc-cc0-3.d0*cx/7.d0

        r0=0.663436444d0*rs
        r1=a4*r0*g4
        r2=xnu*coeff*g3
        r3=dexp(-r1*t2)

        h0=g3*(beta/delta)*dlog(1.d0+delta*q4*t2/q5)
        h1=r3*r2*t2
        h=h0+h1

        t6=t4*t2

        gz=(((1.d0+zeta)**2+eta)**(-con16)-
     &   ((1.d0-zeta)**2+eta)**(-con16))/3.d0

        ccrs=(c2+2.d0*c3*rs)/q7-q6*(c4+2.d0*c5*rs+3.d0*c6*rs2)/q7**2
        rs13=rs*con13
        r4 =rs13*ccrs/coeff
        fac=delta/b+1.d0
        bg=-3.d0*b2*ec*fac/(beta*g4)
        bec=b2*fac/(beta*g3)
        q8=q5*q5+delta*q4*q5*t2
        q9=1.d0+2.d0*b*t2

        h0b=-beta*g3*b*t6*(2.d0+b*t2)/q8
        h0rs=-rs13*ecrs*h0b*bec

        h0z=3.d0*gz*h0/g+h0b*(bg*gz+bec*eczeta)
        h0t=2.d0*beta*g3*q9/q8

        h1rs=r3*r2*t2*(-r4+r1*t2*con13)
        h1z=gz*r3*r2*t2*(3.d0-4.d0*r1*t2)/g
        h1t=2.d0*r3*r2*(1.d0-r1*t2)

        hrs=h0rs+h1rs
        hz=h0z+h1z
        ht=h0t+h1t


        comm=h+hrs-hz*zeta
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
