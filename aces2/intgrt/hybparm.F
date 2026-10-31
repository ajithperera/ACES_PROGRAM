c The hybrid functionals known by name (FUNC=name, KSPOT=name).
c
c For hybrid ihyb (1..fun_num_hyb) this returns its abbreviation (the
c name in the input), its full name, the coefficients of the exchange
c and correlation functionals it is made of, and the range-separation
c parameter omega (bohr**-1, 0 if the hybrid is not range-separated).
c
c Range-separated hybrids split the exchange operator 1/r12 as
c    [alpha + beta*erf(omega*r12)]/r12      : exact exchange
c                                             (alpha hf + beta hf_lr)
c    [1 - alpha - beta*erf(omega*r12)]/r12  : Becke exchange
c                                             ((1-alpha-beta) becke
c                                              + beta becke_sr)
c
c    b3lyp      Becke III LYP
c    cam-b3lyp  Yanai, Tew and Handy, CPL 393, 51 (2004)
c    cam-qtp-00 Verma and Bartlett, JCP 140, 18A534 (2014)
c    cam-qtp-01 Jin and Bartlett, JCP 145, 034107 (2016)
c    cam-qtp-02 Haiduke and Bartlett, JCP 148, 184106 (2018)
c    lc-qtp     Haiduke and Bartlett, JCP 148, 184106 (2018)
c    qtp17      Jin and Bartlett, JCP 149, 064111 (2018)
c
c This table is the only place these values are kept: INITINTGRT,
c the hybrid energies in NUMINTEFF/NUMINTINT and RSOMEGA (omega for
c the integral programs VMOL and XERDINTS) all read it.

      subroutine hybparm(ihyb,abb,nam,cx,cc,omega)

      implicit none
#include <aces.h>
#include <dftfunc.com>

      integer ihyb
      character*(*) abb,nam
      M_REAL
     &    cx(fun_num_exch),cc(fun_num_corr),omega

      integer i
      M_REAL
     &    alpha,beta

      do i=1,fun_num_exch
        cx(i)=0.d0
      end do
      do i=1,fun_num_corr
        cc(i)=0.d0
      end do
      omega=0.d0
      alpha=0.d0
      beta=0.d0

      if (ihyb.eq.fun_hyb_b3lyp) then
         abb='b3lyp'
         nam='Becke III LYP'
         cx(fun_exch_hf)=b3lypa
         cx(fun_exch_lda)=1.d0-b3lypa-b3lypb
         cx(fun_exch_becke)=b3lypb
         cc(fun_corr_vwn)=1.d0-b3lypc
         cc(fun_corr_lyp)=b3lypc
      else if (ihyb.eq.fun_hyb_camb3lyp) then
         abb='cam-b3lyp'
         nam='CAM-B3LYP'
         alpha=0.19d0
         beta=0.46d0
         omega=0.33d0
         cc(fun_corr_vwn)=0.19d0
         cc(fun_corr_lyp)=0.81d0
      else if (ihyb.eq.fun_hyb_camqtp00) then
         abb='cam-qtp-00'
         nam='CAM-QTP-00'
         alpha=0.54d0
         beta=0.37d0
         omega=0.29d0
         cc(fun_corr_vwn)=0.20d0
         cc(fun_corr_lyp)=0.80d0
      else if (ihyb.eq.fun_hyb_camqtp01) then
         abb='cam-qtp-01'
         nam='CAM-QTP-01'
         alpha=0.23d0
         beta=0.77d0
         omega=0.31d0
         cc(fun_corr_vwn)=0.20d0
         cc(fun_corr_lyp)=0.80d0
      else if (ihyb.eq.fun_hyb_camqtp02) then
         abb='cam-qtp-02'
         nam='CAM-QTP-02'
         alpha=0.28d0
         beta=0.72d0
         omega=0.335d0
         cc(fun_corr_lyp)=1.d0
      else if (ihyb.eq.fun_hyb_lcqtp) then
         abb='lc-qtp'
         nam='LC-QTP'
         alpha=0.d0
         beta=1.d0
         omega=0.475d0
         cc(fun_corr_lyp)=1.d0
      else if (ihyb.eq.fun_hyb_qtp17) then
         abb='qtp17'
         nam='QTP17'
         cx(fun_exch_hf)=0.62d0
         cx(fun_exch_lda)=0.38d0
         cc(fun_corr_vwnrpa)=0.20d0
         cc(fun_corr_lyp)=0.80d0
      else
         write(*,*) '@HYBPARM: no hybrid functional ',ihyb
         call errex
      end if

      if (omega.ne.0.d0) then
         cx(fun_exch_hf)=alpha
         cx(fun_exch_hflr)=beta
         cx(fun_exch_becke)=1.d0-alpha-beta
         cx(fun_exch_becke_sr)=beta
      end if

      return
      end
