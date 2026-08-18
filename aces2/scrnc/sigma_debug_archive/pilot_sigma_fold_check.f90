! Pilot code: independent ground-truth check of the two-propagator
! frequency-folding fix applied to scrnc_sigma_driver.F (see
! SCRNC_DEBUG_LOG.md "continued 16"). No ACES2/chemistry involved --
! this validates ONLY the quadrature/folding arithmetic itself, using
! a synthetic analytic surrogate for Wc(i*omega') so an independent
! numerical integral can serve as ground truth.
!
! Model:  Sigma_c(i*omega) = -(1/2pi) Int_{-inf}^{inf} dw'
!                              [1/(i*omega+i*w'-eps_m)] * g(w')
! with g(w') = 1/(1+(w'/b)**2)  (even, smooth, decaying -- a stand-in
! for Wc_reduced(p,q,m; i*omega') at fixed p=q=m).
!
! Ground truth: substitute w' = b*tan(theta), theta in (-pi/2,pi/2),
! which maps g(w')*dw' -> b*dtheta exactly (removes the 1/w'^2 tail
! analytically), then integrate the resulting smooth finite-interval
! integrand with a fine composite Simpson rule. This is a DIFFERENT
! substitution and a DIFFERENT quadrature rule than the production
! code's Gauss-Legendre + rational map -- genuinely independent check,
! not the code testing itself.
!
! Two candidate folded sums are built with the SAME grid-building
! routines as production (Scrnc_sigma_freq_grid/gaussl, copied
! verbatim below, same duplicated-per-module convention as the real
! file):
!   ONE-TERM (pre-fix, buggy):  only the (omega+omega') propagator
!   TWO-TERM (post-fix):        both (omega+omega') and (omega-omega')
!
! Expect: TWO-TERM converges to the Simpson ground truth as N grows;
! ONE-TERM does not (converges to something else, i.e. a real bug).

program pilot_sigma_fold_check
  implicit double precision (a-h,o-z)
  parameter (maxfreq=200)
  dimension xnode(maxfreq), wnode(maxfreq)
  double precision, dimension(3) :: eps_list
  double precision, dimension(3) :: omega_list
  double precision :: b, omega_c, pi

  pi = 4.0d0*datan(1.0d0)
  b  = 1.5d0
  omega_c = 1.0d0

  eps_list    = (/ -0.5d0, 0.3d0, -20.6d0 /)
  omega_list  = (/  0.0d0, 0.3d0,   0.75d0 /)

  write(6,'(a)') '=== Pilot check: negative-frequency folding fix ==='
  write(6,'(a)') 'g(w'') = 1/(1+(w''/b)**2), b=1.5, Omega_c=1.0'
  write(6,'(a)') ''

  do ie = 1, 3
    eps_m = eps_list(ie)
    do io = 1, 3
      omega = omega_list(io)

      call exact_ref(omega, eps_m, b, sre_exact, sim_exact)

      write(6,'(a,1PD12.4,a,1PD12.4)') &
        'eps_m=', eps_m, '   omega=', omega
      write(6,'(a,2(1x,1PD16.8))') &
        '  ground truth   (Re,Im) = ', sre_exact, sim_exact

      do n = 8, 32, 8
        call scrnc_sigma_freq_grid(n, omega_c, xnode, wnode, maxfreq)
        call folded_sum(n, xnode, wnode, omega, eps_m, b, .true., &
                         sre1, sim1)
        call folded_sum(n, xnode, wnode, omega, eps_m, b, .false., &
                         sre2, sim2)
        write(6,'(a,I3,a,2(1x,1PD16.8),a,2(1x,1PD16.8))') &
          '  N=', n, &
          '  one-term(Re,Im)=', sre1, sim1, &
          '   two-term(Re,Im)=', sre2, sim2
      end do
      write(6,'(a,2(1x,1PD10.2))') &
        '  |two-term - exact| at N=32 = ', &
        dabs(sre2-sre_exact), dabs(sim2-sim_exact)
      write(6,'(a,2(1x,1PD10.2))') &
        '  |one-term - exact| at N=32 = ', &
        dabs(sre1-sre_exact), dabs(sim1-sim_exact)
      write(6,'(a)') ''
    end do
  end do

contains

  ! Folded sum, mirrors scrnc_sigma_driver.F's fold loop exactly
  ! (one_term=.true. -> old buggy behavior, .false. -> fixed).
  subroutine folded_sum(n, xnode, wnode, omega, eps_m, b, one_term, &
                         sre, sim)
    implicit double precision (a-h,o-z)
    dimension xnode(n), wnode(n)
    logical :: one_term
    sre = 0.0d0
    sim = 0.0d0
    do k = 1, n
      g = 1.0d0/(1.0d0+(xnode(k)/b)**2)
      dim  = omega + xnode(k)
      dre  = -eps_m
      dmag2 = dre*dre + dim*dim
      tre =  dre/dmag2
      tim = -dim/dmag2
      if (.not. one_term) then
        dim2 = omega - xnode(k)
        dmag22 = dre*dre + dim2*dim2
        tre = tre + dre/dmag22
        tim = tim - dim2/dmag22
      end if
      factor = -wnode(k)/(2.0d0*pi_local())
      sre = sre + factor*tre*g
      sim = sim + factor*tim*g
    end do
  end subroutine folded_sum

  function pi_local() result(p)
    double precision :: p
    p = 4.0d0*datan(1.0d0)
  end function pi_local

  ! Independent ground truth via w'=b*tan(theta) substitution +
  ! fine composite Simpson rule over theta in (-pi/2+d, pi/2-d).
  subroutine exact_ref(omega, eps_m, b, sre, sim)
    implicit double precision (a-h,o-z)
    integer, parameter :: nth = 400000  ! even, for Simpson
    double precision :: pi, delta, th, h, wre, wim
    double precision :: fre, fim, wp, denre, denim, denmag2
    integer :: i

    pi = 4.0d0*datan(1.0d0)
    delta = 1.0d-9
    th_lo = -pi/2.0d0 + delta
    th_hi =  pi/2.0d0 - delta
    h = (th_hi-th_lo)/dble(nth)

    sre = 0.0d0
    sim = 0.0d0
    do i = 0, nth
      th = th_lo + h*dble(i)
      wp = b*dtan(th)
      ! integrand (before the outer -1/2pi and the b*dtheta jacobian):
      !   1/(i*(omega+wp)-eps_m)
      denre = -eps_m
      denim = omega+wp
      denmag2 = denre*denre+denim*denim
      fre =  denre/denmag2
      fim = -denim/denmag2
      ! Simpson weight
      if (i == 0 .or. i == nth) then
        wgt = 1.0d0
      else if (mod(i,2) == 1) then
        wgt = 4.0d0
      else
        wgt = 2.0d0
      end if
      sre = sre + wgt*fre
      sim = sim + wgt*fim
    end do
    sre = sre*h/3.0d0
    sim = sim*h/3.0d0
    ! jacobian b*dtheta already folded into h; multiply by b, then by
    ! the outer -1/(2pi) prefactor of Sigma_c
    sre = -b*sre/(2.0d0*pi)
    sim = -b*sim/(2.0d0*pi)
  end subroutine exact_ref

end program pilot_sigma_fold_check

! ---------------------------------------------------------------
! Verbatim copies of the production grid routines (same duplicated-
! per-module convention as scrnc_sigma_driver.F itself).
subroutine scrnc_sigma_freq_grid(n, omega_c, xnode, wnode, maxfreq)
  implicit double precision (a-h,o-z)
  dimension xnode(maxfreq), wnode(maxfreq)
  dimension xleg(maxfreq), wleg(maxfreq)

  call scrnc_sigma_gaussl(-1.0d0, 1.0d0, n, xleg, wleg, maxfreq)

  do k = 1, n
    xnode(k) = omega_c*(1.0d0+xleg(k))/(1.0d0-xleg(k))
    wnode(k) = wleg(k)*omega_c*2.0d0/ &
               ((1.0d0-xleg(k))*(1.0d0-xleg(k)))
  end do
end subroutine scrnc_sigma_freq_grid

subroutine scrnc_sigma_gaussl(x1, x2, n, x, w, maxfreq)
  implicit double precision (a-h,o-z)
  dimension x(maxfreq), w(maxfreq)
  parameter (eps=3.d-14)
  m = (n+1)/2
  xm = 0.5d0*(x2+x1)
  xl = 0.5d0*(x2-x1)
  do i = 1, m
    z = cos(3.141592654d0*(i-0.25d0)/(n+0.5d0))
    do
      p1 = 1.d0
      p2 = 0.d0
      do j = 1, n
        p3 = p2
        p2 = p1
        p1 = ((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
      end do
      pp = n*(z*p1-p2)/(z*z-1.d0)
      z1 = z
      z = z1-p1/pp
      if (dabs(z-z1) <= eps) exit
    end do
    x(i)     = xm-xl*z
    x(n+1-i) = xm+xl*z
    w(i)     = 2.d0*xl/((1.d0-z*z)*pp*pp)
    w(n+1-i) = w(i)
  end do
end subroutine scrnc_sigma_gaussl
