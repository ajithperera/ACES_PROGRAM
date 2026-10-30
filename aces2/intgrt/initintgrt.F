      subroutine initintgrt(jam,intnumradpts,islaterON)

c This program initializes parameters from the INTGRT namelist.

      implicit none
#include <aces.h>
#include <intgrtflags.com>
#include <dftfunc.com>
      integer err,ptr,num,i,cutoff
      logical jam
C  prakash 23/04
C  scfksexact for exact ks calculation
C  end 23/04
      
       integer islaterON

       integer intnumradpts
      character*20 list(100),ele
      logical printdef, streq
      character*10 abbx
      character*50 namx
      M_REAL
     &    coef,thresh_power,rsomega,wh,wfunc,wpot,
     &    cxh(fun_num_exch),cch(fun_num_corr)

        printdef=.true.


C prakash
C  int_ks and int_ks_finaliter doesn't make any sence for time
C  so need to be commented
c      if (int_ks .and. .not.int_ks_finaliter) printdef=.false.
C end prakash



      call callstack_push('INITINTGRT')
      call nl_init('INTGRT',err,printdef)

C prakash
C  nl_init try to find *INTGRT in ZMAT if it is there it will read
C required information if not found then default value would be
C passed
C


c Namelist entries:
c --------------------------------------------------------------------------
c PRINT_NL   : whether to print the values read in in the namelist
c              Default: true
c --------------------------------------------------------------------------
c Print flags
c
c              All print flags have one of the following values:
c                never   - never print out the property
c                dft     - only print out the property during a DFT calculation
c                ks      - print it out durint the final iteration of KS
c                          or during a DFT calculation
c                always  - print it out during all KS iterations as well
c                          as during a DFT calculation
c
c PRINT_SCF  : whether or not to print a summary of SCF (or correlated in
c              the case of a relaxed density) results
c              Default: dft
c PRINT_INT  : whether to print out a summary of integration parameters
c              Default: dft
c PRINT_SIZE : whether to print relative size information
c              Default: never
c PRINT_ATOM : whether to print out integration information for each atom
c              Default: never
c
c PRINT_MOS  : whether to print the MOs from the SCF or correlated calculation
c              Default: never
c PRINT_OCC  : whether to print the occupancy of each MO
c              Default: never
c --------------------------------------------------------------------------
c Integration flags
c
c POTRADPTS  : the number of radial points to use in the integration around
c              each atomic center in KS calc (default: 50)
c RADTYP     : the type of radial integration available are:
c                Handy           : (default)
c                Gauss-Legendre  :
c PARTPOLY   : partition polyhedra
c                equal     : equally regardless of atomic center
c                bsrad     : by Bragg-Slater radii (default)
c                dynamic   : dynamic (minimum in density) sizing
c RADSCAL    : how to do radial scaling
c                none      : no radial scaling
c                slater    : use Slater's rules for radial scaling (default)
c PARTTYP    : type of partitioning
c                rigid     : rigid partitioning
c                fuzzy     : fuzzy partitioning (default)
c FUZZYITER  : the number of iterations to use in fuzzy partitioning
c              (default: 4)
c RADLIMIT   : the length of the radial integration (default 3.00)
c
c ENEGRID    : the angular grid to use by default at each radial point
c              in an energy calculation (default: 4)
c POTGRID    : the angular grid to use by default at each radial point
c              in a Kohn-Sham calculation of the potential (default: 4)
c ENETYPE    : the angular grid type to use by default at each radial point
c              in an energy calculation
c                lebedev    : default
c POTTYPE    : the angular grid type to use by default at each radial point
c              in a Kohn-Sham calculation of the potential
c                lebedev    : default
c --------------------------------------------------------------------------
c Functionals
c
c FUNC       : The type of functional (and potential) to use in any
c              calculation other than simple HFDFT.  In the case of HFDFT,
c              all energy functionals are calculated.  In the case of
c              Kohn-Sham, only one potential will be calculated.  In the
c              case of an energy search, only one energy is stored in the
c              JOBARC file.
c
c              This keyword is of one of the forms:
c                FUNC=func
c                FUNC=exch,corr
c                FUNC=func,coef,func,coef [,func,coef]*
c              where func is the name of one of the functionals below and
c              coef is the coefficient to use.  In the final form, at least
c              two functionals must be included in the list.
c
c              Recognized correlation functionals are:
c                vwn     : Vosko-Wilk-Nusair, used as VWN5 in Gaussian94 
c                lyp     : Lee-Yang-Parr
c              pbe_cor   : Perdew-Burke-Ernzerhof
c              pw91_cor  : Perdew-Wang 91
c               vwn_rpa  : Vosko-Wilk-Nusair fitted to RPA, with the
c                          plain f(zeta) spin interpolation (VWN_RPA
c                          of libxc, used by qtp17)
c
c                wl      : Wilson-Levy
c                wi      : Wilson-Ivanov
c                wl and wi do not have their corresponding potentials         
c                wi has form for open-shell species only
c
c              Exchange functionals are:
c                lda     : LDA (Slater) exchange
c                becke   : Becke exchange
c               pbe_ex   : Perdew-Burke-Ernzerhof
c               pw91_ex  : Perdew-Wang 91
c               becke_sr : short-range Becke exchange, the Becke
c                          functional attenuated as in Iikura, Tsuneda,
c                          Yanai and Hirao, JCP 115, 3540 (2001)
c                hf      : Hartree-Fock-like exchange
c                hf_lr   : long-range Hartree-Fock-like exchange, with
c                          erf(rs_omega*r12)/r12 integrals
c
c              Hybrid functionals (their coefficients and references
c              are in hybparm.F) are:
c                b3lyp   : Becke III LYP
c              cam-b3lyp : 0.19 hf + 0.46 hf_lr + 0.35 becke
c                          + 0.46 becke_sr, 0.19 vwn + 0.81 lyp,
c                          rs_omega=0.33
c              cam-qtp-00: 0.54 hf + 0.37 hf_lr + 0.09 becke
c                          + 0.37 becke_sr, 0.20 vwn + 0.80 lyp,
c                          rs_omega=0.29
c              cam-qtp-01: 0.23 hf + 0.77 hf_lr + 0.77 becke_sr,
c                          0.20 vwn + 0.80 lyp, rs_omega=0.31
c              cam-qtp-02: 0.28 hf + 0.72 hf_lr + 0.72 becke_sr,
c                          lyp, rs_omega=0.335
c                lc-qtp  : hf_lr + becke_sr, lyp, rs_omega=0.475
c                qtp17   : 0.62 hf + 0.38 lda, 0.20 vwn_rpa + 0.80 lyp
c
c              A hybrid functional may only be entered as FUNC=hybrid.
c
c              The default functional is none.
c
c RS_OMEGA   : the range-separation parameter omega (bohr**-1) used by
c              becke_sr and hf_lr. Default: the omega of the
c              range-separated hybrid named in FUNC or KSPOT, else 0
c              (becke_sr and hf_lr then need it given). VMOL reads it
c              too (see vmol/rsomega.F) to write the long-range
c              integral files IIIILR, ...
c========================================================================
c                         More Flags                                    =
c========================================================================
c
c     numacc         Numerical accelerator. Default: on
c
c     cutoff         Cutoff value for the numerical integration.
c                    Default: 12 Densities < 10**(-12) are assumed 0.
c     exact_ex       Exact LOCAL EXCHANGE. Default: false
c     tdks           Time-dependnet KS is requested. Default: false 

      list(1)='never'
      list(2)='dft'
      list(3)='ks'
      list(4)='always'
      call nl_ele('print_scf' ,list,4,int_prt_dft  ,int_printscf)
      call nl_ele('print_int' ,list,4,int_prt_dft  ,int_printint)
      call nl_ele('print_atom',list,4,int_prt_never,int_printatom)
      call nl_ele('print_size',list,4,int_prt_never,int_printsize)
      call nl_ele('print_mos' ,list,4,int_prt_never,int_printmos)
      call nl_ele('print_occ' ,list,4,int_prt_never,int_printocc)

      call nl_int('potradpts',50,int_numradpts)
      call nl_int('OEPSLATER',5,islaterON)

c  prakash june 16 about intpack for oep
c we need int_overlp info
      call nl_int('overind',3, int_overlp)
c end prakash
      intnumradpts=int_numradpts
C
      list(1)='Handy'
      list(2)='Gauss-Legendre'
      list(3)='Becke'
      list(4)='TA'
      list(5)='multiEX'
      call nl_ele('radtyp',list,5,int_radtyp_handy,int_radtyp)

      list(1)='equal'
      list(2)='bsrad'
      list(3)='dynamic'
      call nl_ele('partpoly',list,3,int_partpoly_bsrad,int_partpoly)

      list(1)='none'
      list(2)='Slater'
      list(3)='Bragg'
      list(4)='TA'
      list(5)='multiEX'
      call nl_ele('radscal',list,5,int_radscal_slater,int_radscal)

      list(1)='rigid'
      list(2)='fuzzy'
      call nl_ele('parttyp',list,2,int_parttyp_fuzzy,int_parttyp)

      call nl_int('fuzzyiter',4,int_fuzzyiter)

      call nl_real('radlimit',3.d0,int_radlimit)
      call nl_log('numacc',.true.,num_acc_ks)
      call nl_log('exact_ex',.false.,ks_exact_ex)


      call nl_log('tdks',.false.,int_tdks)

      call nl_int('enegrid',4,int_defenegrid)
      call nl_int('potgrid',4,int_defpotgrid)
      list(1)='lebedev'
      call nl_ele('enetype',list,1,int_gridtype_leb,int_defenetype)
      call nl_ele('pottype',list,1,int_gridtype_leb,int_defpottype)


C prakash
C KS_maxmem input
C if user provide the information then it will change otherwise it will
C remain as default i.e. 1500000

       call nl_int('ksmem',1500000,int_ksmem)

c end prakash

c Parse the functional to use:

      do i=1,fun_num_exch
        coef_exch(i)=0.d0
      end do
      do i=1,fun_num_corr
        coef_corr(i)=0.d0
      end do

      nam_exch(fun_exch_lda)  = 'LDA (Slater, Xalpha)'
      nam_exch(fun_exch_becke)= 'Becke'
      nam_exch(fun_exch_pbe)= 'Perdew-Burke-Ernzerhof'
      nam_exch(fun_exch_pw91)= 'Perdew-Wang 91'
      nam_exch(fun_exch_becke_sr)= 'Becke short-range'
      nam_exch(fun_exch_hflr) = 'Exact long-range'
      nam_exch(fun_exch_hf)   = 'Exact Nonlocal Exchange'

      nam_corr(fun_corr_vwn)  = 'Vosko-Wilk-Nusair'
      nam_corr(fun_corr_lyp)  = 'Lee-Yang-Parr'
      nam_corr(fun_corr_pbe)  = 'Perdew-Burke-Ernzerhof'
      nam_corr(fun_corr_pw91) = 'Perdew-Wang 91'
      nam_corr(fun_corr_wl)   = 'Wilson-Levy'
      nam_corr(fun_corr_wi)   = 'Wilson-Ivanov'
      nam_corr(fun_corr_vwnrpa) = 'VWN RPA'

      abb_exch(fun_exch_lda)  = 'lda'
      abb_exch(fun_exch_becke)= 'becke'
      abb_exch(fun_exch_pbe)  = 'pbe_ex'
      abb_exch(fun_exch_pw91) = 'pw91_ex'
      abb_exch(fun_exch_becke_sr)= 'becke_sr'
      abb_exch(fun_exch_hflr) = 'hf_lr'
      abb_exch(fun_exch_hf)   = 'hf'

      abb_corr(fun_corr_vwn)  = 'vwn'
      abb_corr(fun_corr_lyp)  = 'lyp'
      abb_corr(fun_corr_pbe)  = 'pbe_cor'
      abb_corr(fun_corr_pw91) = 'pw91_cor'
      abb_corr(fun_corr_wl)   = 'wl'
      abb_corr(fun_corr_wi)   = 'wi'
      abb_corr(fun_corr_vwnrpa) = 'vwn_rpa'

      do i=1,fun_num_hyb
        call hybparm(i,abb_hyb(i),nam_hyb(i),cxh,cch,wh)
      end do
      wfunc=0.d0
      wpot=0.d0

      call nl_str('func','none',nam_func)

      ptr=0
      call list_next_ele(nam_func,ptr,ele,err)

      if (err.ne.0) then
        write (*,910)
        call errex
        stop
      endif

      if (streq(ele,'none',.true.)) then
         int_dft_fun=fun_dft_none
         go to 100
      else
         int_dft_fun=fun_special
      end if  

      call list_next_num(nam_func,ptr,coef,err)

  910 format ('ERROR: invalid functional.')
  920 format ('ERROR: functional ',A10,'name used twice.')
  921 format ('ERROR: functional name unrecognized.')
  922 format ('ERROR: invalid functional coefficient.')
  923 format ('ERROR: 2 or more functionals required.')

c If we're working with the form: FUNC=func,coef,func,coef [,func,coef]*

      if (err.eq.0) then
        num=0
   10   continue

c       At this point, ele contains the name of the functional and coef
c       contains it's coefficient.

c       Check exchange functionals...

        do i=1,fun_num_exch
          if (streq(ele,abb_exch(i),.true.)) then
            if (coef_exch(i).ne.0.d0) then
              write (*,920) abb_exch(i)
              call errex
              stop
            endif
            coef_exch(i)=coef
            num=num+1
            goto 20
          endif
        end do

c       Check correlation functional

        do i=1,fun_num_corr
          if (streq(ele,abb_corr(i),.true.)) then
            if (coef_corr(i).ne.0.d0) then
              write (*,920) abb_corr(i)
              call errex
              stop
            endif
            coef_corr(i)=coef
            num=num+1
            goto 20
          endif
        end do

c       If neither, die.

        write (*,921)
        call errex
        stop

   20   continue

c       Get the next functional and coefficient.  If they don't exist,
c       make sure there have been at least 2.

        call list_next_ele(nam_func,ptr,ele,err)

        if (err.eq.0) then
          call list_next_num(nam_func,ptr,coef,err)
          if (err.eq.0) goto 10
          write (*,922)
          call errex
          stop
        else if (num.lt.2) then
          write (*,923)
          call errex
          stop
        endif

c We're working with the form FUNC=func or FUNC=corr,exch or FUNC=hybrid

      else
        num=0
c       num=0  : nothing entered
c           1  : exchange entered
c           2  : correlation entered

        do i=1,fun_num_exch
          if (streq(ele,abb_exch(i),.true.)) then
            coef_exch(i)=1.d0
            num=num+1
            goto 60
          endif
        end do

        do i=1,fun_num_corr
          if (streq(ele,abb_corr(i),.true.)) then
            coef_corr(i)=1.d0
            num=num+2
            goto 60
          endif
        end do

   60   continue
        if (num.eq.0) go to 70
         call list_next_ele(nam_func,ptr,ele,err)
          if (err.ne.0.d0) go to 100

        do i=1,fun_num_exch
          if (streq(ele,abb_exch(i),.true.)) then
             if (num.eq.2) then
               coef_exch(i)=1.d0
               goto 100
             else
               write (*,910)
               call errex
               stop
             end if 
          end if
        end do


        do i=1,fun_num_corr
          if (streq(ele,abb_corr(i),.true.)) then
             if (num.eq.1) then
               coef_corr(i)=1.d0
               goto 100
             else
               write (*,910)
               call errex
               stop
             end if
          end if
        end do

           write (*,910)
           call errex
           stop

   70   continue

c       Hybrid functionals...

        do i=1,fun_num_hyb
          if (streq(ele,abb_hyb(i),.true.)) then
             int_dft_fun=i
             call hybparm(i,abbx,namx,coef_exch,coef_corr,wfunc)
             goto 100
          end if
        end do
        write (*,910)
        call errex
        stop
       endif
  100 continue



c KSPOT      : Availiable potentials : 
c                vwn     : Vosko-Wilk-Nusair 
c                lyp     : Lee-Yang-Parr
c
c                lda     : LDA (Slater) exchange 
c                becke   : Becke exchange
c                hf      : Nonlocal HF-like (default)
c
c          The format and the sructure are the same as those 
c          for the key word "func".

c     Initialize the coefficients needed in hybrid schemes.
c     Note that hybrid scheme is a more general definition and the KS
c     method with multiplicative potentials is one of the possible choices.

      do i=1,fun_num_exch
        coef_pot_exch(i)=0.d0
      end do
      do i=1,fun_num_corr
        coef_pot_corr(i)=0.d0
      end do
 
      call nl_str('kspot','hf',nam_kspot)
c prakash
C jam will control weather we have to num integration or not

      jam=streq(nam_kspot,'hf',.true.)
C
      ptr=0
      call list_next_ele(nam_kspot,ptr,ele,err)
      if (err.ne.0) then
        write (*,1010)
        call errex
        stop
      endif
      int_kspot=fun_special
      call list_next_num(nam_kspot,ptr,coef,err)

 1010 format ('ERROR: invalid potential.')
 1020 format ('ERROR: potential ', A7, ' name used twice.')
 1021 format ('ERROR: potential name unrecognized.')
 1022 format ('ERROR: invalid potential coefficient.')
 1023 format ('ERROR: 2 or more potentials required.')
 1024 format ('ERROR: potential not implemented.')
c If we're working with the form: kspot=pot,coef,pot,coef [,pot,coef]*
c We could have any combination of exch and corr. 
c Each pot appears only once.

      if (err.eq.0) then
        num=0
  110   continue

c       At this point, ele contains the name of the potential and coef
c       contains it's coefficient.

c       Check exchange potentials...

        do i=1,fun_num_exch
          if (streq(ele,abb_exch(i),.true.)) then
            if (coef_pot_exch(i).ne.0.d0) then
              write (*,1020) abb_exch(i)
              call errex
              stop
            endif
            coef_pot_exch(i)=coef
            num=num+1
            goto 120
          endif
        end do

c       Check correlation potentials...

       do i=1,fun_num_corr
          if (streq(ele,abb_corr(i),.true.)) then
            if (coef_pot_corr(i).ne.0.d0) then
              write (*,1020) abb_corr(i)
              call errex
              stop
            endif
            coef_pot_corr(i)=coef
            num=num+1
            goto 120
          endif
        end do

c       If neither, die.

        write (*,1021)
        call errex
        stop

  120   continue

c       Get the next potential and coefficient.  If they don't exist,
c       make sure there have been at least 2.

        call list_next_ele(nam_kspot,ptr,ele,err)

        if (err.eq.0) then
          call list_next_num(nam_kspot,ptr,coef,err)
          if (err.eq.0) goto 110
          write (*,1022)
          call errex
          stop
        else if (num.lt.2) then
          write (*,1023)
          call errex
          stop
        endif

c We're working with the form kspot=pot, (exch or corr),
c or kspot=exch,corr (or corr,exch), or kspot=standard hybrid

      else
        num=0

c       num=0  : nothing entered
c           1  : exchange entered
c           2  : correlation entered

        do i=1,fun_num_exch
          if (streq(ele,abb_exch(i),.true.)) then
            coef_pot_exch(i)=1.d0
            num=num+1
            goto 160
          endif
        end do

        do i=1,fun_num_corr
          if (streq(ele,abb_corr(i),.true.)) then
            coef_pot_corr(i)=1.d0
            num=num+2
            goto 160
          endif
        end do

  160   continue
        if (num.eq.0) go to 170
         call list_next_ele(nam_kspot,ptr,ele,err)
         if (err.ne.0) go to 200

        do i=1,fun_num_exch
          if (streq(ele,abb_exch(i),.true.)) then
             if (num.eq.2) then
               coef_pot_exch(i)=1.d0
               goto 200
             else
               write (*,1010)
               call errex
               stop
             end if
          end if
        end do


        do i=1,fun_num_corr
          if (streq(ele,abb_corr(i),.true.)) then
             if (num.eq.1) then
               coef_pot_corr(i)=1.d0
               goto 200
             else
               write (*,1010)
               call errex
               stop
             end if
          end if
        end do

        write (*,1010)
        call errex
        stop 

  170   continue

c       Hybrid functionals.  B3LYP has a potential path of its own,
c       the others are built term by term (fun_special).

        do i=1,fun_num_hyb
          if (streq(ele,abb_hyb(i),.true.)) then
             int_kspot=fun_special
             if (i.eq.fun_hyb_b3lyp) int_kspot=fun_hyb_b3lyp
             call hybparm(i,abbx,namx,coef_pot_exch,coef_pot_corr,wpot)
             goto 200
          end if
        end do
        write (*,1010)
        call errex
        stop
       endif

  200 continue

c Save the value of coef_pot_exch(fun_exch_hf) also as
c coef_pot_nonlocal. Needed for Vscf.

        if (ks_exact_ex) then 
           coef_pot_nonlocal=1.d0
           coef_pot_lrx=0.d0
        else
           coef_pot_nonlocal=coef_pot_exch(fun_exch_hf)
           coef_pot_lrx=coef_pot_exch(fun_exch_hflr)
        end if

c Range separation: omega defaults to the one of the hybrid named in
c func or kspot (vmol/rsomega.F applies the same rule).

      if (wfunc.ne.0.d0 .and. wpot.ne.0.d0 .and. wfunc.ne.wpot) then
         write (*,1032)
         call errex
         stop
      end if
      rsomega=wfunc
      if (rsomega.eq.0.d0) rsomega=wpot
      call nl_real('rs_omega',rsomega,int_rs_omega)
      if (int_rs_omega.lt.0.d0) then
         write (*,1030)
         call errex
         stop
      end if
      if (int_rs_omega.eq.0.d0 .and.
     &    (coef_exch(fun_exch_hflr).ne.0.d0 .or.
     &     coef_exch(fun_exch_becke_sr).ne.0.d0 .or.
     &     coef_pot_exch(fun_exch_hflr).ne.0.d0 .or.
     &     coef_pot_exch(fun_exch_becke_sr).ne.0.d0)) then
         write (*,1031)
         call errex
         stop
      end if
 1030 format ('ERROR: rs_omega must not be negative.')
 1031 format ('ERROR: becke_sr and hf_lr need rs_omega.')
 1032 format ('ERROR: func and kspot name hybrids with different ',
     &        'rs_omega.')

c If kspot=wl or wi 
        if ((coef_pot_corr(fun_corr_wl).ne.0.d0).or.
     &    (coef_pot_corr(fun_corr_wi).ne.0.d0)) then
          write (*,1024)
          call errex
          stop
        endif

      call nl_int('cutoff',12,cutoff)
      thresh_power=dble(cutoff)
      thresh_power=dabs(thresh_power)

       if (thresh_power.lt.8) then
         write (*,5010)
      end if

 5010 format(' !!! WARNING !!! LOW NUMERICAL ACCURACY! ')

      thresh=1.0d0/(10.0d0**thresh_power) 


      call nl_term
      call callstack_pop
      return
      end
