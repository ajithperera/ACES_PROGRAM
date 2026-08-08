










            subroutine nocc2occup

C   this routine will occupy nocc to occup for dft works








c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler























cYAU - ACES3 stuff . . . we hope - #include <aces.par>




         implicit none


      integer           iuhf
      common /iuhf_com/ iuhf
      save   /iuhf_com/




      integer maxirrep,num2comb,max2comb
      parameter (maxirrep=8)
      parameter (num2comb=22)
      parameter (max2comb=25)

c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end
      integer        nirrep, numbasir(8),
     &               irpsz1(36),irpsz2(28),irpds1(36),irpds2(56),
     &               old_irpoff(9), irrorboff(9), dirprd(8,8),
     &               old_iwoff1(37), old_iwoff2(29),
     &               inewvc(maxbasfn), idxvec(maxbasfn),
     &               irrtrilen(9), irrtrioff(8),
     &               irrsqrlen(9), irrsqroff(8)
      common /symm2/ nirrep, numbasir,
     &               irpsz1,    irpsz2,    irpds1,    irpds2,
     &               old_irpoff,    irrorboff,    dirprd,
     &               old_iwoff1,     old_iwoff2,
     &               inewvc,           idxvec,
     &               irrtrilen,    irrtrioff,
     &               irrsqrlen,    irrsqroff
      save   /symm2/

      integer             occup(8,2),totocc(2),totocca,totoccb,
     &                    maxirrtri,maxirrsqr,irrtritot,irrsqrtot
      common /sym_ks_com/ occup,     totocc,   totocca,totoccb,
     &                    maxirrtri,maxirrsqr,irrtritot,irrsqrtot
      save   /sym_ks_com/


         integer i,jj
         integer NOCC(16)

         common /POPUL/NOCC



         do i=1,8
             occup(i,1)=NOCC(i)
         end do
 
          do jj=1,8
             occup(jj,2)=NOCC(8+jj)
          end do
         
         return
         end




