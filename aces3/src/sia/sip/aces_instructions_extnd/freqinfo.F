c---------------------------------------------------------------|
c  This subroutine is a superinstruction for the SIAL code      |
c  It provides value of frequencies for dynamic polarizability  |
c             calculations.                                     |
c  Added by Denis Bokhan, Moscow Lomonosov State University,    |
c              Moscow, RUSSIAN FEDERATION                       |
c---------------------------------------------------------------|
      subroutine freqinfo(x,nindex, type, bval,eval, bdim, edim,
     *           x2, nindex2, type2, bval2, eval2, bdim2, edim2)


      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)
      double precision x2(*)
      integer nindex2, type2(*), bval2(*), eval2(*)
      integer bdim2(*), edim2(*)
c----------------------------------------------------------------

      call return_freq(x,x2 )
c----------------------------------------------------------------
      return
      end

C----------------------------------------------------------------
C----------------------------------------------------------------
      subroutine return_freq(x1,x2)
      implicit none
      include 'parallel_info.h'
      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
C      include 'freq.h'
C     freq.h's /FR/ common (FREQ, IMFREQ) can't be included directly here --
C     it would clash with int_gen_parms.h's own freq/imfreq names (the
C     ASV-parsed values read below). Declare local aliases bound to the same
C     /FR/ common block (matched by block name, not variable name) so the
C     real ASV frequency actually reaches eps_divide*_fr / eps_divide*_imfr
C     in energy_denominator_freq.F / energy_denominator_imfreq.F, which read
C     FREQ/IMFREQ from /FR/ -- previously left permanently at 0.0 since
C     nothing ever wrote to /FR/, making the +FREQ/-FREQ (and IMFREQ) branches
C     in those routines identical.
      double precision freq_fr, imfreq_fr
      common /FR/ freq_fr, imfreq_fr

      integer a1,a2,b1,b2,c1,c2,d1,d2,i1,i2,j1,j2,k1,k2
      integer i,j,k,l,a,b,l1,l2,g,kk1,ll1
      integer type1,type2,type3,type4,iatom
      double precision x1,x2

      x1 = FREQ
      x2 = IMFREQ
      freq_fr = FREQ
      imfreq_fr = IMFREQ
c--------------------------------------------------------------------------------
      return
      end

