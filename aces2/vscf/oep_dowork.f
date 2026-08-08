










      subroutine oep_dowork(iter, iuhf, nbas, nocc, icore, maxmem,  
     &                      usdmem, scfh1e, scfevl, scfevc, scfden, 
     &                      scfham, scfenr, dmax, itertype, 
     &                      transfer,slat,Zslater,islaterON)
c
      implicit none
c
c
c This common block defines the OEP variables of general scope, mostly
c flags and their default values.
c
c Igor Schweigert, Feb 2004
c $Id: oep.com,v 1.8 2008/06/06 18:09:17 taube Exp $
c 
c
      logical
     &     OEP1, OEP2, OEP1ON, OEP2ON, OEP2ALWAYS, OEP2ITR, 
     &     OEPDFF, OEPDFW, OEPDWW, OEPFDF, OEPWDW, BOEPR12,
     &     R12ON
c     
      integer
     &     OEP1ONTHR, OEP1ONTHRDEF, OEP2ONTHR, OEP2ONTHRDEF,
     &     OEP2OFFTHR, OEP2OFFTHRDEF, OEP2UPDATE, OEP2UPDATEDEF,
     &     OEP2ALWAYSTHR, OEP2ALWAYSTHRDEF, NAUX,
     &     OEPSVDTHR, OEPSVDTHRDEF, OEP1HYBRIDINT, OEP1HYBRIDDEF,
     &     R12ONTHR, R12ONTHRDEF
c
      double precision OEP1HYBRID
c     
      character*(50)
     &     OEPH0, OEPH0DEF, OEPAXP, OEPAXPDEF
     
c      character*1
c     &     asv_null  
c     
      common /OEPCONFIG/
     &     OEP1, OEP2, OEP1ON, OEP2ON, OEP2ALWAYS, OEP2ITR,
     &     BOEPR12, R12ON, OEPDFF, OEPDFW, OEPDWW, OEPFDF, OEPWDW,
     &     OEP1ONTHR, OEP2ONTHR, OEP2OFFTHR, OEP2ALWAYSTHR, 
     &     OEPSVDTHR, NAUX, OEP2UPDATE, R12ONTHR,OEP1HYBRID, 
     &     OEPH0, OEPAXP
c     
      parameter (OEP1ONTHRDEF  = 1)
      parameter (OEP2ONTHRDEF  = 2)
      parameter (OEP2OFFTHRDEF = 20)
      parameter (OEP2ALWAYSTHRDEF = 6)
      parameter (OEP2UPDATEDEF = 1)
      parameter (OEPSVDTHRDEF       = 8)
      parameter (OEP1HYBRIDDEF = 0)
      parameter (R12ONTHRDEF  = 2)
c
      data
     &     OEPH0DEF, OEPAXPDEF   /'OEP', 'NONE'/
c     
      


c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end



c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
c symm2.com : begin

c This is initialized in vscf/symsiz.

c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end
      integer nirrep,      nbfirr(8),   irpsz1(36),  irpsz2(28),
     &        irpds1(36),  irpds2(56),  irpoff(9),   ireps(9),
     &        dirprd(8,8), iwoff1(37),  iwoff2(29),
     &        inewvc(maxbasfn),         idxvec(maxbasfn),
     &        itriln(9),   itriof(8),   isqrln(9),   isqrof(8),
     &        mxirr2
      common /SYMM2/ nirrep, nbfirr, irpsz1, irpsz2, irpds1, irpds2,
     &               irpoff, ireps,  dirprd, iwoff1, iwoff2, inewvc,
     &               idxvec, itriln, itriof, isqrln, isqrof, mxirr2
c symm2.com : end
c 
      character*3 itertype
c
c     The input parameters
c
      logical addoep2, transfer
c 
      integer
     &     nbas, nocc(16), iuhf,
     &     icore (*) 

      logical slat,Zslater
      integer islaterON
      double precision
     &     scfh1e (itriln (nirrep+1)),
     &     scfevl (nbas), scfevc (nbas, nbas),
     &     scfden (itriln (nirrep+1)),
     &     scfham (itriln (nirrep+1)), scfenr
c
c     Local variables
c
      integer
     &     i000, i010, coretop, usdmem, maxmem, iter, luout,
     &     auxaa, auxph, hfxaa, hfxph, 
     &     resinv, oep1au, oep2au, refevl
c
      integer lnbuf, luint
      parameter (lnbuf=600)
      parameter (luint=10, luout=6)
c
      double precision homo, OEPTOL, au2ev, dmax, tol, etot
c
c     Define the one- and two-electron integrals array sizes.
c     
C#include <oepintind.com>
C#include <oepintind.com>
c     
c     OEP SCF block, Igor Schweigert, Sep 2003,
c     $Id: vscfoep.FPP,v 1.13 2008/06/06 18:10:43 taube Exp $
c     
c     A SCF iteration starts with constructing the Fock matrix. If the
c     OEP method is used, we will solve the OEP equation, build the OEP
c     Hamiltonian matrix and bypass the standard Fock matrix
c     construction and the total energy calculation.
c     
c     Note that we need some initial guess for eigenvalues and
c     eigenvectors. The INITGES routine does not store the 
c     eigenvectors. So if INITGES is used let's perform the HF SCF
c     iteration once as a guess for eigenvectors.
c     
c     Note also that we bypass the "reaction field" contributions, which
c     is actually the first thing done in SCF iterations. The reason for
c     that is if the OEP block is the first in a SCF iteration then we
c     still have eigenvalues and eigenvectors in ICORE (I050) and ICORE
c     (I060) correspondingly after SCFIT. So we don't need to create
c     extra arrays to store them.
c     
c     Note also that all the incompatabilities have already been checked
c     in the OEP preSCF block. So the fact that we get here means that
c     this is RHF calculation (IUHF=0), the VMOL integral package is
c     used (PETITE_LIST=.FALSE.), and the Fock matrix is build out of
c     core (AOFIL=.TRUE.)
c     
c     The OEP SCF block consists of the following steps: <> Figure out
c     what kind of iteration is the current one (all the games with the
c     convergence should be played here.) <> If it is a OEP1 or OEP2
c     iteration: <<>> Set the top of the core behind the SCF
c     eigenvectors. <<>> Solve the corresponding OEP equations and build
c     the OEP Hamiltonian. <<>> Report the iteration type, total energy,
c     and density difference from the previous iteration. <<>> Remember
c     the type of the iteration just performed. <<>> Proceed with the
c     rest of the SCF iteration by jumping to DMPSCF.
C
      r12on=.false.
      call oeptimer ('VSCF: SCFITER', 'ON')

      if(zslater) then 
        if( dmax .le. 10.d0**(-islaterON) .and. iter .gt. 1) then
          if( slat ) then
             call slater
             slat=.false.
          end if 
       end if
      end if
C
      if( .not. slat) then
      if (OEP1.or.OEP2.and..not.OEP1ON) then
         if (dmax.le.10.d0**(-OEP1ONTHR).and.iter.gt.1)
     &        OEP1ON = .TRUE.
      end if
       if (OEP2.and..not.OEP2ON) then
         if (dmax.le.10.d0**(-OEP2ONTHR).and.iter.gt.1) then
            OEP2ON = .TRUE.
            OEP2ITR = .TRUE.
         endif 
      elseif (OEP2.and.OEP2ON) then
         OEP2ITR = .FALSE.
         if (dmax.le.10.d0**(-OEP2ALWAYSTHR))
     &        OEP2ALWAYS = .TRUE.
         if (dmax.le.10.d0**(-OEP2OFFTHR).and..not.OEP2ALWAYS)
     &        OEP2UPDATE = 9999
         if (OEP2ALWAYS.or.iter - (iter/OEP2UPDATE)*OEP2UPDATE .eq. 0)
     &        OEP2ITR = .TRUE.
         if (dmax.le.0.0000) then
            tol = OEPTOL
            OEP2ITR = .TRUE. 
         endif
      endif
c
      i000 = 1

      if (OEP1ON .or. OEP2ON) then
           
         i010 = i000 + (iuhf+1) * isqrln (nirrep+1) * iintfp
         if (.not.OEP2ITR) then
            call oep1iter (
     &           nbas, nocc (1),  
     &           icore (i000), maxmem, usdmem,  
     &           scfh1e, scfevl, scfevc, scfden, scfham, scfenr)
         else
            call oep2iter (
     &           nbas, nocc (1), 
     &           icore (i000), maxmem, usdmem, 
     &           scfh1e, scfevl, scfevc, scfden, scfham, scfenr,
     &           r12on)
         endif
         etot = scfenr
         if (.not.OEP2ON) then
            itertype = 'x  '
         elseif (OEP2ITR) then
            itertype = 'xc '
         else
            itertype = 'xc'//char (39)
         endif
         write (luout, 5950) iter-1, itertype, etot, dmax 
C
         transfer = .true.
      endif 
      end if        
 5950    format (
     &        t3,i6,1x,a3,1x,f20.10,10x,d20.10)
      
      return
      end

