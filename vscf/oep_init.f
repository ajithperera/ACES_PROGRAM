










      subroutine oep_init(petite_list, iuhf, iwork, imemleft,loslater)

      implicit double precision (a-h, o-z)

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
      
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end


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

      Dimension iwork(imemleft)
      character*3 itertype
      logical petite_list, zmat_exsist, oepint_exsist
      logical loslater

      OEP1    = .FALSE.
      OEP2    = .FALSE.
      OEP1ON  = .FALSE.
      OEP2ON  = .FALSE.
      OEP2ITR = .FALSE.
      R12ON   = .FALSE.


      if (iflags2 (153).eq.3) then
         OEP1 = .TRUE.
         write (6, 5910)
      endif
      if (iflags2 (153).eq.4) then
         OEP2 = .TRUE.
         write (6, 5920)
      endif
      if (OEP1 .or. OEP2) then
c
         if (iuhf.ne.0) then
            write (6, 5901)
            call errex
         endif
         if (petite_list) then
            write (6, 5903)
            call errex
         endif

         i020 = i010 + iintfp * itriln(nirrep+1)
         i030 = i020 + 10000
         i040 = i030 + 500*20
         i050 = i040 + 500*20
         i060 = i050 + 500*1
         if ((i060-i000).gt.maxmem)
     &      call nomem (
     &        'Parse the OEP ASV list', '{OEPINI} <-- VSCF',
     &        i060-i000, maxmem)

         inquire (file = 'ZMAT', exist = zmat_exsist)
         if (.not. zmat_exsist) then
            write (6, 5930)
            call errex
         endif
         luoep = 100
         open (unit=luoep, file='ZMAT', form='formatted', status='old')
         call oepini (
     &        luoep, 10000, 500, 20,
     &        iwork (i020), iwork (i030), iwork (i040), iwork (i050),
     &        loslater)
         close (unit=luoep,status='keep')

         inquire (file = 'OEPINT', exist = oepint_exsist)
         if (.not. oepint_exsist) then
            write (6, 5933)
            call errex
         endif
         open (unit=luoep,file='OEPINT',form='unformatted',status='old')
         call locate (luoep, 'NAUX    ')
         read (unit=luoep) naux
         close (unit=luoep)
         write (6, 5936) naux
c
         if (OEP2) then
            OEPTOL = tol
            tol = 1.d-99
         endif
         itertype (1:3) = '   '
c
 5901    format(
     &        t3,'@VSCF-F: The OEP method for an unrestriced',/,
     &        t3,'reference is not implemented yet.')
 5903    format(
     &        t3,'@VSCF-F: The OEP method for any integral package',/,
     &        t3,'other than VMOL is not implemented yet.')
 5910    format(
     &        t3,'@VSCF-I: The first-order OEP method will be used.',/)
 5915    format(
     &        t3,'@VSCF-I: The exchange OEP-FIA method will be used.',/,
     &        t3,'The local exchange potential will be found from ',/,
     &        t3,'the condition that its occupied-virtual matrix ',/,
     &        t3,'elements are equal to those of non-local exchange',/,
     &        t3,'potential.',/)
5917    format(
     &        t3,'@VSCF-I: The exchange OEP-WFIA method will be used.',
     &        /,t3,'The local exchange potential will be found from ',/,
     &        t3,'the condition that its occupied-virtual matrix ',/,
     &        t3,'elements weighted with the energy denominators are',/,
     &        t3,'equal to those of non-local exchange potential.',/)
 5920    format(
     &        t3,'@VSCF-I: The second-order OEP method will be used.',/)
 5930    format (
     &        t3,'@VSCF-F: The ZMAT file is not found.',/)
 5933    format (
     &        t3,'@VSCF-F: The OEPINT file is not found.',/)
 5936    format (
     &        t3,'@VSCF-I: The auxiliary basis set consists of ',i4,/
     &        t3,'functions.',/)
 5950    format (
     &        t3,i6,1x,a3,1x,f20.10,10x,d20.10)
 5960    format(
     &        t3,1x,a10,f20.10,6x,a3,1x,d20.10/)
 5964    format(
     &        t3,'@VSCF-I: The HOMO condition on the OEP exchange',/,
     &        t3,'potential is violated by ',f10.5,' a.u. (',f10.5,
     &        ' eV).', /)
 5965    format(
     &        t3,'@VSCF-I: The OEP eigenvalues will be shifted by',/,
     &        t3,f10.5,' a.u. to satisfy the HOMO condition on',/,
     &        t3,'the OEP exchange potential.',/)
c
      endif
C
      return
      end

