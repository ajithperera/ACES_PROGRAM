










      subroutine oep_post(iuhf, nbas, ibufln, luint, nocc, icore, 
     &                    maxmem, scfevl, scfden, etot, dmax, 
     &                    itertype)
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
      integer
     &     nbas, iuhf, maxmem, nocc(16), ibufln, icore(*), i000,
     &     i010, i020, i030, i040, i050, i060, ione, luint
       integer i
      double precision au2ev, etot, dmax, homo, 
     &                 scfevl (nbas), scfden (itriln (nirrep+1))

      character*1 asv_null, itertype

      parameter (au2ev = 27.2113957D0)
c 
c     OEP post-SCF block: <> Print out the converged OEP energy. <> If
c     local exchange potential is used: Report by how much the HOMO
c     condition on the OEP exchange potential is violated. Shift the OEP
c     eigenvalues to satisfy the HOMO condition if no auxiliary exchange
c     potential is used. <> Allocate memory for scratch arrays necessary
c     to calculate the Fock matrix. <> Check if there is enough memory
c     available. <> Calculate the Fock matrix and write it over the OEP
c     Hamiltonian, in case a conventional correlation method will be
c     used with OEP orbitals. <> Calculate the components of the total
c     energy and report them. (WARNING: this requires reading the
c     integral file again!)
c
c     Note that orbital energies for ROHF are calculated somewhere
c     between the final Fock matrix construction and energy
c     calculation. So it is bypassed if the OEP method is invoked. If
c     one wants to do ROHF-like OEP method, one has to fix this.
c
      asv_null=achar(0)
C
            do i=1,nbas
               write(6,*) i,scfevl(i)
            end do
C
      if (OEP1 .or. OEP2) then
c
         write (6, 5960) 'E (OEP) = ', etot, itertype
c
         call getrec(20, 'JOBARC', 'OEPXHOMO', iintfp, homo)
         write (6, 5964) homo, homo*au2ev
         if (OEPAXP (1:5) .eq. 'NONE'//asv_null) then
             call daxpy (nbas, 1.d0, homo, 0, scfevl, 1)
             write (6, 5965) homo
         endif
            write(6,*) 'after hmo cont'
            do i=1,nbas
               write(6,*) i,scfevl(i)
            end do 
         call putrec(20, 'JOBARC', 'SCFEVLA0', nbas*iintfp, scfevl)
C
         ione = 1
         i010 = ione
         i020 = i010 + (iuhf+1) * isqrln (nirrep+1) * iintfp
         i030 = i020 + (iuhf+1) * isqrln (nirrep+1) * iintfp
         i040 = i030 + mod((i020-i000),iintfp)
         i050 = i040 + ibufln * iintfp
         i060 = i050 + ibufln 
         call oepenrcom (
     &        iflags2 (153), nbas, nocc (1), luint, ibufln,
     &        icore (i050), icore (i040),
     &        scfevl, scfden) 
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
C
       
C --- Make the standard from the main program vscf.F
c
CSSS         call mkrhff(
CSSS     &        icore(i040), icore(i030), icore(i010),
CSSS     &        icore(i055), icore(i056),
CSSS     &        itriln(nirrep+1), nbas, nbfirr,
CSSS     &        icore(i057), ibufln, luint, .true.)
c
cSSS         call oepfiachk (
cSSS     &        nbas, nocc (1), naux,
cSSS     &        icore (i030), maxmem-i080+i000, i080-i000,
cSSS     &        icore (oepauxao),
cSSS     &        scfevl, icore (i010), icore (i020))
c
CSSS         goto 5999
      endif

      return
      end

