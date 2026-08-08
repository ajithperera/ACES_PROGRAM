










      subroutine molden_corsorb(ener,iocc,orb,orb_r,iang,scr,nao,nmo,
     &                          maxcor,iuhf,iexx,iunit,iroot)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      double precision ener(nmo),orb(nao,nmo),scr(maxcor),
     &                 orb_r(nao)

      integer iocc(nmo), iang(nao)
      integer  idrppop(8),idrpvrt(8)
      logical ACESIII



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



c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end

      parameter (one=1.0D0)
      parameter (zilch=0.0D0)
      character*2 iroot
      character*5 sptype(2)
      character*8 string
      character*8 cscfener(2)
      data sptype   /'Alpha','Beta '/
      data cscfener /'SCFEVLA0','SCFEVLB0'/
      call aces_com_info
      call aces_com_syminf
      call aces_com_sym
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      i000=1
      i010=i000+nao*nmo
      i020=i010+nmo**2

      if(i020.gt.maxcor/iintfp)
     &  call insmem('molden_corsorb.F',i020*iintfp,maxcor)
C
C Note that corresponding orbitals are already transformed to the 
C ZMAT ordered AO basis.
C
      do ispin=1,iuhf+1

         If (Ispin .EQ. 1) String = "CRORBA"//iroot
         If (Ispin .EQ. 2) String = "CRORBB"//iroot
 
         Write(6,*) String
         call getrec(20,'JOBARC',String,nmo*nao*iintfp,orb)
C
         call getrec(20, "JOBARC", 'OCCUPYA0', 1, POP(1,1))
         Vrt(1,1) = NMO - POP(1,1)
         Write(6,*) pop(1,1), NMO
         If (iuhf.eq.1) then
            call getrec(20, "JOBARC", 'OCCUPYB0', 1, POP(2,1))
            Vrt(2,1) = NMO - POP(2,1)
         Endif
C
C Note the these orbital energies and occupations numbers are 
C based on Koopman's 
C
      call getrec(-1,'JOBARC',cscfener,nmo*iintfp,ener)
      do irrep =1, nirrep
        nocco(ispin)=nocco(ispin)+pop(irrep,ispin)
      enddo

      Enddo

      do 80 imo=1,nmo
          write(iunit,90)' Ene=   ',ener(imo)
90        format(A,F8.4)
          write(iunit,100)' Spin= ' // sptype(ispin)
100       format(A)
          write(iunit,110)'  Occup=   ',iocc(imo),'.000000'
110       format(A,I1,A)
C
C Reorder the F functions (perhaps g too ..), Peter Szalay identified
C the problem and proivided the reorder routines. 
C Ajith Perera, 06/2011.
C
          call getrec(-1, 'JOBARC', 'ANMOMBF0', nao, iang)
          call reorder(orb(1, imo), orb_r, iang, nao)
          do 120 iao=1,nao
CSSS             write(iunit,130)iao,'  ',orb(iao,imo)
            write(iunit,130)iao,'  ',orb_r(iao)
130         format(I4,A,F12.6)
120       continue
80      continue 
10    continue

      return
      end




