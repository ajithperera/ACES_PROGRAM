













































































































































































































      Subroutine Scrnc_form_lrspn_rcc_lists(Work,Memleft,Iuhf,Nperts)

      Implicit Double Precision (A-H, O-Z)
      Integer AAAA_LENGTH_IJAB,BBBB_LENGTH_IJAB,AABB_LENGTH_IJAB
      Integer PHA_Length,PHB_Length,HHA_length,HHB_Length
      Integer PPA_Length,PPB_Length
      Integer Ttypel,Ttyper
      Integer Length(8)
      Logical Ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,Drccd

      Dimension Work(Memleft)

      Common /Extrap/Maxexp,Nreduce,Ntol,Nsizec
      Common /Eominfo/Ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,Drccd



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



c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 

c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end

       Imode = 0

C 490 (Ispin=1,2) keeps Mbar(a,i) elements 

      MAX_PHA_Length = 0
      MAX_PHB_Length = 0
      MAX_PPA_Length = 0
      MAX_PPB_Length = 0
      MAX_HHA_Length = 0
      MAX_HHB_Length = 0

      Do Irrep = 1, Nirrep 
         MAX_PHA_Length = Max(MAX_PHA_Length,Irpdpd(Irrep,9))
         MAX_PHB_Length = Max(MAX_PHB_Length,Irpdpd(Irrep,10))
         MAX_PPA_Length = Max(MAX_PPA_Length,Irpdpd(Irrep,19))
         MAX_PPB_Length = Max(MAX_PPB_Length,Irpdpd(Irrep,20))
         MAX_HHA_Length = Max(MAX_HHA_Length,Irpdpd(Irrep,21))
         MAX_HHB_Length = Max(MAX_HHB_Length,Irpdpd(Irrep,22))
      Enddo 
      
      CALL Updmoi(1,MAX_PHA_length,1,480,0,0)
      CALL Updmoi(1,MAX_PHA_length,1,482,0,0)

      Call Aces_list_memset(1,480,0)
      Call Aces_list_memset(1,482,0)

      If (Iuhf .NE. 0) Then
         CALL Updmoi(1,MAX_PHB_length,2,480,0,0)
         CALL Updmoi(1,MAX_PHB_length,2,482,0,0)

         Call Aces_list_memset(2,480,0)
         Call Aces_list_memset(2,482,0)

      Endif
C 
      Do Irrep = 1, Nirrep 
          Length(Irrep) = Irpdpd(Irrep,9)
          If (Iuhf .NE. 0) Then
             Length(Irrep) = Length(Irrep) + Irpdpd(Irrep,9)
          Endif
          Call Updmoi(Nperts,Length(Irrep),Irrep,373,0,0)
          Call Updmoi(Nperts,Length(Irrep),Irrep,374,0,0)
          Call Updmoi(Nperts,Length(Irrep),Irrep,375,0,0)
          Call Updmoi(Nperts,Length(Irrep),Irrep,376,0,0)
      Enddo
         
      Return
      End 
 
