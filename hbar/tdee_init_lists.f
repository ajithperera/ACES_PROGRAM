










      Subroutine tdee_init_lists(work,Memleft,Iuhf)

      Implicit Integer(A-Z)

      Dimension Work(Memleft)
      Dimension Target_lists(2),Source_lists(2)
      Logical Source, Target



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
   
      Irrepx = 1

      pha_length = Irpdpd(Irrepx,9)
      phb_length = Irpdpd(Irrepx,10)
      hha_length = Irpdpd(Irrepx,21)
      hhb_length = Irpdpd(Irrepx,22)
      ppa_length = Irpdpd(Irrepx,19)
      ppb_length = Irpdpd(Irrepx,20)

      Max_s_length = Max(pha_length,phb_length,hha_length,hhb_length,
     +                   ppa_length,ppb_length)
     
      Call DZero(Work,Max_s_length)

      Call Putlst(Work,1,1,1,1,90)
      Call Putlst(Work,1,1,1,2,90)

C  Other auxilary doubles lists

      If (Iuhf .EQ. 0) Then
        Call Zerolist(Work,Memleft,34)
        Call Zerolist(Work,Memleft,37)
        Call Zerolist(Work,Memleft,39)
        Call Zerolist(Work,Memleft,42)
        Call Zerolist(Work,Memleft,43)
        Call Zerolist(Work,Memleft,44)
        Call Zerolist(Work,Memleft,46)
      Else
        Call Zerolist(Work,Memleft,34)
        Call Zerolist(Work,Memleft,35)
        Call Zerolist(Work,Memleft,36)
        Call Zerolist(Work,Memleft,37)
        Call Zerolist(Work,Memleft,38)
        Call Zerolist(Work,Memleft,39)
        Call Zerolist(Work,Memleft,40)
        Call Zerolist(Work,Memleft,41)
        Call Zerolist(Work,Memleft,42)
        Call Zerolist(Work,Memleft,43)
        Call Zerolist(Work,Memleft,44)
        Call Zerolist(Work,Memleft,45)
        Call Zerolist(Work,Memleft,46)
      Endif 

      Return
      End

    
