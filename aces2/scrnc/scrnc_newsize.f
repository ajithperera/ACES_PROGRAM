
























































































































































































































      Subroutine Scrnc_newsize(Work,Memleft,Iuhf,Irrepx,Nsizec)

      Implicit Double Precision (A-H, O-Z)

      Dimension Work(Memleft)



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

      Nsizec = Irpdpd(Irrepx,9)
      If (Iuhf .EQ. 0) Then
         Nsizec = Nsizec  + Idsymsz(Irrepx,13,14)
      Else
        Nsizec = Nsizec + Irpdpd(Irrepx,10)
        Nsizec = Nsizec + Idsymsz(Irrepx,1,3)
        Nsizec = Nsizec + Idsymsz(Irrepx,2,4)
        Nsizec = Nsizec + Idsymsz(Irrepx,13,14)
      Endif

      Return
      End
