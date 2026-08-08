













































































































































































































      Subroutine Check_1body_dens(Doo,Dvv,Dvo,Ioo,Ivv,Ivo,Work,
     &                            Maxcor,Nbas,IUhf) 

      Implicit Integer (A-Z)



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



c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 

c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end

      Double Precision Work(Maxcor),Ioo(1),Ivv(1),Ivo(1)
      Double Precision Doo(1),Dvv(1),Dvo(1)
      Double Precision Eoo_aa,Eoo_bb,Evv_aa,Evv_bb,Evo_aa,Evo_bb
      Double Precision scr1(Maxbasfn),Ddot 
      Integer Iscr2(Maxbasfn)
      Logical UHF
C
C This is a routine written for debugging purposes. What this
C is to add delta to one electron terms (fock matrices etc.).

      UHF    = .False.
      UHF    = (Iuhf .EQ. 1)
      Irrepx = 1


C The DOO check 

C The DVV check 

C The DVO check 

C The IOO check 

C The IVV check 

C The IVO check 
      Return
      End 

