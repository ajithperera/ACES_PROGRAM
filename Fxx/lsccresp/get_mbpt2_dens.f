










      Subroutine Get_mbpt2_dens(T2amps,Doo,Dvv,Dov,Work,Maxcor,T2ln,
     &                          T2ln_aa,T2ln_bb,T2ln_ab,Iuhf)

      Implicit Double Precision (A-H,O-Z)
      Integer T2ln,T2ln_aa,T2ln_bb,T2ln_ab



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
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end

      Dimension T2amps(T2ln)
      Dimension Work(Maxcor)
      Dimension Doo(Nfmi(1)+Iuhf*Nfmi(2))
      Dimension Dvv(Nfea(1)+Iuhf*Nfea(2))
      Dimension DOv(Nt(1)+Iuhf*Nt(2))

      Call Dzero(Doo,Nfmi(1)+Iuhf*Nfmi(2))
      Call Dzero(Dvv,Nfea(1)+Iuhf*Nfea(2))

      Call Densoo_mp2(Doo,T2amps,Work,Maxcor,T2ln,T2ln_aa,T2ln_bb,
     +                T2ln_ab,Iuhf)

      Call Densvv_mp2(Dvv,T2amps,Work,Maxcor,T2ln,T2ln_aa,T2ln_bb,
     +                T2ln_ab,Iuhf)

      Return
      End 
