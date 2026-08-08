













































































































































































































      Subroutine Rcc_expand(Win,Wout,Ncol_squar_ij,Nrow_squar_ab,
     +                      ncol_tring_ij,Nrow_tring_ab,
     +                      Irrep,Ispin,Work,Left)

      Implicit Double Precision (A-H, O-Z)
      Dimension Work(Left),Win(ncol_tring_ij*Nrow_tring_ab)
      Dimension Wout(Ncol_squar_ij*Nrow_squar_ab)



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

c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end

      Irrepl = Irrep
      Irrepr = Irrep

      I000 = 1
      I010 = I000 + Nrow_squar_ab*Ncol_squar_ij
      I020 = I010 + Nrow_squar_ab*Ncol_squar_ij
      Iend = I020
      IF (Iend .Gt. Left) Call Insmem("rcc_expand",Iend,Left)

C D(A<=B,I<=J) -->D(AB,I<=J)

      Call Symexp6(Irrepl,Vrt(1,ispin),Vrt(1,ispin),
     +             Nrow_squar_ab,Nrow_tring_ab,
     +             ncol_tring_ij,work(I000),Win,
     +             work(I010))
      Call Transp(work(I000),work(I010),ncol_tring_ij,
     +            nrow_squar_ab)
C D(IJ,A<=B) -->D(iJ,AB)

      Call Symexp6(Irrepr,Pop(1,ispin),Pop(1,ispin),
     +             Ncol_squar_ij,Ncol_tring_ij,
     +             nrow_squar_ab,work(I000),work(I010),
     +             Wout)
C D(IJ,AB) -->D(AB,IJ)

       Call Transp(work(I000),Wout,Nrow_squar_ab,
     +                  Ncol_squar_ij)
      Return
      End

