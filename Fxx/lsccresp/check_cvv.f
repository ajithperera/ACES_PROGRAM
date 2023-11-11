










      Subroutine check_cvv(Cvv,Iuhf)

      Implicit Double Precision (A-H,O-Z)

c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end


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



c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end
      
      Dimension Cvv(Nfea(1)+Iuhf*Nfea(2))
     
      Do Ispin = 1, Iuhf+1 
         Ibegin = 1 + Nfea(1)*(Ispin-1)
         Ioff = 0
         Write(6,"(a)") " The CVV matrix"
      Do Irrep = 1, Nirrep
         Ndim = Vrt(Irrep,Ispin) 
         Call output(Cvv(Ioff+Ibegin),1,Ndim,1,Ndim,Ndim,Ndim,1)
         Ioff = Ioff + Ndim*Ndim
      Enddo 
      Enddo

      Return
      End 
