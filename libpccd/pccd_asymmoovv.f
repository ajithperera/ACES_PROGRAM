










      Subroutine Pccd_asymmoovv(Oo,Vv,Oot,Vvt,Lenoo,Lenvv,Nocc,
     +                          Nvrt)

      Implicit Double Precision(A-H,O-Z)

      Dimension Oo(Lenoo)
      Dimension Oot(Lenoo)
      Dimension Vv(Lenvv)
      Dimension Vvt(Lenvv)



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
     
      Data Onem,Ione /-1.0D0,1/

      Irrepx = Ione
      Ispin  = Ione
      Ioff   = Ione
      Call Dzero(Vvt,Lenvv)
      Call Dzero(OOt,Lenoo)

      Nrow = Nocc
      Ncol = Nocc
      Call Transp(Oo(Ioff),Oot(Ioff),Ncol,Nrow)
      Call Daxpy(Ncol*Nrow,Onem,Oot(Ioff),1,Oo(Ioff),1)

      Nrow = Nvrt
      Ncol = Nvrt
      Call Transp(Vv(Ioff),Vvt(Ioff),Ncol,Nrow)
      Call Daxpy(Ncol*Nrow,Onem,Vvt(Ioff),1,Vv(Ioff),1)


      write(6,*)
      Write(6,"(a)") "The OO-MO antisymm. gradient matrices"
      call output(OO,1,nocc,1,nocc,nocc,nocc,1)
      Write(6,"(a)") "The VV-MO antisymm. gradient matrices"
      call output(VV,1,nvrt,1,nvrt,nvrt,nvrt,1)
      Return
      End 


