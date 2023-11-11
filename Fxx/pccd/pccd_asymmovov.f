      Subroutine Pccd_asymmovvo(Vo,Ov,Tmp1,Tmp2,Lenvo,Nocc,Nvrt)

      Implicit Double Precision(A-H,O-Z)
      Logical Symmetry

      Dimension Vo(Lenvo)
      Dimension Ov(Lenvo)
      Dimension Tmp1(Lenvo)
      Dimension Tmp2(Lenvo)



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

      Common /Symm/Symmetry

      Data One,Onem,Ione /1.0D0,-1.0D0,1/

      Irrepx = Ione
      Ispin  = Ione
      Ioff   = Ione

      Write(6,"(a)") "The VO-MO gradient matrices:@-Entry"
      call output(Vo,1,nvrt,1,nocc,nvrt,nocc,1)
      Write(6,"(a)") "The OV-MO gradient matrices"
      call output(Ov,1,nocc,1,nvrt,nocc,nvrt,1)

      Nrow = Nocc
      Ncol = Nvrt
      Call Dcopy(Ncol*Nrow,Vo(Ioff),1,Tmp2(Ioff),1)

      Call Transp(Ov,Tmp1,Nvrt,Nocc)
      Call Daxpy(Ncol*Nrow,Onem,Tmp1(Ioff),1,Vo(Ioff),1)

      Call Transp(Tmp2,Tmp1,Nocc,Nvrt)
      Call Daxpy(Ncol*Nrow,Onem,Tmp1(Ioff),1,Ov(Ioff),1)
      
      Write(6,"(a)") "The VO-MO antisymm. gradient matrices",
     +               "@-Exit"
      call output(Vo,1,nvrt,1,nocc,nvrt,nocc,1)
      Write(6,"(a)") "The OV-MO antisymm. gradient matrices"
      call output(Ov,1,nocc,1,nvrt,nocc,nvrt,1)
      

      Return
      End 


