      Subroutine Pccd_Frmful(Htau,Hoo,Hvv,Hvo,Hov,Work,Maxcor,Nocc,
     +                       Nvrt,Nbas,String,Sym_packed)

       Implicit Double Precision(A-H,O-Z)
       Character*7 String
       Logical Sym_packed 
       Logical Symmetry

       Dimension Htau(Nbas,Nbas)
       Dimension Hoo(Nocc*Nocc)
       Dimension Hvv(Nvrt*Nvrt)
       Dimension Hvo(Nvrt*Nocc)
       Dimension Hov(Nvrt*Nocc)
       Dimension Work(Maxcor)



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
      Common/Symm/Symmetry

      Data Ione, Izero /1,0/

      Irrepx   = Ione
      If (Sym_packed) Then 
          Lenoo    = Irpdpd(Irrepx,21)
          Lenvv    = Irpdpd(Irrepx,19)
          Lenvo    = Irpdpd(Irrepx,9)
      Else
          Lenoo    = Nocc*Nocc
          Lenvv    = Nvrt*Nvrt
          Lenvo    = Nvrt*Nocc
      Endif 

      Call Pccd_putblocks(Htau,Hoo,Hvv,Hvo,Hov,Work,Maxcor,Lenoo,
     +                    Lenvv,Lenvo,Nocc,Nvrt,Nbas,String,
     +                    Sym_packed)

      Return 
      End 


       
