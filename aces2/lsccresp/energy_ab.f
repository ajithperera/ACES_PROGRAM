










      Subroutine Energy_ab(T1aa,T1bb,T2ab,W2ab,Ndim_a,Ndim_b,Nrow,Ncol,
     +                     Listw,Irrep,Tau,E,ES)

      Implicit Double Precision(A-H,O-Z)
      Logical tau

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
     
      Dimension T1aa(Ndim_a)
      Dimension T1bb(Ndim_b)
      Dimension T2ab(Nrow,Ncol)
      Dimension W2ab(Nrow,Ncol) 
      Dimension E(3),ES(3)
      Data One /1.0D0/

      Call Getlst(W2ab,1,Ncol,1,Irrep,Listw)

      ES(3) = ES(3) + Ddot(Nrow*Ncol,W2ab,1,T2ab,1)

      If (Tau) then
         Call Ftau(T2ab,T1aa,T1bb,Nrow,Ncol,Pop(1,1),Pop(1,2),
     +             Vrt(1,1),Vrt(1,2),Irrep,3,One)
      Endif 

      call checksum("Wab  :",W2ab,Ncol*Nrow)
      call checksum("T2ab :",T2ab,Ncol*Nrow)
      E(3)  = E(3) + Ddot(Nrow*Ncol,W2ab,1,T2ab,1)

      Return
      End
      
