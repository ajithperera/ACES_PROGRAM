










      Subroutine Energy_aa(T1aa,T1bb,T2aa,W2aa,Ndim_a,Ndim_b,Ncol,Nrow,
     +                     Listw,Irrep,Ispin,Tau,E,ES)

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
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end
     
      Dimension T1aa(Ndim_a)
      Dimension T1bb(Ndim_b)
      Dimension T2aa(Nrow,Ncol)
      Dimension W2aa(Nrow,Ncol) 
      Dimension E(3),ES(3)
      Data One /1.0D0/

      Call Getlst(W2aa,1,Ncol,1,Irrep,Listw)

      If (Ispin .Eq. 1) Then
         ES(1) = ES(1) + Ddot(Nrow*Ncol,W2aa,1,T2aa,1)
      ELseif (Ispin .Eq. 2) then
         ES(2) = ES(2) + Ddot(Nrow*Ncol,W2aa,1,T2aa,1)
      Endif 

      If (Tau) then
         Call Ftau(T2aa,T1aa,T1bb,Nrow,Ncol,Pop(1,Ispin),
     +             Pop(1,Ispin),Vrt(1,Ispin),Vrt(1,Ispin),
     +             Irrep,Ispin,One) 
      Endif 

      Print*, "After Tau"
      call checksum("Waa  :",W2aa,Ncol*Nrow)
      call checksum("T2aa :",T2aa,Ncol*Nrow)

      If (Ispin .Eq. 1) Then
         E(1)  = E(1) + Ddot(Nrow*Ncol,W2aa,1,T2aa,1)
      ELseif (Ispin .Eq. 2) then
         E(2)  = E(2) + Ddot(Nrow*Ncol,W2aa,1,T2aa,1)
      Endif 
      
      Return
      End
      
