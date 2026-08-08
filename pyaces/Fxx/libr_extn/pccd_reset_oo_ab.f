      Subroutine Pccd_reset_oo_ab(T2_in,Nij,Irrepab)

      Implicit Double Precision(A-H,O-Z)

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
c symoff.com : begin
      integer  Ioff_oo(8,2),Ioff_vv(8,2)
      common /symoff/ Ioff_oo,Ioff_vv
c symoff.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end

      Dimension T2_in(Nij)

      Data Ione,Dnull,Inull /1,0.0D0,0/

      Ioff = Ione
      Do Irrep = 1, Nirrep
         Irrep_b = Irrep
         Irrep_a = Dirprd(Irrep_b,Irrepab)
         Ncol = Pop(Irrep_b,1)
         Nrow = Pop(Irrep_a,1)
         If (Irrep_a .Eq. Irrep_b) Then 
             Call Pccd_set2_zero(T2_in(Ioff),Nrow,Ncol)
             Ioff = Ioff + Nrow*Ncol
         Else
             Call Dzero(T2_in(Ioff),Nrow*Ncol)
             Ioff = Ioff + Nrow*Ncol
         Endif 
      Enddo

      Return
      End

 
