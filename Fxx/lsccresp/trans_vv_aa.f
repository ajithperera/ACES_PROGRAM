










      Subroutine Trans_vv_aa(T2_in,T2_out,Cvv,Nab,Irrepab,Iuhf,Ispin,
     +                       Type)

      Implicit Double Precision(A-H,O-Z)
      Character*1 A,B,C,D
      Character*3 Type 

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

      Dimension Cvv(Nfea(1)+Iuhf*Nfea(2))
      Dimension T2_in(Nab),T2_out(Nab)

      Data Ione /1/
      Data Done,Dnull /1.0D0,0.0D0/ 


      If (Type .Eq. "C2N") Then
          A = "T"
          B = "N"
          C = "N"
          D = "N"
      elseif (Type .Eq. "N2C") Then
          A = "N"
          B = "N"
          C = "N"
          D = "T"
      Endif 

      Koff = Ione
      Do Irrep = 1, Nirrep
         Irrep_b = Irrep
         Irrep_a = Dirprd(Irrep_b,Irrepab)
 
         Ncol = Vrt(Irrep_b,Ispin)
         Nrow = Vrt(Irrep_a,Ispin)

         Ioff = Ioff_vv(Irrep_a,Ispin) 
         Joff = Ioff_vv(Irrep_b,Ispin) 

         Call Xgemm(A,B,Nrow,Ncol,Nrow,Done,Cvv(Ioff),Nrow,
     +              T2_in(Koff),Nrow,Dnull,T2_out(Koff),Nrow)

         Call Xgemm(C,D,Nrow,Ncol,Ncol,Done,T2_out(Koff),Nrow,
     +              Cvv(Joff),Ncol,Dnull,T2_in(Koff),Nrow)

         Koff = Koff + Nrow*Ncol
      Enddo

      Return
      End

 
