










      Subroutine Form_vv_ab(T2cc,T2mp,T2,Iwork,Imaxcor,Nab,Irrepab)

      Implicit Double Precision(A-H,O-Z)
      Integer Aend,Bend

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
c active_space.com : begin
      Parameter(Max_xp=100)
      integer Active_oo(8,2),Active_vv(8,2)
      integer Pactive_oo(100,8,2),Pactive_vv(100,8,2)
      integer Ioff_active_oo(8,2),Ioff_active_vv(8,2)
      integer Pioff_active_oo(100,8,2),Pioff_active_vv(100,8,2)
      Double Precision Oo_threshold,Vv_threshold
      common /actvsp_info/Active_oo,Pactive_oo,Active_vv,Pactive_vv,
     +                    Ioff_active_oo,Pioff_active_oo,
     +                    Ioff_active_vv,Pioff_active_vv,
     +                    Oo_threshold,Poo_threshold,
     +                    Vv_threshold,Pvv_threshold,Eta_val(Max_xp),
     +                    E_k(Max_xp),E_ks(Max_xp)
c active_space.com: end

      Dimension Iwork(Imaxcor)
      Dimension T2cc(Nab) 
      Dimension T2mp(Nab) 
      Dimension T2(Nab) 

      Data Ione /1/
    
      Ioff = Ione

      Do Irrep_r = 1, Nirrep

         Irrep_b = Irrep_r
         Irrep_a = Dirprd(Irrep_b,Irrepab)

         Nrow_a  = Vrt(Irrep_a,1)
         Ncol_b  = Vrt(Irrep_b,2)

         Aend    = Active_vv(Irrep_a,1)
         Bend    = Active_vv(Irrep_b,2)
       
         Call Form_vv(T2cc(Ioff),T2mp(Ioff),T2(Ioff),Iwork,Imaxcor,
     +                Nrow_a,Ncol_b,Aend,Bend)

         Ioff = Ioff + Nrow_a*Ncol_b

      Enddo
     
      Return 
      End 
