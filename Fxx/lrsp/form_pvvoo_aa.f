










      Subroutine Form_pvvoo_aa(T2cc,T2mp,T2,Iwork,Imaxcor,Nij,Nab,
     +                         Irrepij,Irrepab,Ispin,Ipert)

      Implicit Double Precision(A-H,O-Z)
      Integer A,B,Aend,Bend

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

      Dimension T2cc(Nab,Nij),T2mp(Nab,Nij),T2(Nab,Nij)

      Data Ione,Inull /1,0/
      Data Done,Dnull /1.0D0,0.0D0/

C Note the both type of T2 amplitudes are stored as AB,IJ form 

      Ioff    = Ione 
      Nij_b   = Nij

      Call Dcopy(Nab*Nij,T2mp,1,T2,1)

      Do Irrep_r = 1, Nirrep

         Irrep_j = Irrep_r
         Irrep_i = Dirprd(Irrep_j,Irrepij)

         Ncol_j  = Pop(Irrep_j,Ispin)
         Nrow_i  = Pop(Irrep_i,Ispin)
        
         Iend    = Nrow_i - Pactive_oo(Ipert,Irrep_i,Ispin) + Ione
         Jend    = Ncol_j - PactiVe_oo(Ipert,Irrep_j,Ispin) + Ione

C Select active ij pairs

         Do J = Ncol_j, Jend, -1
            Do I = Nrow_i, Iend, -1 

               Nij_a =  (I-Nrow_i)+(J-Ncol_j)*Nrow_i
               Nij_c  = Nij_b + Nij_a
   
               Call Form_pvv_aa(T2cc(Ioff,Nij_c),T2mp(Ioff,Nij_c),
     +                          T2(Ioff,Nij_c),Iwork,Imaxcor,Nab,
     +                          Irrepab,Ipert,Ispin)
            Enddo
         Enddo
         Nij_b = Nij_b - Ncol_j*Nrow_i
      Enddo

      Call Dcopy(Nab*Nij,T2,1,T2cc,1)

      Return
      End
