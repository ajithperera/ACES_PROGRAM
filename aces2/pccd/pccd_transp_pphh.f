










      Subroutine Pccd_transp_pphh(Pphh,Hhpp,Nsize)

      Implicit Double Precision(A-H,O-Z)
      Dimension Pphh(Nsize)
      Dimension Hhpp(Nsize)

c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end

      Data One,Onem,Dnull,Ione /1.0D0,-1.0D0,0.0D0,1/

      Irrepx = 1

C H(ab,ij) -> H(ij,ab)

      Ioff = Ione 
      Do Irrep_ij = 1, Nirrep
         Irrep_ab = Dirprd(Irrep_ij,Irrepx)

         Nrow_ab = Irpdpd(Irrep_ab,15)
         Ncol_ij = Irpdpd(Irrep_ij,14)
         Call Transp(Pphh(Ioff),Hhpp(Ioff),Ncol_ij,Nrow_ab)
         Ioff = Ioff + Nrow_ab*Ncol_ij

      Enddo 

      Return
      End
