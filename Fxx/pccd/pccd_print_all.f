      Subroutine Pccd_print_all(Work,Nsize,Type,Irrepx)

      Implicit Double Precision(A-H,O-Z)
      Dimension Work(Nsize) 
      Character*4 Type

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

      Ioff = 1
      Do Irrep_r = 1, Nirrep
         Irrep_l = Dirprd(Irrep_r,Irrepx)
         If (Type .EQ. "ijkl") Then
            Nrow = Irpdpd(Irrep_l,21)
            Ncol = Irpdpd(Irrep_r,21)
            Call output(Work(Ioff),1,Nrow,1,Ncol,Nrow,Ncol,1)
	    Ioff = Ioff + Nrow*Ncol 
         Endif
         If (Type .EQ. "abcd") Then
            Nrow = Irpdpd(Irrep_l,19)
            Ncol = Irpdpd(Irrep_r,19)
            Call output(Work(Ioff),1,Nrow,1,Ncol,Nrow,Ncol,1)
	    Ioff = Ioff + Nrow*Ncol 
         Endif
      Enddo

      Return
      End 

      

