










      Subroutine Pccd_check_hess(Work,Maxcor)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Logical Print

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

      Irrepx = 1
      Print  = .False.
      Write(6,*)
      List_h = 213
      Nsize = Idsymsz(Irrepx,21,21)
      Call Getall(Work,Nsize,Irrepx,List_h)
      Call checksum("H(Ij,Kl):",Work,Nsize)
      If (Print) Call Pccd_print_all(Work,Nsize,"ijkl",Irrepx)
    
      Print  = .False.
      List_h = 234
      Nsize = Idsymsz(Irrepx,19,19)
      Call Getall(Work,Nsize,Irrepx,List_h)
      Call checksum("H(Ab,cD):",Work,Nsize)
      If (Print) Call Pccd_print_all(Work,Nsize,"abcd",Irrepx)

      List_h = 217
      Nsize = Idsymsz(Irrepx,14,15)
      Call Getall(Work,Nsize,Irrepx,List_h)
      Call checksum("H(Ij,Ab):",Work,Nsize)
   
      List_h = 216
      Nsize = Idsymsz(Irrepx,15,14)
      Call Getall(Work,Nsize,Irrepx,List_h)
      Call checksum("H(Ab,iJ):",Work,Nsize)

C These are not needed for orbital optimizations

      Return 
      End 
