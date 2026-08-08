










      Subroutine Pccd_zero_hess(Work,Maxcor)

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
      Print  = .True.
      Write(6,*)
      List_h = 213
      Nsize = Idsymsz(Irrepx,14,14)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)
    
      Print  = .False.
      List_h = 234
      Nsize = Idsymsz(Irrepx,15,15)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)

      List_h = 217
      Nsize = Idsymsz(Irrepx,14,15)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)
   
      List_h = 216
      Nsize = Idsymsz(Irrepx,15,14)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)

      List_h = 210
      Nsize = Idsymsz(Irrepx,14,18)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)

      List_h = 211
      Nsize = Idsymsz(Irrepx,18,14)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)

      List_h = 230
      Nsize = Idsymsz(Irrepx,15,11)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)

      List_h = 231
      Nsize = Idsymsz(Irrepx,11,15)
      Call Dzero(Work,Nsize)
      Call Putall(Work,Nsize,Irrepx,List_h)

      Return 
      End 
