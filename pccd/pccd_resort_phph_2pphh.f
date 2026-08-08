










      Subroutine Pccd_resort_phph_2pphh(Work,Maxcor,List_s,List_t)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

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
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end

      Data One,Onem,Dnull,Ione /1.0D0,-1.0D0,0.0D00,1/
      Print*, "in",List_s,List_t

      Ispin   = Ione
      Irrepx  = Ione
      Nsize_s = Idsymsz(Irrepx,11,11)
      Nsize_t = Idsymsz(Irrepx,15,14)

      I000 = Ione
      I010 = I000 + Nsize_s
      I020 = I010 + Nsize_t
      Iend = I020 + Nocco(Ispin)*Nvrto(Ispin)
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_resort_phph_2pphh",
     +                                   Iend,Maxcor)
      Call Getall(Work(I000),Nzise_s,Irrepx,List_s)
      Call SStgen(Work(I000),Work(I010),Nsize_s,Vrt(1,Ispin),
     +             Pop(1,Ispin),Vrt(1,Ispin),Pop(1,Ispin),Work(I020),
     +             Irrepx,"1324")
      Call Getall(Work(I010),Nzise_t,Irrepx,List_t)

      Return
      End
       
     
