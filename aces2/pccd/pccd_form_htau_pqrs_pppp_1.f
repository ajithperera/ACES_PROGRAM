










      Subroutine Pccd_form_htau_pqrs_pppp_1(Work,Maxcor,List_v,
     +                                      List_g,List_h)
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

      Data One,Onem,Dnull,Ione,Two,Four,Twom/1.0D0,-1.0D0,0.0D0,1,
     +                                       2.0D0,4.0D0,-2.0D0/
      Irrepx = Ione
      Ispin  = Ione

      Nsize_t = Idsymsz(Irrepx,15,15)
      I000 = Ione
      I010 = I000 + Nsize_t
      Call Dzero(Work(I000),Nsize_t)

C <be|cf>*G(ae,df) + <eb|cf>*G(ea,df) -> H(bc,ad)
C <df|ae>*G(be,cf) + <df|ea>*G(be,fc) -> H(ad,bc)

      Nsize = Idsymsz(Irrepx,15,15)

C Spin-adapt and then <be|cf> -> <bc|ef>

      I020  = I010 + Nsize
      I030  = I020 + Nsize
      Iend  = I030 + Max(Nvrto(Ispin)*Nvrto(Ispin),
     +                   Nocco(Ispin)*Nocco(Ispin),
     +                   Nocco(Ispin)*Nvrto(Ispin))
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pppp_1",
     +                                   Iend,Maxcor)
      Joff = I020
      Do Irrep_cf = 1, Nirrep
         Irrep_be = Dirprd(Irrep_cf,Irrepx)

         Nrow_be = Irpdpd(Irrep_be,15)
         Ncol_cf = Irpdpd(Irrep_cf,15)

         I040 = I030 + Max(Nrow_be,Ncol_cf)
         Iend = I040 + Max(Nrow_be,Ncol_cf)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pppp_1",
     +                                      Iend,Maxcor)
         Call Getlst(Work(Joff),1,Ncol_cf,1,Irrep_cf,List_v)
         Call Spinad3(Irrep_be,Vrt(1,1),Nrow_be,Ncol_cf,Work(Joff),
     +                Work(I030),Work(I040))
         Joff = Joff + Nrow_be*Ncol_cf
      Enddo

      Call Sstgen(Work(I020),Work(I010),Nsize,Vrt(1,Ispin),
     +            Vrt(1,Ispin),Vrt(1,Ispin),Vrt(1,Ispin),Work(I030),
     +            Irrepx,"1324")

C G(ea,fd) -> G(ef,ad)

      I030  = I020 + Nsize
      I040  = I030 + Nsize
      Iend  = I040 + Max(Nvrto(Ispin)*Nvrto(Ispin),
     +                   Nocco(Ispin)*Nocco(Ispin),
     +                   Nocco(Ispin)*Nvrto(Ispin))
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pppp_1",
     +                                      Iend,Maxcor)
      Call Getall(Work(I030),Nsize,Irrepx,List_g)
      Call Sstgen(Work(I030),Work(I020),Nsize,Vrt(1,Ispin),
     +            Vrt(1,Ispin),Vrt(1,Ispin),Vrt(1,Ispin),Work(I040),
     +            Irrepx,"1324")

C <bc|ef>*G(ef,ad) -> H(bc,ad)

      Ioff = I000
      Joff = I010
      Koff = I020
      Do Irrep_ad = 1, Nirrep
         Irrep_ef = Dirprd(Irrep_ad,Irrepx)
         Irrep_bc = Dirprd(Irrep_ef,Irrepx)
         
         Nrow_ef = Irpdpd(Irrep_ef,15)
         Nrow_bc = Irpdpd(Irrep_bc,15)
         Ncol_ad = Irpdpd(Irrep_ad,15)
  
         Iend = I030 + Nrow_bc*Ncol_ad

         Icheck = Min(Nrow_ef,Nrow_bc,Ncol_ad)
         If (Icheck .Gt. 0) Then
            Call Dgemm("N","N",Nrow_bc,Ncol_ad,Nrow_ef,Twom,Work(Joff),
     +                  Nrow_bc,Work(Koff),Nrow_ef,One,Work(Ioff),
     +                  Nrow_bc)
         Endif 

C H(bc,ad) + H(bc,ad)^t

         Call Transp(Work(Ioff),Work(I030),Ncol_ad,Nrow_bc)
         Call Daxpy(Ncol_ad*Nrow_bc,One,Work(I030),1,Work(Ioff),1)

         Ioff = Ioff + Nrow_bc*Ncol_ad
         Joff = Joff + Nrow_bc*Nrow_ef
         Koff = Koff + Ncol_ad*Nrow_ef
      Enddo 

      call checksum("VVVV-C  :",Work(I000),Nsize_t)

      I020 = I010 + Nsize_t
      Call Getall(Work(I010),Nsize_t,Irrepx,List_h)
      Call Daxpy(Nsize_t,One,Work(I010),1,Work(I000),1)
      Call Putall(Work(I000),Nsize_t,Irrepx,List_h)

      Return
      End 


