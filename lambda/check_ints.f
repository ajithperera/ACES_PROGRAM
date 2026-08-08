










      Subroutine Check_ints(Work,Length,Iuhf,null)
      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Length)

      Logical null 
      Data onem,One /-1.0,1.0/
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 


      Irrepx = 1
      If (Ispar .and. Coulomb) Then
         nrow=Isytyp(1,119)
         ncol=Isytyp(2,119)
         Nsize = Idsymsz(Irrepx,Isytyp(1,119),Isytyp(2,119))
         I000 = 1
         I010 = I000 + Nsize
         call getall(Work(I000),Nsize,irrepx,119)
         Call Checksum("119 :",Work(I000),Nsize)
CSS         call output(Work(I000),1,nrow,1,ncol,nrow,ncol,1)
         Nsize = Idsymsz(Irrepx,Isytyp(1,18),Isytyp(2,18))
         I020 = I010 + Nsize
          call getall(Work(i010),Nsize,irrepx,18)
CSS         call output(Work(I010),1,nrow,1,ncol,nrow,ncol,1)
          Call Checksum("21 :",Work(I010),Nsize)
          Call Daxpy(Nsize,Onem,Work(I010),1,Work(I000),1)
          Call Checksum("21 vs 119 :",Work(I000),Nsize)
      Endif 

      Length_23=IDSYMSZ(IRREPX,ISYTYP(1,23),ISYTYP(2,23))
      Length_24=IDSYMSZ(IRREPX,ISYTYP(1,24),ISYTYP(2,24))
      Length_25=IDSYMSZ(IRREPX,ISYTYP(1,25),ISYTYP(2,25))
      Length_26=IDSYMSZ(IRREPX,ISYTYP(1,26),ISYTYP(2,26))
      Length_17=IDSYMSZ(IRREPX,ISYTYP(1,17),ISYTYP(2,17))
      Length_18=IDSYMSZ(IRREPX,ISYTYP(1,18),ISYTYP(2,18))
      Call Getall(Work, Length_23, Irrepx, 23)
      Call checksum("List-23",Work,Length_23)
      Call Getall(Work, Length_24, Irrepx, 24)
      Call checksum("List-24",Work,Length_24)
      Call Getall(Work, Length_25, Irrepx, 25)
      Call checksum("List-25",Work,Length_25)
      Call Getall(Work, Length_26, Irrepx, 26)
      Call checksum("List-26",Work,Length_26)
      if (null) then 
         Write(6,*) "Zeroing 17 and 18"
         Call getall(Work, Length_17, Irrepx, 17)
         Call Putall(Work, Length_17, Irrepx, 147)
         Call getall(Work, Length_18, Irrepx, 18)
         Call Putall(Work, Length_18, Irrepx, 148)

         call dzero(Work,length_17)
         Call Putall(Work, Length_17, Irrepx, 17)
         call dzero(Work,length_18)
         Call Putall(Work, Length_18, Irrepx, 18)
      else 
         Call getall(Work, Length_17, Irrepx, 147)
         Call Putall(Work, Length_17, Irrepx, 17)
         Call Getall(Work, Length_18, Irrepx, 148)
         Call Putall(Work, Length_18, Irrepx, 18)
      endif 
      call checksum("17",work,Length_17)
      call checksum("18",work,Length_18)
       
      Return 
      end 
       

