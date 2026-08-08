










      Subroutine Add_external_t3_to_h12(Work,Length,Iuhf,Non_hf)

      Implicit Double Precision (A-H,O-Z)
      Logical Non_hf,Uhf

      Dimension Work(Length)

      Data Ione,onem,One /1,-1.0D0,1.0D0/

c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
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


       Write(6,"(a)") "-----Entered add_external_t3_to_h12-----"
       Write(6,"(a,l)") " The NON-HF is set to: ", Non_hf

      Irrepx = 1
      Uhf    = (Iuhf .NE. 0)
      Imode = 0

      Lenhp_aa = Irpdpd(Irrepx,9)
      Lenhp_bb = Irpdpd(Irrepx,10)
      Call Updmoi(Irrepx,Lenhp_aa,5,93,0,0)
      Maxlen = Max(Lenhp_aa,Lenhp_bb)
      I000 = Ione
      Iend = I000 + Maxlen
      If (Iend .Gt. Length)Call Insmem ("add_external_t3_to_h12",Maxlen,
     +                                  Length)
      Call Getlst(Work(I000),1,1,1,3,90)
      Call Putlst(Work(I000),1,1,1,5,93)

      If (Iuhf .Ne. 0) Then
         Call Updmoi(Irrepx,Lenhh_bb,6,93,0,0)
         Call Getlst(Work(I000),1,1,1,4,90)
         Call Putlst(Work(I000),1,1,1,6,93)
      Endif 

      Call Inipck(Irrepx,1,3,147,Imode,0,1)
      Call Inipck(Irrepx,13,14,149,Imode,0,1)
      If (Uhf) Call Inipck(Irrepx,2,4,148,Imode,0,1)

C The <ab|ij> contributions from T3 are temporarliy stored in
C list 61-63 (T2 residual list) copy them to 147-149

      Length_61=Idsymsz(IRREPX,ISYTYP(1,61),ISYTYP(2,61))
      If (Iuhf .NE. 0) Length_62=Idsymsz(IRREPX,ISYTYP(1,62),
     +                                   ISYTYP(2,62))
      Length_63=Idsymsz(IRREPX,ISYTYP(1,63),ISYTYP(2,63))

      Maxln = Max(Length_61,Length_62,Length_63)
      I000 = Ione
      I010 = I000 + Maxln
      Iend = I010 + Maxln
      If (Iend .Gt. Length) Call Insmem("add_external_t3_to_h12",
     +                                   Iend,Length)

      If (Uhf) Then
         Call Getall(Work(I000),Length_61,Irrepx,61)
         Call Putall(Work(I000),Length_61,Irrepx,147)
         Call Getall(Work(I000),Length_62,Irrepx,62)
         Call Putall(Work(I000),Length_62,Irrepx,148)
      Endif

      Call Getall(Work(I000),Length_63,Irrepx,63)
      Call Putall(Work(I000),Length_63,Irrepx,149)

C The <ij|ka> and <ab|ci> contributions from T3 are 3 are temporarliy 
C stored in 107-110 and 127-130. Add them to respective integrals lists
C in  7-10 and 27-30.

      If (Uhf) then
         Length_7=IDSYMSZ(IRREPX,ISYTYP(1,7),ISYTYP(2,7))
         Length_8=IDSYMSZ(IRREPX,ISYTYP(1,8),ISYTYP(2,8))
         Length_9=IDSYMSZ(IRREPX,ISYTYP(1,9),ISYTYP(2,9))
         Length_10=IDSYMSZ(IRREPX,ISYTYP(1,10),ISYTYP(2,10))

         Maxln = Max(Length_7,Length_8,Length_9,Length_10)
         I000 = Ione
         I010 = I000 + Maxln
         Iend = I010 + Maxln
         If (Iend .Gt. Length) Call Insmem("add_external_t3_to_h12",
     +                                      Iend,Length)
  
         Call Getall(Work(I000),Length_7,Irrepx,107)
         Call Getall(Work(I010),Length_7,Irrepx,7)
         Call Daxpy(Length_7,One,Work(I000),1,Work(I010),1)
         Call Putall(Work(I010),Length_7,Irrepx,7)
         
         Call Getall(Work(I000),Length_8,Irrepx,108)
         Call Getall(Work(I010),Length_8,Irrepx,8)
         Call Daxpy(Length_8,One,Work(I000),1,Work(I010),1)
         Call Putall(Work(I010),Length_8,Irrepx,8)

         Call Getall(Work(I000),Length_9,Irrepx,109)
         Call Getall(Work(I010),Length_9,Irrepx,9)
         Call Daxpy(Length_9,One,Work(I000),1,Work(I010),1)
         Call Putall(Work(I010),Length_9,Irrepx,9)

         Call Getall(Work(I000),Length_10,Irrepx,110)
         Call Getall(Work(I010),Length_10,Irrepx,10)
         Call Daxpy(Length_10,One,Work(I000),1,Work(I010),1)
         Call Putall(Work(I010),Length_10,Irrepx,10)
      Else
         Length_7=IDSYMSZ(IRREPX,ISYTYP(1,7),ISYTYP(2,7))
         Length_10=IDSYMSZ(IRREPX,ISYTYP(1,10),ISYTYP(2,10))

         Maxln = Max(Length_7,Length_10)
         I000 = Ione
         I010 = I000 + Maxln
         Iend = I010 + Maxln
         If (Iend .Gt. Length) Call Insmem("add_external_t3_to_h12",
     +                                      Iend,Length)

         Call Getall(Work(I000),Length_7,Irrepx,107) 
         Call Getall(Work(I010),Length_7,Irrepx,7)
         Call Daxpy(Length_7,One,Work(I000),1,Work(I010),1)
         Call Putall(Work(I010),Length_7,Irrepx,7)

         Call Getall(Work(I000),Length_10,Irrepx,110)
         Call Getall(Work(I010),Length_10,Irrepx,10)
         Call Daxpy(Length_10,One,Work(I000),1,Work(I010),1)
         Call Putall(Work(I010),Length_10,Irrepx,10)

      Endif

      If (Uhf) then
         Length_27=IDSYMSZ(IRREPX,ISYTYP(1,27),ISYTYP(2,27))
         Length_28=IDSYMSZ(IRREPX,ISYTYP(1,28),ISYTYP(2,28))
         Length_29=IDSYMSZ(IRREPX,ISYTYP(1,29),ISYTYP(2,29))
         Length_30=IDSYMSZ(IRREPX,ISYTYP(1,30),ISYTYP(2,30))

         Maxln = Max(Length_27,Length_28,Length_29,Length_30)
         I000 = Ione
         I010 = I000 + Maxln
         Iend = I010 + Maxln
         If (Iend .Gt. Length) Call Insmem("add_external_t3_to_h12",
     +                                      Iend,Length)

         Call Getall(Work(I000),Length_27,Irrepx,127)
         Call Getall(Work(I010),Length_27,Irrepx,27)
         Call Daxpy(Length_27,One,Work(I000),1,Work(I010),1)
         Call Putall(Work(I010),Length_27,Irrepx,27)

         Call Getall(Work(I000),Length_28,Irrepx,128)
         Call Getall(Work(I010),Length_28,Irrepx,28)
         Call Daxpy(Length_28,One,Work(I000),1,Work(I010),1)
         Call Putall(Work(I010),Length_28,Irrepx,28)

         Call Getall(Work(I000),Length_29,Irrepx,129)
CSSS         write(6,"(a)") "List-129-UHF"
CSSS         write(6,"(4(1x,F15.9))") (Work(i000+i-1),i=1,Length_29)
         Call Getall(Work(I010),Length_29,Irrepx,29)
CSSS         write(6,"(a)") "List-29-UHF"
CSSS         write(6,"(4(1x,F15.9))") (Work(i010+i-1),i=1,Length_29)
         Call Daxpy(Length_29,One,Work(I000),1,Work(I010),1)
         Call Putall(Work(I010),Length_29,Irrepx,29)

         Call Getall(Work(I000),Length_30,Irrepx,130)
CSSS         write(6,"(a)") "List-130-UHF"
CSSS         write(6,"(4(1x,F15.9))") (Work(i000+i-1),i=1,Length_30)
         Call Getall(Work(I010),Length_30,Irrepx,30)
CSSS         write(6,"(a)") "List-30-UHF"
CSSS         write(6,"(4(1x,F15.9))") (Work(i010+i-1),i=1,Length_30)
         Call Daxpy(Length_30,One,Work(I000),1,Work(I010),1)
         Call Putall(Work(I010),Length_30,Irrepx,30)

      Else
CSSS         Length_27=IDSYMSZ(IRREPX,ISYTYP(1,27),ISYTYP(2,27))
         Length_30=IDSYMSZ(IRREPX,ISYTYP(1,30),ISYTYP(2,30))

         Maxln = Length_30
         I000 = Ione
         I010 = I000 + Maxln
         Iend = I010 + Maxln
         If (Iend .Gt. Length) Call Insmem("add_external_t3_to_h12",
     +                                      Iend,Length)

CSSS         Call Getall(Work(I000),Length_27,Irrepx,127)
CSSS         Call Getall(Work(I010),Length_27,Irrepx,27)
CSSS         Call Daxpy(Length_27,One,Work(I000),1,Work(I010),1)
CSSS         Call Putall(Work(I010),Length_27,Irrepx,27)
 
         Call Getall(Work(I000),Length_30,Irrepx,130)
CSSS         write(6,"(a)") "List-130-RHF"
CSSS         write(6,"(4(1x,F15.9))") (Work(i000+i-1),i=1,Length_30)
         Call Getall(Work(I010),Length_30,Irrepx,30)
CSSS         write(6,"(a)") "List-30-RHF"
CSSS         write(6,"(4(1x,F15.9))") (Work(i010+i-1),i=1,Length_30)
         Call Daxpy(Length_30,One,Work(I000),1,Work(I010),1)
         Call Putall(Work(I010),Length_30,Irrepx,30)
      Endif
      Write(6,*) 

      Return
      End
