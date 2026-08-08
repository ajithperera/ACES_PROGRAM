










      Subroutine Add_t3in2t2_abij(Work,Length,Iuhf)

      Implicit Double Precision (A-H,O-Z)
      Logical Uhf

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


       Write(6,"(a)") "-----Entered add_t3int2_abij-----"

      Irrepx = 1
      Uhf    = (Iuhf .NE. 0)
C Tai contributions are temporarily stored in 10,12 coulmns of list 
C 90. Add them to the current increment. 

      Lenhp_aa = Irpdpd(Irrepx,9)
      Lenhp_bb = Irpdpd(Irrepx,10)
      Maxln = Max(Lenhp_aa,Lenhp_bb)
      I000 = Ione
      I010 = I000 + Maxln
      Iend = I010 + Maxln
      If (Iend .Gt. Length)Call Insmem ("add_external_t3_to_h12",Maxln,
     +                                  Length)
      Call Getlst(Work(I000),1,1,1,5,93)
      call checksum("T1-AA   :",Work(I000),lenhp_aa)
      Call Getlst(Work(I010),1,1,1,3,90)
      Call Daxpy(Lenhp_aa,One,Work(I000),1,WOrk(I010),1)
      Call Putlst(Work(I010),1,1,1,3,90)

      If (Iuhf .Ne. 0) Then
         Call Getlst(Work(I000),1,1,1,6,93)
      call checksum("T1-BB   :",Work(I000),lenhp_bb)
         Call Getlst(Work(I010),1,1,1,4,90)
         Call Daxpy(Lenhp_bb,One,Work(I000),1,WOrk(I010),1)
         Call Putlst(Work(I010),1,1,1,4,90)
      Endif

C The <ab|ij> contributions from T3 are temporarliy stored in
C list 61-63 (T2 residual list). Add them to <ab||ij> list

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
         Call Getall(Work(I000),Length_61,Irrepx,147)
      call checksum("T2-AAAA :",Work(I000),Length_61)
         Call Getall(Work(I010),Length_61,Irrepx,61)
         Call Daxpy(Length_61,One,Work(I000),1,Work(I010),1)
         Call Putall(Work(I010),Length_61,Irrepx,61)

         Call Getall(Work(I000),Length_62,Irrepx,148)
      call checksum("T2-BBBB :",Work(I000),Length_62)
         Call Getall(Work(I010),Length_62,Irrepx,62)
         Call Daxpy(Length_62,One,Work(I000),1,Work(I010),1)
         Call Putall(Work(I010),Length_62,Irrepx,62)
      Endif

      Call Getall(Work(I000),Length_63,Irrepx,149)
      call checksum("T2-ABAB :",Work(I000),Length_63)
      Call Getall(Work(I010),Length_63,Irrepx,63)
      Call Daxpy(Length_63,One,Work(I000),1,Work(I010),1)
      Call Putall(Work(I010),Length_63,Irrepx,63)

C This block is not needed since all the redundent lists are 
C generated from ABAB (list 63).


      Return
      End
