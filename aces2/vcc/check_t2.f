










      Subroutine Check_T2_VCC(Work,Length,Iuhf)
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


C T2-ABIJ
      Irrepx = 1
      Length_61=IDSYMSZ(IRREPX,ISYTYP(1,61),ISYTYP(2,61))
      If (Iuhf .NE. 0) Length_62=IDSYMSZ(IRREPX,ISYTYP(1,62),
     +                                   ISYTYP(2,62))
      Length_63=IDSYMSZ(IRREPX,ISYTYP(1,63),ISYTYP(2,63))

      Write(6,*) 
      Call Getall(Work, Length_61, Irrepx, 44)
      Call checksum("T2-AAAA :",Work,Length_61)
      If (Iuhf .ne. 0) Then
      Call Getall(Work, Length_62, Irrepx, 45)
      Call checksum("T2-BBBB :",Work,Length_62)
      Endif 
      Call Getall(Work, Length_63, Irrepx, 46)
      Call checksum("T2-ABAB :",Work,Length_63)

C T1-AI
      Length_aa =Irpdpd(Irrepx,9)
      call getlst(Work,1,1,1,1,90)
      Call checksum("T1-AA   :",Work,Length_aa)
      If (Iuhf .Ne. 0) Then 
      Length_bb =Irpdpd(Irrepx,10)
      call getlst(Work,1,1,1,2,90)
      Call checksum("T1-BB   :",Work,Length_bb)
      Endif 
      Write(6,*) 

      Return
      End

