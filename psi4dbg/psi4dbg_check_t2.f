










      Subroutine Psi4dbg_check_T2(Work,Length,Iuhf)
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
      Length_44=IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
      If (Iuhf .NE. 0) Length_45=IDSYMSZ(IRREPX,ISYTYP(1,45),
     +                                   ISYTYP(2,45))
      Length_46=IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))

      Write(6,*) 
      Call Getall(Work, Length_44, Irrepx, 44)
      Call checksum("T2-AAAA:",Work,Length_44)
      If (Iuhf .ne. 0) Then
      Call Getall(Work, Length_45, Irrepx, 45)
      Call checksum("T2-BBBB:",Work,Length_45)
      Endif 
      Call Getall(Work, Length_46, Irrepx, 46)
      Call checksum("T2-ABAB:",Work,Length_46)

      Write(6,*) 

      Return
      End

