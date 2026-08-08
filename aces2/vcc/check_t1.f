










      Subroutine Check_T1(Work,Length,Iuhf)
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


C T1-AI
      Irrepx = 1
      Length_aa =Irpdpd(Irrepx,9)
      call getlst(Work,1,1,1,1,90)
      Call checksum("T1-AA   :",Work,Length_aa)
      call output(work,1,Length_aa,1,1,Length_aa,1,1)
      If (Iuhf .Ne. 0) Then 
      Length_bb =Irpdpd(Irrepx,10)
      call getlst(Work,1,1,1,2,90)
      Call checksum("T1-BB   :",Work,Length_bb)
      Endif 
      Write(6,*) 

      Return
      End

