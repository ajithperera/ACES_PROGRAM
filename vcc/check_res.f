










      Subroutine Check_Res(Work,Length,Iuhf)
      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Length)

      Logical null 
      Data onem,One /-1.0,1.0/
      Integer Dissyt 
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



      INTEGER FSDPDAN,FSDPDNA,FSDPDAA,FSDPDIN,FSDPDNI,FSDPDII,FSDPDAI,
     $   FSDPDIA
      COMMON /FSSYMPOP/ FSDPDAN(8,22),FSDPDNA(8,22),FSDPDAA(8,22),
     &                  FSDPDIN(8,22),FSDPDNI(8,22),FSDPDII(8,22),
     &                  FSDPDAI(8,22),FSDPDIA(8,22)
      Irrepx = 1

      
      DO 100 IRREP=1,NIRREP
         NUMSYT=FSDPDAA(IRREP,ISYTYP(2,63))
         IF(NUMSYT.EQ.0) GOTO 100
         DISSYT=FSDPDII(IRREP,ISYTYP(1,63))
         IF(DISSYT.EQ.0) GOTO 100
           write(*,'(A,I1,3X,A,I2,3X,A,I2,3X,A,I2,3X,A,I2)')
     $          "SPIN=",ISPIN,"IRREP=",IRREP,
     $          "NUMSYT=",NUMSYT,"DISSYT=",DISSYT
            CALL FSGET(Work,1,NUMSYT,1,IRREP,63,'IIAA')
            CALL OUTPUT(Work,1,DISSYT,1,NUMSYT,DISSYT,NUMSYT,1)
 100  CONTINUE

      Return
      End

