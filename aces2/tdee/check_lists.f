










      Subroutine Check_lists(Work,Maxcor,Irrepx,Iuhf,Iside)
      Implicit Integer(A-Z)

      Double Precision Work(Maxcor)
      Logical UHF



c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end



c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end

      Write(6,*)
      Length=idsymsz(Irrepx,isytyp(1,304),isytyp(2,304))
      call getall(Work(1),length,irrepx,304)
      call checksum("List-304",Work(1),length,s)
      If (iuhf .gt.0) Then
      Length=idsymsz(Irrepx,isytyp(1,305),isytyp(2,305))
      call getall(Work(1),length,irrepx,305)
      call checksum("List-305",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,306),isytyp(2,306))
      call getall(Work(1),length,irrepx,306)
      call checksum("List-306",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,308),isytyp(2,308))
      call getall(Work(1),length,irrepx,308)
      call checksum("List-308",Work(1),length,s)
      endif 
      Length=idsymsz(Irrepx,isytyp(1,307),isytyp(2,307))
      call getall(Work(1),length,irrepx,307)
      call checksum("List-307",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,309),isytyp(2,309))
      call getall(Work(1),length,irrepx,309)
      call checksum("List-309",Work(1),length,s)

      Length=idsymsz(Irrepx,isytyp(1,350),isytyp(2,350))
      call getall(Work(1),length,irrepx,350)
      call checksum("List-350",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,351),isytyp(2,351))
      call getall(Work(1),length,irrepx,351)
      call checksum("List-351",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,352),isytyp(2,352))
      call getall(Work(1),length,irrepx,352)
      call checksum("List-352",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,353),isytyp(2,353))
      call getall(Work(1),length,irrepx,353)
      call checksum("List-353",Work(1),length,s)

      If (iside .eq. 2) Then
      Length=idsymsz(Irrepx,isytyp(1,360),isytyp(2,360))
      call getall(Work(1),length,irrepx,360)
      call checksum("List-360",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,361),isytyp(2,361))
      call getall(Work(1),length,irrepx,361)
      call checksum("List-361",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,362),isytyp(2,362))
      call getall(Work(1),length,irrepx,362)
      call checksum("List-362",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,363),isytyp(2,363))
      call getall(Work(1),length,irrepx,363)
      Endif 
      Return 
      End
