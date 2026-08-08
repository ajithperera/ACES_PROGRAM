










      Subroutine Check_vee_lists(Work,Maxcor,Irrepx,Iuhf,Iside)
      Implicit Integer(A-Z)

      Double Precision Work(Maxcor)
      Logical UHF
      Logical David,Block_david,Lanczos



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
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end

      Common /Rootsearch/David,Block_david,Lanczos

      pha_length = Irpdpd(Irrepx,9)
      phb_length = Irpdpd(Irrepx,10)
      hha_length = Irpdpd(Irrepx,21)
      hhb_length = Irpdpd(Irrepx,22)
      ppa_length = Irpdpd(Irrepx,19)
      ppb_length = Irpdpd(Irrepx,20)

      Write(6,*) 
      If (iside .Eq. 1) then 

      Call Getlst(Work,1,1,1,3,490)
      Call checksum("list-490-AA      :",Work,Pha_length,s)
      If (Iuhf .Ne. 0) Then
         Call Getlst(Work,1,1,1,4,490)
         Call checksum("list-490-BB      :",Work,Phb_length,s)
      Endif

      Elseif (Iside .Eq. 2) Then

      If (Lanczos) Then
      Call Getlst(Work,1,1,1,3,493)
      Call checksum("list-493-AA      :",Work,Pha_length,s)
      If (Iuhf .Ne. 0) Then
         Call Getlst(Work,1,1,1,4,493)
         Call checksum("list-493-BB      :",Work,Phb_length,s)
      Endif 
      Else 
      Call Getlst(Work,1,1,1,3,490)
      Call checksum("list-490-AA      :",Work,Pha_length,s)
      If (Iuhf .Ne. 0) Then
         Call Getlst(Work,1,1,1,4,490)
         Call checksum("list-490-BB      :",Work,Phb_length,s)
      Endif 
      Endif 
      Endif 

      If (iside .eq. 1) Then 
      If (iuhf .ne. 0) Then
      Length=idsymsz(Irrepx,isytyp(1,461),isytyp(2,461))
      call getall(Work(1),length,irrepx,461)
      call checksum("List-461;Z(AB,IJ):",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,462),isytyp(2,462))
      call getall(Work(1),length,irrepx,462)
      call checksum("List-462;Z(ab,ij):",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,463),isytyp(2,463))
      call getall(Work(1),length,irrepx,463)
      call checksum("List-463;Z(Ab,Ij):",Work(1),length,s)
      Else
      Length=idsymsz(Irrepx,isytyp(1,463),isytyp(2,463))
      call getall(Work(1),length,irrepx,463)
      call checksum("List-463;Z(Ai,Bj):",Work(1),length,s)
      Endif 

      elseif (iside .eq. 2) then

      If (Lanczos) Then
      If (iuhf .ne. 0) Then
      Length=idsymsz(Irrepx,isytyp(1,467),isytyp(2,467))
      call getall(Work(1),length,irrepx,467)
      call checksum("List-467;Z(AB,IJ):",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,468),isytyp(2,468))
      call getall(Work(1),length,irrepx,468)
      call checksum("List-468;Z(ab,ij):",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,469),isytyp(2,469))
      call getall(Work(1),length,irrepx,469)
      call checksum("List-469;Z(Ab,Ij):",Work(1),length,s)
      Else
      Length=idsymsz(Irrepx,isytyp(1,469),isytyp(2,469))
      call getall(Work(1),length,irrepx,469)
      call checksum("List-469;Z(Ai,Bj):",Work(1),length,s)
      Endif 
      Else
      If (iuhf .ne. 0) Then
      Length=idsymsz(Irrepx,isytyp(1,461),isytyp(2,461))
      call getall(Work(1),length,irrepx,461)
      call checksum("List-461;Z(AB,IJ):",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,462),isytyp(2,462))
      call getall(Work(1),length,irrepx,462)
      call checksum("List-462;Z(ab,ij):",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,463),isytyp(2,463))
      call getall(Work(1),length,irrepx,463)
      call checksum("List-463;Z(Ab,Ij):",Work(1),length,s)
      Else
      Length=idsymsz(Irrepx,isytyp(1,463),isytyp(2,463))
      call getall(Work(1),length,irrepx,463)
      call checksum("List-463;Z(Ai,Bj):",Work(1),length,s)
      Endif
      Endif
      Endif 

      Return 
      End
