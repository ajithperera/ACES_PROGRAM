










      Subroutine tdee_iterate_pre(work,Memleft,Iuhf,Irrepx,Iside)

      Implicit Integer(A-Z)

      Dimension Work(Memleft)
      Logical Init 



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



c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 

c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
   
C Source and target singels

      pha_length = Irpdpd(Irrepx,9)
      phb_length = Irpdpd(Irrepx,10)
      hha_length = Irpdpd(Irrepx,21)
      hhb_length = Irpdpd(Irrepx,22)
      ppa_length = Irpdpd(Irrepx,19)
      ppb_length = Irpdpd(Irrepx,20)
    
      Write(6,*) 
      Write(6,"(a,l)") " Entered tdee_iterate-pre:"
      Write(6,*) 

      If (Iside .EQ.1) Then
         Call Getlst(Work,1,1,1,1,490)
         Call checksum("list-490-AA      :",Work,Pha_length,s)
         Call Putlst(Work,1,1,1,1,390)
         If (Iuhf .EQ. 1) Then
            Call Getlst(Work,1,1,1,2,490)
            Call checksum("list-490-BB      :",Work,Phb_length,s)
            Call Putlst(Work,1,1,1,2,390)
         Endif 
      Else if (Iside .EQ. 2) Then
         Call Getlst(Work,1,1,1,1,490)
         Call checksum("list-490-AA      :",Work,Pha_length,s)
         Call Putlst(Work,1,1,1,1,392)
         If (Iuhf .EQ. 1) Then
            Call Getlst(Work,1,1,1,2,490)
            Call checksum("list-490-BB      :",Work,Phb_length,s)
            Call Putlst(Work,1,1,1,2,392)
         Endif
      Endif 

      If (Iside .EQ. 1) Then
          If (Iuhf .EQ. 0) Then
              Length = idsymsz(Irrepx,isytyp(1,446),isytyp(2,446))
              Call Getall(Work,Length,Irrepx,446)
              Call checksum("List-446;Z(Ai,Bj):",Work,Length,s)
              Call Putall(Work,Length,Irrepx,316)
           Else
              Length = idsymsz(Irrepx,isytyp(1,444),isytyp(2,444))
              Call Getall(Work,Length,Irrepx,444)
              Call checksum("List-444;Z(AB,IJ):",Work,Length,s)
              Call Putall(Work,Length,Irrepx,314)
              Length = idsymsz(Irrepx,isytyp(1,445),isytyp(2,445))
              Call Getall(Work,Length,Irrepx,445)
              Call checksum("List-445;Z(ab,ij):",Work,Length,s)
              Call Putall(Work,Length,Irrepx,315)
              Length = idsymsz(Irrepx,isytyp(1,446),isytyp(2,446))
              Call Getall(Work,Length,Irrepx,446)
              Call checksum("List-446;Z(Ai,Bj):",Work,Length,s)
              Call Putall(Work,Length,Irrepx,316)
           Endif

      Else if (Iside .EQ.2) Then
          If (Iuhf .EQ. 0) Then
              Length = idsymsz(Irrepx,isytyp(1,446),isytyp(2,446))
              Call Getall(Work,Length,Irrepx,443)
              Call checksum("List-443;Z(Ai,Bj):",Work,Length,s)
              Call Putall(Work,Length,Irrepx,326)
           Else
              Length = idsymsz(Irrepx,isytyp(1,444),isytyp(2,444))
              Call Getall(Work,Length,Irrepx,444)
              Call checksum("List-444;Z(AB,IJ):",Work,Length,s)
              Call Putall(Work,Length,Irrepx,324)
              Length = idsymsz(Irrepx,isytyp(1,445),isytyp(2,445))
              Call Getall(Work,Length,Irrepx,445)
              Call checksum("List-445;Z(ab,ij):",Work,Length,s)
              Call Putall(Work,Length,Irrepx,325)
              Length = idsymsz(Irrepx,isytyp(1,446),isytyp(2,446))
              Call Getall(Work,Length,Irrepx,446)
              Call checksum("List-446;Z(Ai,Bj):",Work,Length,s)
              Call Putall(Work,Length,Irrepx,326)
           Endif
      Endif 
 
      Return
      End
 
      
