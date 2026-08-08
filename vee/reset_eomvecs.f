










      Subroutine Reset_eomvecs(Work,Maxcor,IUhf,Irrepx,Iside,pCCD)
     
      Implicit Double Precision(A-H,O-Z)
      Integer Dissiz,Pha_length
      Logical pCCD
      Dimension Work(Maxcor)
      
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end


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



c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end


      Data Ione /1/

      Pha_length = Irpdpd(Irrepx,9)

      List_r1 = 490
      List_l1 = 490

      Print*, "Pha_length", Pha_length
      If (Iside .Eq. 1) Then
         Call Getlst(Work,1,1,1,3,List_r1)
         call checksum("R1",Work,pha_length,s)
         Call Dzero(Work,Pha_length)
         Call Putlst(Work,1,1,1,3,List_r1)
      Else if (Iside .Eq. 2) Then
         Call Getlst(Work,1,1,1,3,List_l1)
         Call Dzero(Work,Pha_length)
         Call Putlst(Work,1,1,1,3,List_l1)
      Endif 

CSSS      If (.Not. pCCD) Return 
      RETURN

      List_r2 = 463
      List_l2 = List_r2

      Nsize_ab = Idsymsz(Irrepx,Isytyp(1,List_r2),Isytyp(2,List_r2))

      Do Irrep = 1, Nirrep

         Irrepr = Irrep
         Irrepl = Dirprd(Irrepr,Irrepx) 
            
          Numdis = Irpdpd(Irrepr,14)
          Dissiz = Irpdpd(Irrepl,15)

          I000 = Ione
          I010 = I000 + Numdis*Dissiz
          Iend = I010 + Numdis*Dissiz
          If (Iend .Gt. Maxcor) Call Insmem("transform_4index_ab",
     +                                       Iend,Maxcor)
C R2(Ab,Ij)
              Call Getlst(Work(I000),1,Numdis,1,Irrepr,List_r2)

          Do I = 1, Numdis

             Ioff = (I-1)*Dissiz 
             Call pccd_reset_vv_ab(Work(I000+Ioff),Dissiz,Irrepl)
          Enddo 

C R2(Ab,Ij)->R2(Ij,Ab)

          Call Transp(Work(I000),Work(I010),Numdis,Dissiz)

          do I = 1, Dissiz

             Joff = (I-1)*Numdis 
              Call pccd_reset_oo_ab(Work(I010+Joff),Numdis,Irrepr)
          Enddo
            
C R2(Ij,Ab)->R2(Ab,Ij)

          Call Transp(Work(I010),Work(I000),Dissiz,Numdis)

          Call Putlst(Work(I000),1,Numdis,1,Irrepr,List_r2)     

      Enddo 

      Return
      End
            
