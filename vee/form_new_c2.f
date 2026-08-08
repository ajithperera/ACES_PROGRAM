










      Subroutine form_new_c2(Work,Length,Iuhf,Irrepx,Listz0,Listd0)
      
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Length)



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
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end

C ABAB Block

      e_ab = 0.0D0
      Ioff = 1 
      Do Irrepr = 1, Nirrep
         Irrepl = Dirprd(Irrepr,Irrepx)

         Listz = Listz0 + 2
         Listd = Listd0 + 3
         Nrow  = Irpdpd(Irrepl,Isytyp(1,Listz))
         Ncol  = Irpdpd(Irrepr,Isytyp(2,Listz))
         Nsize = Nrow * Ncol

         I000  = 1
         I010  = I000  + Nsize
         I020  = I010  + Nsize
         IScr1 = I020  + Nsize
         Iscr2 = Iscr1 + Nsize
         Iend  = Iscr2 + Nsize 
         If (Iend .Ge. Length) Call Insmem("form_new_c2",
     +                                      Iend,Length)
         Call Getlst(Work(I000),1,Ncol,1,Irrepr,Listz)
         Call Getlst(Work(I010),1,Ncol,1,Irrepr,Listd)

      write(6,"(a,3i4)") "List-z,Ncol,Nrow        :",Listz,Ncol,Nrow
      call checksum("DC2_to_C2(AB):",Work(I000),nrow*ncol,s)
      call checksum("DC2_to_D2(AB):",Work(I010),nrow*ncol,s)
      Call Vecprd(Work(I000),Work(I010),Work(I020),Nsize)
      e_ab = e_ab + Ddot(Nsize,Work(I000),1,Work(I020),1)
         Call Vecprd(Work(I000),Work(I010),Work(I000),Nsize)
C         If (Iuhf .Eq. 0) Then
C            Call Spinad1(Irrepr,Pop,Nrow,Work(Ioff),Work(Iscr1),
C     +                   Work(Iscr2))
C            Ioff = Ioff + Nrow*Ncol
C         Endif
         Call Putlst(Work(I000),1,Ncol,1,Irrepr,Listz)

      Enddo
      e =  e_ab
      write(6,*)
      Write(6, "(a,1x,F15.9)") "(C1*W)^2/D     = ", e
      write(6,*)
      If (Iuhf .eq. 0) Return 

C AAAA and BBBB blocks

         e_aa = 0.0D0
         e_bb = 0.0D0
      Do Ispin = 1, 1+Iuhf
         Do Irrepr = 1, Nirrep
            Irrepl = Dirprd(Irrepr,Irrepx)

            Listz = Listz0 + Ispin - 1
            Listd = Listd0 + Ispin 
            Nrow  = Irpdpd(Irrepl,Isytyp(1,Listz))
            Ncol  = Irpdpd(Irrepr,Isytyp(2,Listz))
            Nsize = Nrow * Ncol

            I000 = 1
            I010 = I000 + Nsize 
            I020 = I010 + Nsize
            Iend = I020
            If (Iend .Ge. Length) Call Insmem("form_new_c2",
     +                                         Iend,Length)

            Call Getlst(Work(I000),1,Ncol,1,Irrepr,Listz)
            Call Getlst(Work(I010),1,Ncol,1,Irrepr,Listd)

      write(6,"(a,3i4)") "List-z,Ncol,Nrow        :", Listz,Ncol,Nrow
      if (ispin.eq.1) call checksum("DC2_to_C2(AA):",Work(I000),
     +                               nrow*ncol,s)
      if (ispin.eq.1) call checksum("DC2_to_D2(AA):",Work(I010),
     +                               nrow*ncol,s)
      if (ispin.eq.2) call checksum("DC2_to_C2(BB):",Work(I000),
     +                               nrow*ncol,s)
      if (ispin.eq.2) call checksum("DC2_to_D2(BB):",Work(I010),
     +                               nrow*ncol,s)
      if (ispin .eq. 1) then
         Call Vecprd(Work(I000),Work(I010),Work(I020),Nsize)
         e_aa = e_aa + Ddot(Nsize,Work(I000),1,Work(I020),1)
      else 
         Call Vecprd(Work(I000),Work(I010),Work(I020),Nsize)
         e_bb = e_bb + Ddot(Nsize,Work(I000),1,Work(I020),1)
      endif 
            Call Vecprd(Work(I000),Work(I010),Work(I000),Nsize)
            Call Putlst(Work(I000),1,Ncol,1,Irrepr,Listz)
         Enddo 
      Enddo 

      e = e_aa + e_bb + e_ab
      write(6,*) 
      Write(6, "(a,3(1x,F15.9))") "e_ab,e_aa,e_bb = ", 
     +                              e_ab,e_aa,e_bb
      Write(6, "(a,1x,F15.9)") "(C1*W)^2/D     = ", e
      write(6,*)

      Return
      End
