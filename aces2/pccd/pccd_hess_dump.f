










      Subroutine Pccd_hess_dump(Work,Maxcor,Tag)
C
C Debug (pccd_hess.in, see Pccd_do_orbrots): Tag = "ZERO" zeroes the
C orbital Hessian lists; any other Tag writes each of them to
C hess_<Tag>_<list>.dat: a header line (list, row type, column type,
C rows, columns) and the irrep-1 block, rows fastest. RHF, SYMMETRY=OFF
C intended (one irrep block). Tag "GA" writes the two-particle density
C lists (100-160, those that exist) the same way.
C
      Implicit Double Precision(A-H,O-Z)
      Character*4 Tag
      Character*24 Fname
      Dimension Work(Maxcor),Lists(11)
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
      Data Lists /213,234,216,217,210,211,230,231,223,225,218/

      Irrepx = 1
C Tag "GA": the two-particle density lists 100-160 instead.
C Tag "AL": all lists 1-250 that exist.
      Nl = 11
      If (Tag .Eq. "GA  ") Nl = 61
      If (Tag .Eq. "AL  ") Nl = 250
      Do Il = 1, Nl
         List = Lists(Min(Il,11))
         If (Tag .Eq. "GA  " .Or. Tag .Eq. "AL  ") Then
            List = 99 + Il
            If (Tag .Eq. "AL  ") List = Il
            If (Isytyp(1,List) .Eq. 0 .Or. Isytyp(2,List) .Eq. 0)
     +         Go to 10
         Endif
         Nrow  = Irpdpd(Irrepx,Isytyp(1,List))
         Ncol  = Irpdpd(Irrepx,Isytyp(2,List))
         Nsize = Idsymsz(Irrepx,Isytyp(1,List),Isytyp(2,List))
         If (Nsize .Gt. Maxcor) Call Insmem("Pccd_hess_dump",Nsize,
     +                                      Maxcor)
         If (Tag .Eq. "ZERO") Then
            Call Dzero(Work,Nsize)
            Call Putall(Work,Nsize,Irrepx,List)
         Else
            Call Getall(Work,Nsize,Irrepx,List)
            Write(Fname,"(a5,a,a1,i3,a4)") "hess_",Tag(1:Len_trim(Tag)),
     +                                     "_",List,".dat"
            Open(Unit=43,File=Fname,Status="Unknown")
            Write(43,*) List,Isytyp(1,List),Isytyp(2,List),Nrow,Ncol,
     +                  Nsize
            Write(43,"(4(1x,E25.17))") (Work(I),I=1,Nsize)
            Close(43)
         Endif
 10      Continue
      Enddo
      Return
      End
