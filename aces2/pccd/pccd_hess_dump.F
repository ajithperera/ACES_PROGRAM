      Subroutine Pccd_hess_dump(Work,Maxcor,Tag)
C
C Debug (pccd_hess.in, see Pccd_do_orbrots): Tag = "ZERO" zeroes the
C orbital Hessian lists; any other Tag writes each of them to
C hess_<Tag>_<list>.dat: a header line (list, row type, column type,
C rows, columns) and the irrep-1 block, rows fastest. RHF, SYMMETRY=OFF
C intended (one irrep block).
C
      Implicit Double Precision(A-H,O-Z)
      Character*4 Tag
      Character*24 Fname
      Dimension Work(Maxcor),Lists(11)
#include "sympop.com"
      Data Lists /213,234,216,217,210,211,230,231,223,225,218/

      Irrepx = 1
      Do Il = 1, 11
         List  = Lists(Il)
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
      Enddo
      Return
      End
