#include "flags.h"
#include "disttype.h"

      Subroutine Scrnc_newsize(Work,Memleft,Iuhf,Irrepx,Nsizec)

      Implicit Double Precision (A-H, O-Z)

      Dimension Work(Memleft)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "parcc.com"
#include "flags.com"

      Nsizec = Irpdpd(Irrepx,9)
      If (Iuhf .EQ. 0) Then
         Nsizec = Nsizec  + Idsymsz(Irrepx,13,14)
      Else
        Nsizec = Nsizec + Irpdpd(Irrepx,10)
        Nsizec = Nsizec + Idsymsz(Irrepx,1,3)
        Nsizec = Nsizec + Idsymsz(Irrepx,2,4)
        Nsizec = Nsizec + Idsymsz(Irrepx,13,14)
      Endif

#ifdef _DEBUG_LVL0
      Write(6,"(a,i5,a,4i5)") " NEWSIZE-CHECK Nsizec=",Nsizec,
     +   " hardcoded AA/BB/AB via Idsymsz(1,3)/(2,4)/(13,14)=",
     +   Idsymsz(Irrepx,1,3),Idsymsz(Irrepx,2,4),
     +   Idsymsz(Irrepx,13,14)
      Write(6,"(a,4i5)") " NEWSIZE-CHECK via ISYTYP(list44/45/46)="
     + ,Idsymsz(Irrepx,ISYTYP(1,44),ISYTYP(2,44)),
     +  Idsymsz(Irrepx,ISYTYP(1,45),ISYTYP(2,45)),
     +  Idsymsz(Irrepx,ISYTYP(1,46),ISYTYP(2,46))
#endif

      Return
      End
