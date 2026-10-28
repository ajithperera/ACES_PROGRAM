      Subroutine Pccd_fcmask(Work,Maxcor,List)
C
C Frozen (dropped) MOs with OO_METHOD=RFO (see runpccd/pccd_fc_setup):
C zero every element T2(Ab,Ij) of list List (RHF, A-b/I-j full
C distributions) with a frozen index: I or j among the first NFZO
C occupied MOs of its irrep, A or b among the last NFZV virtual MOs of
C its irrep. Applied to the guess T2 and to each new increment, this
C keeps the frozen pairs out of the pCCD amplitudes (and Lambda, which
C starts from T2), so that the full-space densities, orbital gradient
C and Hessian are those of the active-space pCCD energy with the
C frozen MOs in the reference. No-op without the PCCDNFZO record or
C when nothing is frozen.
C
      Implicit Double Precision(A-H,O-Z)
      Integer Dissiz
      Dimension Work(Maxcor),Nfzo(8),Nfzv(8)
#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "sympop.com"

      Call Getrec(0,"JOBARC","PCCDNFZO",Length,Junk)
      If (Length .Le. 0) Return
      Call Getrec(20,"JOBARC","PCCDNFZO",8,Nfzo)
      Call Getrec(20,"JOBARC","PCCDNFZV",8,Nfzv)
      Nfz = 0
      Do Irr = 1, Nirrep
         Nfz = Nfz + Nfzo(Irr) + Nfzv(Irr)
      Enddo
      If (Nfz .Eq. 0) Return
C Layout as in Pccd_reset: rows A-b (type 15), columns I-j (type 14).

      Do Irrep = 1, Nirrep
         Numdis = Irpdpd(Irrep,14)
         Dissiz = Irpdpd(Irrep,15)
         If (Numdis*Dissiz .Eq. 0) Go to 10
         Ncol = 0
         Nrow = 0
         Do Irb = 1, Nirrep
            Ncol = Ncol + Pop(Irb,1)*Pop(Dirprd(Irb,Irrep),1)
            Nrow = Nrow + Vrt(Irb,1)*Vrt(Dirprd(Irb,Irrep),1)
         Enddo
         If (Ncol .Ne. Numdis .Or. Nrow .Ne. Dissiz) Then
            Write(6,"(a,i4)") " @-Pccd_fcmask: layout mismatch, list",
     +                        List
            Call Errex
         Endif
         If (Numdis*Dissiz .Gt. Maxcor) Call Insmem("Pccd_fcmask",
     +                                   Numdis*Dissiz,Maxcor)
         Call Getlst(Work,1,Numdis,1,Irrep,List)
C Columns: distributions Ij, I fastest within the (irrep I, irrep j)
C block, blocks ordered by the irrep of j.
         Icol = 0
         Do Irj = 1, Nirrep
            Iri = Dirprd(Irj,Irrep)
            Do J = 1, Pop(Irj,1)
               Do I = 1, Pop(Iri,1)
                  Icol = Icol + 1
                  If (I .Le. Nfzo(Iri) .Or. J .Le. Nfzo(Irj))
     +               Call Dzero(Work((Icol-1)*Dissiz+1),Dissiz)
               Enddo
            Enddo
         Enddo
C Rows: Ab, A fastest, blocks ordered by the irrep of b.
         Irow = 0
         Do Irb = 1, Nirrep
            Ira = Dirprd(Irb,Irrep)
            Do Ib = 1, Vrt(Irb,1)
               Do Ia = 1, Vrt(Ira,1)
                  Irow = Irow + 1
                  If (Ia .Gt. Vrt(Ira,1)-Nfzv(Ira) .Or.
     +                Ib .Gt. Vrt(Irb,1)-Nfzv(Irb)) Then
                     Do K = 1, Numdis
                        Work((K-1)*Dissiz+Irow) = 0.0D0
                     Enddo
                  Endif
               Enddo
            Enddo
         Enddo
         Call Putlst(Work,1,Numdis,1,Irrep,List)
 10      Continue
      Enddo
      Return
      End
