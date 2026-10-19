#include "flags.h"
#include "disttype.h"

      Subroutine Pccd_rfo_tguess(Work,Maxcor,Mode)
C
C OO_METHOD=RFO amplitude continuation. Away from equilibrium the pCCD
C equations have several solutions at the same orbitals and, starting
C each macro iteration from the MP2 guess, DIIS lands on one or another
C of them depending on round-off (N2 cc-pVDZ, R=2.5 A: two branches
C 3.4 mEh apart, picked at random for orbital steps of 1e-17). The
C orbital optimizer then sees a discontinuous energy and the trust
C radius collapses. Starting T2 from the amplitudes of the last
C accepted orbitals keeps the solver on one branch. Mode:
C   "LOAD ": if RFO_T2C exists, overwrite T2(AB) (list 46) with it;
C            otherwise (first macro iteration) with T2 from a TGUESS
C            file of an earlier calculation, if present and of the
C            right size (geometry continuation: copy TGUESS with the
C            orbitals, oorot LOAD). Called before the T2 resorts in
C            Pccd_rdriver.
C   "SAVE ": converged T2(AB) -> RFO_T2N (after Pccd_rdriver).
C   "PROMO": RFO_T2N -> RFO_T2C (the RFO step accepted this point) and
C            to the TGUESS file (DMPTGSS layout, RHF: T1 = 0, list 44
C            slot = 0, list 46), so that the file holds the T2 of the
C            last accepted point (the minimum at the end).
C No-op unless OO_METHOD=RFO. RHF (list 46 only).
C
      Implicit Double Precision(A-H,O-Z)
      Character*5 Mode
      Character*80 Fulnam
      Logical Exist
      Dimension Work(Maxcor)
#include "machsp.com"
#include "sym.com"
#include "sympop.com"
#include "flags2.com"

      If (Iflags2(h_IFLAGS2_oo_method) .Ne. 6) Return
      Nsize = Isymsz(Isytyp(1,46),Isytyp(2,46))
C TGUESS layout (Pccd_ddmptgss/DMPTGSS, RHF): T1, T2(A<B,I<J), T2(Ab,Ij)
      N44  = Isymsz(h_DistType_AltB,h_DistType_IltJ)
      Ltg  = Nt(1) + N44 + Nsize
      If (Nsize+Ltg .Gt. Maxcor) Call Insmem("Pccd_rfo_tguess",
     +                                       Nsize+Ltg,Maxcor)

      If (Mode .Eq. "LOAD ") Then
         Call Getrec(0,"JOBARC","RFO_T2C ",Length,Junk)
         If (Length .Le. 0) Then
            Call Gfname("TGUESS  ",Fulnam,Namlen)
            Inquire(File=Fulnam(1:Namlen),Exist=Exist,Size=Isz)
            If (.Not. Exist) Return
C Unformatted sequential: one record, 4-byte markers on both sides.
            If (Isz .Ne. 8*Ltg+8) Then
               Write(6,"(2a)") " @-RFO TGUESS file ignored: size does",
     +              " not match this calculation."
               Return
            Endif
            Call Rdtgss(Work,Ltg,"TGUESS  ",Exist)
            Call Putall(Work(Nt(1)+N44+1),Nsize,1,46)
            Write(6,"(2a)") " @-RFO T2 guess: amplitudes from the",
     +                      " TGUESS file."
            Return
         Endif
         If (Length .Ne. Nsize*Iintfp) Then
            Write(6,"(a)") " @-Pccd_rfo_tguess: RFO_T2C size mismatch."
            Call Errex
         Endif
         Call Getrec(20,"JOBARC","RFO_T2C ",Nsize*Iintfp,Work)
         Call Putall(Work,Nsize,1,46)
         Write(6,"(2a)") " @-RFO T2 guess: amplitudes of the last",
     +                   " accepted orbitals."
      Elseif (Mode .Eq. "SAVE ") Then
         Call Getall(Work,Nsize,1,46)
         Call Putrec(20,"JOBARC","RFO_T2N ",Nsize*Iintfp,Work)
      Elseif (Mode .Eq. "PROMO") Then
         Call Getrec(0,"JOBARC","RFO_T2N ",Length,Junk)
         If (Length .Le. 0) Return
         Call Getrec(20,"JOBARC","RFO_T2N ",Nsize*Iintfp,Work)
         Call Putrec(20,"JOBARC","RFO_T2C ",Nsize*Iintfp,Work)
         Call Dzero(Work(Nsize+1),Nt(1)+N44)
         Call Dcopy(Nsize,Work,1,Work(Nsize+Nt(1)+N44+1),1)
         Call Dmptgss(Work(Nsize+1),Ltg,"TGUESS  ")
      Endif
      Return
      End
