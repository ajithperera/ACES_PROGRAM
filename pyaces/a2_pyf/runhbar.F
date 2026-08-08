      Subroutine Runhbar(Iuhf)

      Implicit None
      Integer Iuhf 
      Integer Length_scf
      Integer Length_cc
      Double Precision Tmp
      Logical Scf_done
      Logical Cc_done

#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"
#include "jobarc.com"


      Call Getrec(0,"JOBARC", "SCFEVCA0",Length_scf,Tmp)
      Call Getrec(0,"JOBARC", "CC_DONE ",Length_cc,Tmp)
      Scf_done = .False. 
      Cc_done  = .False.
      If (Length_scf .Gt. 0) Scf_done = .True.
      If (Length_cc .Gt. 0)  Cc_done = .True.

      If (Scf_done) Then
         Call Init_legacies()
         Call aces_init_chemsys()
         Call Vtran(icore(i0),Icrsiz,Iuhf)
         Call aces_init_chemsys()
         Call Intproc(icore(i0),Icrsiz,Iuhf)
         Call aces_init_chemsys()
         Call Vcc(icore(i0),Icrsiz,Iuhf)
         Call Hbar(icore(i0),Icrsiz,Iuhf)
      Elseif (CC_done) Then
         Call aces_init_chemsys()
         Call Hbar(icore(i0),Icrsiz,Iuhf)
      Else
         Call Init_legacies()
         Call Ints_(Icore(i0),Icrsiz,Iuhf)
         Call V2ja(icore(i0),Icrsiz,Iuhf)
         Call Vscf(Icore(i0),Icrsiz,Iuhf)
         Call aces_init_chemsys()
         Call Vtran(icore(i0),Icrsiz,Iuhf)
         Call aces_init_chemsys()
         Call Intproc(icore(i0),Icrsiz,Iuhf)
         Call aces_init_chemsys()
         Call Vcc(icore(i0),Icrsiz,Iuhf)
         Call aces_init_chemsys()
         Call Hbar_(icore(i0),Icrsiz,Iuhf)
      Endif 

C Do not call aces_fin (use might be intending to run something
C else after external manupilations).

CSSS      Call aces_fin
CSSS      Call C_free(Icore)

      Return
      End

