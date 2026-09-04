      Subroutine Runintgrt(Iuhf)
C
C KS/DFT numerical integration. The classic driver runs xvscf and then
C xintgrt as separate programs (job_control.F), and xintgrt is the stage
C that produces the exchange-correlation contribution and the final KS
C energy -- so with no entry point for it the *.ks suite cases could not
C run here at all.
C
C intgrt.F's program body has been split into intgrt_run, the same
C wrapper-split as joda/vscf/vtran/vprops/dirmp2: everything between
C aces_init_rte/aces_ja_init and aces_ja_fin, which pyaces has already
C done and will do itself.
C
      Implicit None
      Integer Iuhf
      Integer Length
      Double Precision Tmp
      Logical Scf_done
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"
#include "jobarc.com"
C
C The grid work reads the converged orbitals, so run the SCF first if this
C process has not already done one -- same test Runcc and Rundirmp2 use.
C
      Call Getrec(0,"JOBARC","SCFEVCA0",Length,Tmp)
      Scf_done = .False.
      If (Length .Gt. 0) Scf_done = .True.
      If (.Not. Scf_done) Then
         Call Ints_(Icore(i0),Icrsiz,Iuhf)
         Call V2ja(icore(i0),Icrsiz,Iuhf)
         Call Vscf(Icore(i0),Icrsiz,Iuhf)
      End If
      Call Intgrt_run
      Return
      End
